import biom
import pandas as pd
import numpy as np
import xarray as xr
import glob
import cmdstanpy
import birdman.visualization as viz
import matplotlib.pyplot as plt
from birdman import NegativeBinomial

# cmdstanpy.install_cmdstan()

fpath = "metadata_2024Feb02.tsv"
table = biom.load_table("../data/sam/genus-table-exported/feature-table.biom")
metadata = pd.read_csv(
    fpath,
    sep="\t",
    index_col=0
)

diabetes_cat = ['T2D', 'Obesity/T2D', 'Type II Diabetes', 'Type 2 diabetes', 'Type_II_Diabetes']
short_meta = metadata[(metadata['Disease'] == 'Healthy') | (metadata['Disease'].isin(diabetes_cat))]
short_meta['T2D'] = short_meta['Disease'].apply(lambda x: 'Healthy' if x == 'Healthy' else 'T2D')
short_meta = short_meta.loc[short_meta.index.intersection(table.to_dataframe().columns)]
short_table = table.filter(short_meta.index.tolist())
prevalence = short_table.to_dataframe().clip(upper=1).sum(axis=1)
features_to_keep = prevalence[prevalence >= 5].index.tolist()
short_table_filt = short_table.filter(features_to_keep, axis='observation')

# prevalence = table.to_dataframe().clip(upper=1).sum(axis=1)
# features_to_keep = prevalence[prevalence >= 20].index.tolist()
# table_filt = table.filter(features_to_keep, axis="observation")

# short_meta = metadata.loc[table.to_dataframe().columns]
# short_meta = short_meta[(short_meta['Disease_Type'] == 'Healthy') | (short_meta['Disease_Type'] == 'Metabolic')]
# table_filt = table_filt.filter(short_meta.index.tolist())

nb = NegativeBinomial(
    table=short_table_filt,
    formula="T2D",
    metadata=short_meta,
)

print('Starting model compilation')
nb.compile_model()

print('Finished model compilation')
print('Starting model fitting')
nb.fit_model(method="vi", num_draws=100)
print('Finished model fitting')
inference = nb.to_inference()
print(inference)

# break down _beta_alr_to_clr

def _alr_to_clr(x: np.ndarray) -> np.ndarray:
    """Convert ALR coordinates to centered CLR coordinates.

    :param x: Matrix of ALR coordinates (features x draws)
    :type x: np.ndarray

    :returns: Matrix of centered CLR coordinates
    :rtype: np.ndarray
    """
    num_draws = x.shape[1]
    z = np.zeros((1, num_draws))
    x_clr = np.vstack((z, x))
    x_clr = x_clr - x_clr.mean(axis=0).reshape(1, -1)
    return x_clr


def _beta_alr_to_clr(beta: np.ndarray) -> np.ndarray:
    """Convert feature-covariate coefficients from ALR to CLR.

    :param beta: Matrix of beta ALR coordinates (n draws x p covariates x
        d features)
    :type beta: np.ndarray

    :returns: Matrix of beta CLR coordinates (n draws x p covariates x d+1
        features)
    :rtype: np.ndarray
    """
    num_draws, num_covariates, num_features = beta.shape
    beta_clr = np.zeros((num_draws, num_covariates, num_features+1))
    for i in range(num_covariates):  # TODO: vectorize
        beta_slice = beta[:, i, :].T  # features x draws
        beta_clr[:, i, :] = _alr_to_clr(beta_slice).T
    return beta_clr

print('Calculating posterior CLRs')

# break down posterior_alr_to_clr
posterior = inference.posterior
alr_params=["beta_var"]
dim_replacement={"feature_alr": "feature"}
new_labels=nb.feature_names

new_posterior = posterior.copy()
for param in alr_params:
    param_da = posterior[param]
    all_chain_alr_coords = param_da
    all_chain_clr_coords = []

    for i, chain_alr_coords in all_chain_alr_coords.groupby("chain"):
        chain_clr_coords = _beta_alr_to_clr(chain_alr_coords[0])
        all_chain_clr_coords.append(chain_clr_coords)

    all_chain_clr_coords = np.array(all_chain_clr_coords)

    new_dims = [
        dim_replacement[x]
        if x in dim_replacement else x
        for x in param_da.dims
    ]
    # Replace coords with updated labels

    new_coords = dict()
    for dim in param_da.dims:
        if dim in dim_replacement:
            new_name = dim_replacement[dim]
            new_coords[new_name] = new_labels
        else:
            new_coords[dim] = param_da.coords.get(dim).data

    new_param_da = xr.DataArray(
        all_chain_clr_coords,
        dims=new_dims,
        coords=new_coords
    )
    new_posterior[param] = new_param_da

new_posterior = new_posterior.drop_vars(dim_replacement.keys())
inference.posterior = new_posterior

covariate = inference.posterior['covariate'].to_numpy()[1]

print(covariate)

ax = viz.plot_parameter_estimates(
    inference,
    parameter="beta_var",
    coords={"covariate": covariate},
)

plt.savefig("T2D_ranks.pdf")

param_means = inference.posterior['beta_var'].sel(
    **{'covariate': covariate}
    ).mean(['chain', 'draw'])
param_stds = inference.posterior['beta_var'].sel(
    **{'covariate': covariate}
    ).std(['chain', 'draw'])
sort_indices = param_means.argsort().data
param_means = param_means.data[sort_indices]
param_stds = param_stds.data[sort_indices]
param_labels = inference.posterior['feature'].data[sort_indices]

differential_results = pd.DataFrame({'feature': param_labels, 'mean clr': param_means, 'std clr': param_stds})
differential_results.to_csv('T2D_differential_results.csv')

