import biom
import pandas as pd
import numpy as np
import xarray as xr
import glob
import cmdstanpy
import birdman.visualization as viz
import matplotlib.pyplot as plt
from birdman import NegativeBinomial
from typing import Sequence

# check if cmdstan is installed
# cmdstanpy.install_cmdstan()

# the following three functions are copied from https://github.com/biocore/BIRDMAn/tree/main
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


def posterior_alr_to_clr(
    posterior: xr.Dataset,
    alr_params: list,
    dim_replacement: dict,
    new_labels: Sequence
) -> xr.Dataset:
    """Convert posterior draws from ALR to CLR.

    :param posterior: Posterior draws for fitted parameters
    :type posterior: xr.Dataset

    :param alr_params: List of parameters to convert from ALR to CLR
    :type alr_params: list

    :param dim_replacement: Dictionary of updated posterior dimension names
        e.g. {"feature_alr": "feature"}
    :type dim_replacement: dict

    :param new_labels: Coordinates to assign to CLR posterior draws
    :type new_labels: Sequence
    """
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
    return new_posterior

    
def run_birdman(otu_table, metadata, formula, result_fp):
    '''
    Runs the BIRDMAn analysis on the provided otu table (OTU abundances) and metadata file.

    :param: otu_table: a biom.table.Table object of the otu table
    :param: metadata: pandas DataFrame of the metadata
    '''
    filter_meta = metadata.loc[metadata.index.intersection(otu_table.to_dataframe().columns)]
    filter_otu = otu_table.filter(filter_meta.index.tolist(), inplace=False)

    prevalence = filter_otu.to_dataframe().clip(upper=1).sum(axis=1)
    features_to_keep = prevalence[prevalence >= 5].index.tolist()
    filter_otu = filter_otu.filter(features_to_keep, axis='observation', inplace=False)
    
    nb = NegativeBinomial(
        table = filter_otu,
        formula = formula,
        metadata = filter_meta
    )

    nb.compile_model()
    nb.fit_model(method='vi', num_draws=500)

    inference = nb.to_inference()
    inference.posterior = posterior_alr_to_clr(
        posterior = inference.posterior,
        alr_params = ['beta_var'],
        dim_replacement = {'feature_alr': 'feature'},
        new_labels = nb.feature_names
    )

    covariate = inference.posterior['covariate'].to_numpy()[1]

    ax = viz.plot_parameter_estimates(
        inference,
        parameter = 'beta_var',
        coords = {'covariate': covariate}
    )

    plt.savefig(f'graphs/{result_fp}_ranks.pdf')

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
    differential_results.to_csv(f'graphs/{result_fp}.csv')