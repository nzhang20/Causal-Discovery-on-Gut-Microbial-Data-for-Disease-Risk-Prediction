'''
sparsify.py contains functions used to prune and sure screen features and edges before running causal discovery algorithms
'''

import pandas as pd
import numpy as np
import os
import ruamel.yaml


##### MICROBE-MICROBE NETWORK METHODS ######


def run_sparcc(disease, group0, group1):
    '''
    Run SparCC on pairs of microbes for each cohort. Dataframe must be in the format, microbes (rows) x samples (columns). 

    :param: disease: disease of interest ('pcos' or 't2d')
    :param: group0: name of the healthy cohort
    :param: group1: name of the diseased cohort
    '''
    # check that the SparCC folder exists 
    assert os.path.isdir('SparCC')

    sparcc_vars = ['data_input', 'n_iteractions', 'x_iteractions', 'save_corr_file', 'num_simulate_data', 'outpath', 'outfile_pvals']
    sparcc_healthy_vals = [f'data/{disease}/{group0}.T.csv', 20, 10, f'data/{disease}/sparcc_{group0}.csv', 100, f'data/{disease}/pvals/', f'data/{disease}/sparcc_{group0}_pvals_one_sided.csv']
    sparcc_diseased_vals = [f'data/{disease}/{group1}.T.csv', 20, 10, f'data/{disease}/sparcc_{group1}.csv', 100, f'data/{disease}/pvals/', f'data/{disease}/sparcc_{group1}_pvals_one_sided.csv']

    healthy_change = dict(zip(sparcc_vars, sparcc_healthy_vals))
    diseased_change = dict(zip(sparcc_vars, sparcc_diseased_vals))

    yaml = ruamel.yaml.YAML()
    yaml.preserve_quotes = True
    yaml.explicit_start = True
    # run SparCC on HC
    with open('SparCC/configuration.yml', 'r') as stream:
        content = yaml.load(stream)
        
    for key, value in healthy_change.items():
        content[key] = value

    with open('SparCC/configuration.yml', 'w') as stream:
        yaml.dump(content, stream)
        
    os.system('python SparCC/General_Execution.py')

    # run SparCC on PCOS
    with open('SparCC/configuration.yml', 'r') as stream:
        content = yaml.load(stream)
        
    for key, value in diseased_change.items():
        content[key] = value

    with open('SparCC/configuration.yml', 'w') as stream:
        yaml.dump(content, stream)
        
    os.system('python SparCC/General_Execution.py')


def run_glasso(rscript_fp):
    '''
    Runs Graphical LASSO in R.

    :param: rscript_fp: filepath to the GLASSO R file
    '''
    cmd = 'Rscript ' + rscript_fp
    os.system(cmd)


def get_sig_cor_pairs_sparcc(sparcc_fp, sparcc_pval_fp, cohort_fp):
    '''
    Returns a dataframe of the pairs of genera that had a significant correlation (p <= 0.05). 

    :param: sparcc_fp: filepath to the sparcc correlations
    :param: sparcc_pval_fp: filepath to the p-values of the sparcc correlations
    :param: cohort: filepath to the dataframe corresponding to the cohort

    :return: pandas DataFrame of the significant correlations
    '''
    sparcc = pd.read_csv(sparcc_fp)
    sparcc = sparcc.iloc[:, 1:]
    cohort = pd.read_csv(cohort_fp)
    cohort = cohort.iloc[:, 1:]
    sparcc.columns = cohort.columns

    sparcc_pval = pd.read_csv(sparcc_pval_fp)
    sparcc_pval = sparcc_pval.iloc[:, 1:]
    sparcc_pval.columns = cohort.columns

    sparcc_genus_A, sparcc_genus_B, sparcc_rho, sparcc_p = [], [], [], []
    
    for i in range(sparcc_pval.shape[1]):
        for j in range(i, sparcc_pval.shape[1]):
            if sparcc_pval.iloc[i, j] <= 0.05:
                sparcc_genus_A.append(sparcc_pval.columns[i])
                sparcc_genus_B.append(sparcc_pval.columns[j])
                sparcc_rho.append(sparcc.iloc[i, j])
                sparcc_p.append(sparcc_pval.iloc[i, j])

    sparcc_sig = pd.DataFrame(data = {'genus_A': sparcc_genus_A,
                                      'genus_B': sparcc_genus_B,
                                      'rho': sparcc_rho,
                                      'p': sparcc_p})
    return sparcc_sig


def get_sig_cor_pairs_glasso(glasso_hc_fp, glasso_pcos_fp, genus_lst):
    '''
    Returns a dataframe of the pairs of genera that had a nonzero precision after estimation of the sparse precision matrix by Graphical LASSO.
    Runs both cohorts to reduce runtime.

    :param: glasso_hc_fp: filepath to the sparse precision matrix estimated by GLASSO for the HC cohort
    :param: glasso_pcos_fp: filepath to the sparse precision matrix estimated by GLASSO for the PCOS cohort
    :param: genus_lst: list of genera corresponding to the dimensions of the sparse precision matrix

    :return: GLASSO_hc: pandas DataFrame of the pairs of genera that had a nonzero precision for the HC cohort
    :return: GLASSO_pcos: pandas DataFrame of the pairs of genera that had a nonzero precision for the PCOS cohort
    :return: keep_nodes_hc: list of the unique genera in GLASSO_hc
    :return: keep_nodes_pcos: list of the unique genera in GLASSO_pcos
    '''
    GLASSO_hc = pd.read_csv(glasso_hc_fp)
    GLASSO_pcos = pd.read_csv(glasso_pcos_fp)
    GLASSO_hc.columns = genus_lst
    GLASSO_pcos.columns = genus_lst
    
    # get microbe pairs (edges) with nonzero covariance
    keep_nodes_hc = set()
    keep_nodes_pcos = set()
    
    GLASSO_hc_genus_A, GLASSO_hc_genus_B, GLASSO_hc_prec = [], [], []
    GLASSO_pcos_genus_A, GLASSO_pcos_genus_B, GLASSO_pcos_prec = [], [], []
    
    for i in range(GLASSO_hc.shape[1]-1): # both GLASSO_hc and GLASSO_pcos have the same dimension
        for j in range(i+1, GLASSO_hc.shape[1]):
            if GLASSO_hc.iloc[i, j] != 0:
                keep_nodes_hc.add(GLASSO_hc.columns[i])
                keep_nodes_hc.add(GLASSO_hc.columns[j])
                GLASSO_hc_genus_A.append(GLASSO_hc.columns[i])
                GLASSO_hc_genus_B.append(GLASSO_hc.columns[j])
                GLASSO_hc_prec.append(GLASSO_hc.iloc[i, j])
            if GLASSO_pcos.iloc[i, j] != 0:
                keep_nodes_pcos.add(GLASSO_pcos.columns[i])
                keep_nodes_pcos.add(GLASSO_pcos.columns[j])
                GLASSO_pcos_genus_A.append(GLASSO_pcos.columns[i])
                GLASSO_pcos_genus_B.append(GLASSO_pcos.columns[j])
                GLASSO_pcos_prec.append(GLASSO_pcos.iloc[i, j])
    
    GLASSO_hc = pd.DataFrame(data = {'genus_A': GLASSO_hc_genus_A,
                                     'genus_B': GLASSO_hc_genus_B,
                                     'precision': GLASSO_hc_prec})
    
    GLASSO_pcos = pd.DataFrame(data = {'genus_A': GLASSO_pcos_genus_A,
                                       'genus_B': GLASSO_pcos_genus_B,
                                       'precision': GLASSO_pcos_prec})

    keep_nodes_hc = sorted(list(keep_nodes_hc))
    keep_nodes_pcos = sorted(list(keep_nodes_pcos))

    return GLASSO_hc, GLASSO_pcos, keep_nodes_hc, keep_nodes_pcos


def sparse_to_adj(data_sparse, genus_lst):
    '''
    Returns the adjacency matrix of the significant correlations found between microbe pairs using SparCC or Graphical LASSO. 

    :param: data_sparse: pandas DataFrame of the SparCC or Graphical LASSO results for a cohort (HC or PCOS)
    :param: genus_lst: list of genera corresponding to the dimensions of the adjacency matrix 

    :return: numpy adjacency matrix
    '''
    sparse_sig = data_sparse.copy()

    adj_mat = np.zeros((len(genus_lst), len(genus_lst)))

    for i, genus_A in enumerate(genus_lst):
        print(genus_A)
        for j, genus_B in enumerate(genus_lst):
            if any((sparse_sig['genus_A'] == genus_A) & (sparse_sig['genus_B'] == genus_B)):
                adj_mat[i, j] = 1

    return adj_mat
    

##### MICROBE-DISEASE NETWORK METHODS #####


def run_lasso(rscript_fp):
    '''
    Runs Logistic Regression w/ LASSO in R. 

    :param: rscript_fp: filepath to the LASSO R file
    '''
    cmd = 'Rscript ' + rscript_fp
    os.system(cmd)


def prune_lasso(clean, metadata, LASSO_fp):
    '''
    Returns a pruned version of the clean dataset, based on the LASSO covariates in LASSO_fp. 

    :param: clean: pandas DataFrame of the clean dataset
    :param: metadata: pandas DataFrame of the metadata, needed for its columns
    :param: LASSO_fp: filepath to the LASSO covariates txt file

    :return: pandas DataFrame of the pruned dataset
    '''
    LASSO_covariates = pd.read_csv(LASSO_fp)

    # prune
    clean_short = clean[list(metadata.columns) + list(LASSO_covariates.columns[1:])] # the first LASSO ogcol is the intercept

    print('The pruned dataset has the following dimensions: ', clean_short.shape)

    return clean_short