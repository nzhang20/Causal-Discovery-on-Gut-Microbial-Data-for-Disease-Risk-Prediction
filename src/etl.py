'''
etl.py contains functions used to clean the raw dataframes and transforming data into adjacency matrices
'''

import pandas as pd
import numpy as np
import os


def clean_data(raw_fp):
    '''
    Returns a clean dataset.

    :param: raw_fp: filepath to the raw dataset

    :return: pandas DataFrame of the clean dataset
    '''
    pcosyang2024 = pd.read_excel(raw_fp, engine='openpyxl')
    pcosyang2024 = pcosyang2024.T
    pcosyang2024.columns = pcosyang2024.iloc[0, :]
    pcosyang2024 = pcosyang2024.iloc[1:, :]

    unclass_col = ['uclassified', 'unclassified', 'unidentified']
    pcosyang2024['Unclassified'] = pcosyang2024[unclass_col[0]] + pcosyang2024[unclass_col[1]] + pcosyang2024[unclass_col[2]]
    pcosyang2024 = pcosyang2024.drop(columns = unclass_col)

    pcosyang2024['region'] = pcosyang2024['region'].map({'Europe': 0, 'Asia': 1})
    pcosyang2024['group'] = pcosyang2024['group'].map({'HC': 0, 'PCOS': 1})
    pcosyang2024 = pcosyang2024.drop(columns = ['T'])

    for col in pcosyang2024.columns[2:]:
        pcosyang2024[col] = pcosyang2024[col].astype(float)

    return pcosyang2024


def run_lasso(rscript_fp):
    '''
    Runs Logistic Regression w/ LASSO in R. 

    :param: rscript_fp: filepath to the LASSO R file
    '''
    cmd = "Rscript " + rscript_fp
    os.system(cmd)


def prune_data(clean, LASSO_fp):
    '''
    Returns a pruned version of the clean dataset, based on the LASSO covariates in LASSO_fp. 

    :param: clean: clean dataset
    :param: LASSO_fp: filepath to the LASSO covariates txt file

    :return: pandas DataFrame of the pruned dataset
    '''
    LASSO_covariates = pd.read_csv(LASSO_fp)

    # manually fix the '.' to '-' conversion done when reading table into R
    LASSO_ogcol = []
    for col in LASSO_covariates.columns:
        if '.' in col: # having fully checked all original columns, this condition suffices
            LASSO_ogcol.append(col.replace('.', '-'))
        else:
            LASSO_ogcol.append(col)

    # prune
    pcos_clean = clean.copy()
    pcos_short = pcos_clean[['group'] + LASSO_ogcol[1:]] # the first LASSO ogcol is the intercept

    print('The pruned dataset has the following dimensions: ', pcos_short.shape)

    return pcos_short


def sparcc_to_adj(sparcc_fp, genus_lst):
    '''
    Returns the adjacency matrix of the significant correlations found between microbe pairs using SparCC. 

    :param: sparcc_fp: filepath to the SparCC results for a cohort (HC or PCOS)
    :param: genus_lst: list of genera corresponding to the dimensions of the adjacency matrix 

    :return: numpy adjacency matrix
    '''
    sparcc_sig = pd.read_csv(sparcc_fp)

    adj_mat = np.zeros((len(genus_lst), len(genus_lst)))

    for i, genus_A in enumerate(genus_lst):
        for j, genus_B in enumerate(genus_lst):
            if any((sparcc_sig['genus_A'] == genus_A) & (sparcc_sig['genus_B'] == genus_B)):
                adj_mat[i, j] = 1

    return adj_mat