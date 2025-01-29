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

    # study assigned by order presented in paper and self-matching sample size numbers & region
    study_site = {
        1: [0, 19+24],
        2: [19+24, 19+24+48+73],
        3: [19+24+48+73, 19+24+48+73+12+14],
        4: [19+24+48+73+12+14, 19+24+48+73+12+14+12],
        5: [19+24+48+73+12+14+12, 19+24+48+73+12+14+12+20+20],
        6: [19+24+48+73+12+14+12+20+20, 19+24+48+73+12+14+12+20+20+131+68],
        7: [19+24+48+73+12+14+12+20+20+131+68, 19+24+48+73+12+14+12+20+20+131+68+41+47],
        8: [19+24+48+73+12+14+12+20+20+131+68+41+47, 19+24+48+73+12+14+12+20+20+131+68+41+47+24+24],
        9: [19+24+48+73+12+14+12+20+20+131+68+41+47+24+24, 19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98],
        10: [19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98, 19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98+20+20],
        11: [19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98+20+20, 19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98+20+20+37+45],
        12: [19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98+20+20+37+45, 19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98+20+20+37+45+15+18],
        13: [19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98+20+20+37+45+15+18, 19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98+20+20+37+45+15+18+15+33],
        14: [19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98+20+20+37+45+15+18+15+33, 19+24+48+73+12+14+12+20+20+131+68+41+47+24+24+36+98+20+20+37+45+15+18+15+33+17+17]
    }

    pcosyang2024['study_site'] = np.full(pcosyang2024.shape[0], 0)
    for i in range(1, 15):
        pcosyang2024.iloc[study_site[i][0]:study_site[i][1], -1] = i

    return pcosyang2024


def filter_rare_OTUs(data, k=1):
    '''
    Filters the dataset for rare OTUs based on a chosen threshold. 

    :param: data: clean dataset, only containing abundance values
    :param: k: rare OTU abundance threshold (in percentage), default = 1%

    :return: a pandas DataFrame
    '''
    rare_OTUs = data.columns[(data <= k).sum() == data.shape[0]]
    filter_rare = data.drop(columns=rare_OTUs)
    return filter_rare


def split_data(df):
    '''
    Returns two dataframes where each corresponds to each cohort: healthy controls (HC) or PCOS patients (PCOS).

    :param: df: clean & filtered dataset, located in `data/filter_rare.csv`

    :return: a pandas DataFrame for HC, a pandas DataFrame for PCOS
    '''
    hc = df[df['group'] == 0].drop(columns = ['group'])
    pcos = df[df['group'] == 1].drop(columns = ['group'])
    return hc, pcos


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


def sparcc_to_adj(data_sparcc, genus_lst):
    '''
    Returns the adjacency matrix of the significant correlations found between microbe pairs using SparCC. 

    :param: sparcc_fp: pandas DataFrame of the SparCC results for a cohort (HC or PCOS)
    :param: genus_lst: list of genera corresponding to the dimensions of the adjacency matrix 

    :return: numpy adjacency matrix
    '''
    sparcc_sig = data_sparcc.copy()

    adj_mat = np.zeros((len(genus_lst), len(genus_lst)))

    for i, genus_A in enumerate(genus_lst):
        print(genus_A)
        for j, genus_B in enumerate(genus_lst):
            if any((sparcc_sig['genus_A'] == genus_A) & (sparcc_sig['genus_B'] == genus_B)):
                adj_mat[i, j] = 1

    return adj_mat