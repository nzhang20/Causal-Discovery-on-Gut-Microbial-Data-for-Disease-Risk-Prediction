#!/usr/bin/env python

import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pylab as py
import itertools
import glob, os
from causallearn.search.ConstraintBased.CDNOD import cdnod
from causallearn.utils.GraphUtils import GraphUtils
from causallearn.utils.PCUtils.BackgroundKnowledge import BackgroundKnowledge
from causallearn.graph.GraphNode import GraphNode
import pydot
from IPython.display import Image, display

from src.etl import *
from src.eda import *
from src.graph import *


def main(targets):
    if 'data' in targets:
        # clean data
        data = clean_data('data/raw/pcosyang2024.xlsx')
        data.to_csv('data/clean.csv')

        # filter rare OTUs
        with open("config/data-params.json") as fh:
            data_params = json.load(fh)

        non_microbes = ['group', 'region', 'study_site']
        covariates = data[non_microbes]
        gut_16s = data.drop(columns=non_microbes)
        filter_rare = filter_rare_OTUs(gut_16s, data_params['rare_OTU_threshold'])
        filter_rare = pd.concat([covariates, filter_rare], axis=1)
        filter_rare.to_csv('data/filter_rare.csv')
        
        # split for easier access, transpose for SparCC input
        hc, pcos = split_data(filter_rare)
        hc.to_csv('data/hc.csv')
        hc.T.to_csv('data/hc.T.csv')
        pcos.to_csv('data/pcos.csv')
        pcos.T.to_csv('data/pcos.T.csv')

    if 'eda' in targets:
        # generate plots
        print('hi')

    if 'graph' in targets:
        data = pd.read_csv('data/filter_rare.csv', index_col=0)
        hc = pd.read_csv('data/hc.csv', index_col=0)
        pcos = pd.read_csv('data/pcos.csv', index_col=0)

        # 1a) SparCC + our algorithm
        # try Graphical LASSO
        # run SparCC on pruned (by Graphical LASSO) or full dataset
        # run our algorithm

        # print('--- Step 1. Running SparCC ---')
        # run_sparcc()
        # print('--- Finished Step 1 ---')

        # print('--- Step 2. Obtain statistically significant pairs ---')
        # data_sparcc_hc = get_sig_cor_pairs('src/hc_sparcc.csv', 'src/hc_pvals_one_sided.csv', 'data/hc.csv')
        # data_sparcc_pcos = get_sig_cor_pairs('src/pcos_sparcc.csv', 'src/pcos_pvals_one_sided.csv', 'data/pcos.csv')
        # print('--- Finished Step 2 ---')

        # print('--- Step 3. Converting SparCC results to adjacency matrices ---')
        # print('--- HC ---')
        # sparcc_hc_adj = sparcc_to_adj(data_sparcc_hc, list(hc.columns))
        # print('--- PCOS ---')
        # sparcc_pcos_adj = sparcc_to_adj(data_sparcc_pcos, list(pcos.columns))
        # print('--- Finished Step 3 ---')

        # print('--- Step 4. Graphing adjacency matrices ---')
        # graph_networkx(sparcc_hc_adj, list(hc.columns), 'sparcc_hc')
        # graph_networkx(sparcc_pcos_adj, list(pcos.columns), 'sparcc_pcos')
        # print('--- Finished Step 4 ---')
        
        # print('--- Step 5. Running our algorithm ---')
        # print('--- HC ---')
        # sparcc_hc_ouralg = run_ouralg(sparcc_hc_adj, hc, list(hc.columns), fisherz)
        # print('--- PCOS ---')
        # sparcc_pcos_ouralg = run_ouralg(sparcc_pcos_adj, pcos, list(pcos.columns), fisherz)
        # print('--- Finished Step 5 ---')

        # print('--- Step 6. Graphing resulting adjacency matrices ---')
        # graph_networkx(sparcc_hc_ouralg, list(hc.columns), 'sparcc_hc_ouralg')
        # graph_networkx(sparcc_pcos_ouralg, list(pcos.columns), 'sparcc_pcos_ouralg')
        # print('--- Finished Step 6 ---')

        # 1b) Logistic LASSO + CD-NOD
        run_lasso('src/LASSO.R')
        data_loglasso = prune_data(data, 'src/LASSO_covariates.txt')
        run_cdnod(data_loglasso, 'cdnod')

    if 'predict' in targets:
        print('prediction TBD')

    if 'clean' in targets:
        for f in glob.glob('data/*.csv'):
            os.remove(f)
        for f in glob.glob('plots/*.png'):
            os.remove(f)
        for f in glob.glob('graphs/*.png'):
            os.remove(f)

if __name__ == '__main__':
    targets = sys.argv[1:]
    if 'all' in targets:
        main(['data', 'eda', 'graph', 'predict'])
    else:
        main(targets)