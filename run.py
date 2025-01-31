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
from src.sparsify import *


def main(targets):
    with open("src/data-params.json") as fh:
            data_params = json.load(fh)
        
    if 'data' in targets:
        # clean data
        data = clean_data('data/raw/pcosyang2024.xlsx')
        data.to_csv('data/clean.csv')

        # filter rare OTUs
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

        check_1a = input('Would you like to generate graphs for the microbe-microbe networks? (Y/N): ')
        if check_1a.lower() == 'y':
            
            # 1a) SparCC or Graphical LASSO + our algorithm
            if data_params['sparse_microbe_microbe'].lower() == 'sparcc':
                print('--- Step 1. Running SparCC ---')
                run_sparcc()
                print('--- Finished Step 1 ---')

                print('--- Step 2. Obtain statistically significant pairs ---')
                data_sparse_hc = get_sig_cor_pairs_sparcc('data/sparcc_hc.csv', 'data/sparcc_hc_pvals_one_sided.csv', 'data/hc.csv')
                data_sparse_pcos = get_sig_cor_pairs_sparcc('data/sparcc_pcos.csv', 'data/sparcc_pcos_pvals_one_sided.csv', 'data/pcos.csv')
                print('--- Finished Step 2 ---')

                sparse_method = 'SparCC'
                sparse_method_fp = 'sparcc'

            if data_params['sparse_microbe_microbe'].lower().replace(' ', '') == 'graphicallasso':
                print('--- Step 1. Running Graphical LASSO ---')
                run_glasso('src/GLASSO.R')
                print('--- Finished Step 1 ---')

                print('--- Step 2. Obtain statistically significant pairs ---')
                data_sparse_hc, data_sparse_pcos, keep_nodes_hc, keep_nodes_pcos = get_sig_cor_pairs_glasso('data/glasso_hc.csv', 
                                                                                                            'data/glasso_pcos.csv', 
                                                                                                             list(hc.columns))
                print('--- Finished Step 2 ---')

                hc = hc[keep_nodes_hc]
                pcos = pcos[keep_nodes_pcos]

                sparse_method = 'Graphical LASSO'
                sparse_method_fp = 'glasso'

            print(f'--- Step 3. Converting {sparse_method} results to adjacency matrices ---')
            print('--- HC ---')
            hc_adj = sparse_to_adj(data_sparse_hc, list(hc.columns))
            print('--- PCOS ---')
            pcos_adj = sparse_to_adj(data_sparse_pcos, list(pcos.columns))
            print('--- Finished Step 3 ---')
    
            print('--- Step 4. Graphing adjacency matrices ---')
            graph_networkx(hc_adj, list(hc.columns), f'{sparse_method_fp}_hc')
            graph_networkx(pcos_adj, list(pcos.columns), f'{sparse_method_fp}_pcos')
            print('--- Finished Step 4 ---')
            
            print('--- Step 5. Running our algorithm ---')
            print('--- HC ---')
            hc_ouralg = run_ouralg(hc_adj, hc, list(hc.columns), fisherz)
            print('--- PCOS ---')
            pcos_ouralg = run_ouralg(pcos_adj, pcos, list(pcos.columns), fisherz)
            print('--- Finished Step 5 ---')
    
            print('--- Step 6. Graphing resulting adjacency matrices ---')
            graph_networkx(hc_ouralg, list(hc.columns), f'{sparse_method_fp}_hc_ouralg')
            graph_networkx(pcos_ouralg, list(pcos.columns), f'{sparse_method_fp}_pcos_ouralg')
            print('--- Finished Step 6 ---')

        check_1b = input('Would you like to generate graphs for the microbe-disease network? (Y/N): ')
        if check_1b.lower() == 'y':
            # 1b) Logistic LASSO + CD-NOD
            run_lasso('src/LASSO.R')
            data_loglasso = prune_lasso(data, 'data/lasso_covariates.txt')
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