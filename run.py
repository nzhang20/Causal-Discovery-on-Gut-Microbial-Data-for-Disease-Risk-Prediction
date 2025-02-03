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

    disease = data_params['disease']
    if disease == 'pcos':
        non_microbes = ['group', 'region', 'study_site']
        group0 = 'hc'
        group1 = 'pcos'
        cov1 = ['Europe', 'Asia']
        cov2 = list(range(1, 15))
    if disease == 't2d':
        non_microbes = ['IRIS', 'Gender', 'Ethnicity']
        group0 = 'IS'
        group1 = 'IR'
        cov1 = ['Male', 'Female']
        cov2 = ['C', 'A', 'B', 'H']
        
    if 'data' in targets:
        # make directory for disease
        if not os.path.exists(f'data/{disease}'): os.mkdir(f'data/{disease}')
        # clean data
        if disease == 'pcos':
            data = clean_data_pcos('data/raw/pcosyang2024.xlsx')

        if disease == 't2d':
            data = clean_data_t2d('data/raw/S1_Subjects.csv', 'data/raw/S3_SampleList.csv', 'data/raw/gut_16s_abundance.txt')
            
        data.to_csv(f'data/{disease}/clean.csv')

        # filter rare OTUs
        covariates = data[non_microbes]
        gut_16s = data.drop(columns=non_microbes)
        filter_rare = filter_rare_OTUs(gut_16s, data_params['rare_OTU_threshold'])
        filter_rare = pd.concat([covariates, filter_rare], axis=1)
        filter_rare.to_csv(f'data/{disease}/filter_rare.csv')
        
        # split for easier access, transpose for SparCC input
        healthy, diseased = split_data(filter_rare, non_microbes[0])
        healthy.to_csv(f'data/{disease}/{group0}.csv')
        healthy.T.to_csv(f'data/{disease}/{group0}.T.csv')
        diseased.to_csv(f'data/{disease}/{group1}.csv')
        diseased.T.to_csv(f'data/{disease}/{group1}.T.csv')

    if 'eda' in targets:
        # make directory for disease
        if not os.path.exists(f'plots/{disease}'): os.mkdir(f'plots/{disease}')
        if not os.path.exists(f'plots/{disease}/linearity'): os.mkdir(f'plots/{disease}/linearity')
        if not os.path.exists(f'plots/{disease}/normality'): os.mkdir(f'plots/{disease}/normality')
        
        # generate plots
        filter_rare = pd.read_csv(f'data/{disease}/filter_rare.csv', index_col=0)
        indiv_per_cohort(filter_rare, disease, non_microbes[0], [group0.upper(), group1.upper()])
        indiv_per_cohort_cov(filter_rare, disease, non_microbes[0], non_microbes[1], [group0.upper(), group1.upper()], cov1)
        indiv_per_cohort_cov(filter_rare, disease, non_microbes[0], non_microbes[2], [group0.upper(), group1.upper()], cov2)

        # genera only
        gut_16s = filter_rare.drop(columns = non_microbes)
        ok_linearity = input(f'Would you like to generate {len(list(itertools.combinations(gut_16s.columns, 2)))} scatter plots to check for linearity in {len(gut_16s.columns)} variables? (This may take a while.) \n(Y/N): ')
        if ok_linearity.lower() == 'y': check_linearity(gut_16s, disease)
        check_normality(gut_16s, disease)
        plot_corr_heatmap(gut_16s, disease)

    if 'graph' in targets:
        data = pd.read_csv(f'data/{disease}/filter_rare.csv', index_col=0)
        healthy = pd.read_csv(f'data/{disease}/{group0}.csv', index_col=0)
        diseased = pd.read_csv(f'data/{disease}/{group1}.csv', index_col=0)

        check_1a = input('Would you like to generate graphs for the microbe-microbe networks? (Y/N): ')
        if check_1a.lower() == 'y':
            
            # 1a) SparCC or Graphical LASSO + our algorithm
            if data_params['sparse_microbe_microbe'].lower() == 'sparcc':
                print('--- Step 1. Running SparCC ---')
                run_sparcc(disease, group0, group1)
                print('--- Finished Step 1 ---')

                print('--- Step 2. Obtain statistically significant pairs ---')
                data_sparse_healthy = get_sig_cor_pairs_sparcc(f'data/{disease}/sparcc_{group0}.csv', f'data/{disease}/sparcc_{group0}_pvals_one_sided.csv', f'data/{disease}/{group0}.csv')
                data_sparse_diseased = get_sig_cor_pairs_sparcc(f'data/{disease}/sparcc_{group1}.csv', f'data/{disease}/sparcc_{group1}_pvals_one_sided.csv', f'data/{disease}/{group1}.csv')
                print('--- Finished Step 2 ---')

                sparse_method = 'SparCC'
                sparse_method_fp = 'sparcc'

            if data_params['sparse_microbe_microbe'].lower().replace(' ', '') == 'graphicallasso':
                print('--- Step 1. Running Graphical LASSO ---')
                run_glasso('src/GLASSO.R')
                print('--- Finished Step 1 ---')

                print('--- Step 2. Obtain statistically significant pairs ---')
                data_sparse_healthy, data_sparse_diseased, keep_nodes_healthy, keep_nodes_diseased = get_sig_cor_pairs_glasso(f'data/{disease}/glasso_{group0}.csv', f'data/{disease}/glasso_{group1}.csv', list(healthy.columns))
                print('--- Finished Step 2 ---')

                healthy = healthy[keep_nodes_healthy]
                diseased = diseased[keep_nodes_diseased]

                sparse_method = 'Graphical LASSO'
                sparse_method_fp = 'glasso'

            print(f'--- Step 3. Converting {sparse_method} results to adjacency matrices ---')
            print(f'--- {group0.upper()} ---')
            healthy_adj = sparse_to_adj(data_sparse_healthy, list(healthy.columns))
            print(f'--- {group1.upper()} ---')
            diseased_adj = sparse_to_adj(data_sparse_diseased, list(diseased.columns))
            print('--- Finished Step 3 ---')

            # make directory for graph
            if not os.path.exists(f'graphs/{disease}'): os.mkdir(f'graphs/{disease}')
                
            print('--- Step 4. Graphing adjacency matrices ---')
            graph_networkx(healthy_adj, list(healthy.columns), f'{disease}/{sparse_method_fp}_{group0}')
            graph_networkx(diseased_adj, list(diseased.columns), f'{disease}/{sparse_method_fp}_{group1}')
            print('--- Finished Step 4 ---')
            
            print('--- Step 5. Running our algorithm ---')
            print(f'--- {group0.upper()} ---')
            healthy_ouralg = run_ouralg(healthy_adj, healthy, list(healthy.columns), fisherz)
            print(f'--- {group1.upper()} ---')
            diseased_ouralg = run_ouralg(diseased_adj, diseased, list(diseased.columns), fisherz)
            print('--- Finished Step 5 ---')
    
            print('--- Step 6. Graphing resulting adjacency matrices ---')
            graph_networkx(healthy_ouralg, list(healthy.columns), f'{disease}/{sparse_method_fp}_{group0}_ouralg')
            graph_networkx(diseased_ouralg, list(diseased.columns), f'{disease}/{sparse_method_fp}_{group1}_ouralg')
            print('--- Finished Step 6 ---')

        check_1b = input('Would you like to generate graphs for the microbe-disease network? (Y/N): ')
        if check_1b.lower() == 'y':
            # make directory for graph
            if not os.path.exists(f'graphs/{disease}'): os.mkdir(f'graphs/{disease}')
                
            # 1b) Logistic LASSO + CD-NOD
            run_lasso('src/LASSO.R')
            data_loglasso = prune_lasso(data, f'data/{disease}/lasso_covariates.txt')
            run_cdnod(data_loglasso, disease, f'{disease}/cdnod')

    if 'predict' in targets:
        print('prediction TBD')

    if 'clean' in targets:
        data_files = glob.glob(f'data/{disease}/*')
        plots_files = glob.glob(f'plots/{disease}/*')
        graphs_files = glob.glob(f'graphs/{disease}/*')
        for f in data_files + plots_files + graphs_files:
            if os.path.isdir(f): 
                for f1 in glob.glob(f'{f}/*'): 
                    os.remove(f1)
                os.rmdir(f)
            else: os.remove(f)

if __name__ == '__main__':
    targets = sys.argv[1:]
    if 'all' in targets:
        main(['data', 'eda', 'graph', 'predict'])
    else:
        main(targets)