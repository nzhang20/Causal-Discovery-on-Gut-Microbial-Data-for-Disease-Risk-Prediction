#!/usr/bin/env python

import sys
import json
import pandas as pd
import numpy as np
import itertools
import glob, os
from biom.table import Table

from src.etl import *
from src.eda import *
from src.sparsify import *
from src.graph import *
from src.birdman import *
from src.predict import *


def main(targets):

    ##### SETUP #####
    config_file = targets[0]
    with open(config_file) as fh:
        data_params = json.load(fh)

    # load parameters
    otu_table_fp = data_params['otu_table_fp']
    metadata_fp = data_params['metadata_fp']
    disease = data_params['name']
    disease_col = data_params['disease_column']
    group0, group1 = data_params['cohort_names']
    covariates_map = data_params['covariates_map']
    rare_otu_threshold = data_params['rare_otu_threshold']
    transformation = data_params['transformation']

    # setup data dir
    if not os.path.exists(f'data/{disease}'): os.mkdir(f'data/{disease}')

    # check that otu_table and metadata exist
    otu_table, metadata = [], []
    if (not os.path.exists(otu_table_fp)) | (not os.path.exists(metadata_fp)) :
        if otu_table_fp == 'data/t2d/otu_table.csv':
            print('Generating otu table and metadata for T2D.')
            otu_table, metadata = clean_data_t2d('data/raw/S1_Subjects.csv', 'data/raw/S3_SampleList.csv', 'data/raw/gut_16s_abundance.txt')
            covariates_map = {
                'Gender': ['Male', 'Female'],
                'Ethnicity': ['C', 'A', 'B', 'H']
            }
        elif otu_table_fp == 'data/pcos/otu_table.csv':
            print('Generating otu table and metadata for PCOS.')
            otu_table, metadata = clean_data_pcos('data/raw/pcosyang2024.xlsx')
            covariates_map = {
                'region': ['Europe', 'Asia'],
                'study_site': list(range(1, 15))
            }
        else:
            raise Exception('Please provide the otu table and metadata table!')

        otu_table.to_csv(otu_table_fp)
        metadata.to_csv(metadata_fp)

    
    ##### EDA #####
    if 'eda' in targets:
        # make directory for disease
        if not os.path.exists(f'plots/{disease}'): os.mkdir(f'plots/{disease}')
        if not os.path.exists(f'plots/{disease}/linearity'): os.mkdir(f'plots/{disease}/linearity')
        if not os.path.exists(f'plots/{disease}/linearity/{transformation}'): os.mkdir(f'plots/{disease}/linearity/{transformation}')
        if not os.path.exists(f'plots/{disease}/normality'): os.mkdir(f'plots/{disease}/normality')
        if not os.path.exists(f'plots/{disease}/normality/{transformation}'): os.mkdir(f'plots/{disease}/normality/{transformation}')

        # bar charts from metadata
        metadata = pd.read_csv(metadata_fp, index_col=0)
        covariates = [x for x in metadata.columns if x != disease_col]
        indiv_per_cohort(metadata, disease, disease_col, [group0.upper(), group1.upper()])
        for cov in covariates:
            indiv_per_cohort_cov(metadata, disease, disease_col, cov, [group0.upper(), group1.upper()], covariates_map[cov])

        # genera only
        otu_table = pd.read_csv(otu_table_fp, index_col=0)
        filtered_otu_table = filter_rare_otus(otu_table, rare_otu_threshold)
        total_scatter_plots = len(list(itertools.combinations(filtered_otu_table.columns, 2)))
        ask_linearity = f'Would you like to generate {total_scatter_plots} scatter plots to check for linearity in {filtered_otu_table.shape[1]} variables?'
        if total_scatter_plots > 100: ask_linearity += ' (This may take a while.)'
        ask_linearity += '\n(Y/N): '
        if input(ask_linearity).lower() == 'y': check_linearity(filtered_otu_table, disease, transformation)

        ask_normality = f'Would you like to generate {filtered_otu_table.shape[1]} qqplots to check for normality?'
        if filtered_otu_table.shape[1] > 100: ask_normality += ' (This may take a while.)'
        ask_normality += '\n(Y/N): '
        if input(ask_normality).lower() == 'y': check_normality(filtered_otu_table, disease, transformation)
        
        plot_corr_heatmap(filtered_otu_table, disease)

    
    ##### GRAPH MICROBE-MICROBE AND MICROBE-DISEASE NETWORKS #####
    if 'graph' in targets:
        otu_table = pd.read_csv(otu_table_fp, index_col=0)
        filtered_otu_table = filter_rare_otus(otu_table, rare_otu_threshold)

        if transformation == 'clr': filtered_otu_table = clr(filtered_otu_table)
            
        filtered_otu_table.to_csv(f'data/{disease}/filtered_otu_table_{transformation}.csv')
        metadata = pd.read_csv(metadata_fp, index_col=0)
        
        merged = pd.concat([metadata, filtered_otu_table], axis=1)
        healthy = merged[merged[disease_col] == 0] #.drop(columns=[disease_col])
        diseased = merged[merged[disease_col] == 1] #.drop(columns=[disease_col])

        check_1a = input('Would you like to generate graphs for the microbe-microbe networks? (Y/N): ')
        if check_1a.lower() == 'y':
            
            # 1a) SparCC or Graphical LASSO + our algorithm
            if data_params['sparse_microbe_microbe'].lower() == 'sparcc':
                # Set up files
                healthy.T.to_csv(f'data/{disease}/{group0}.T.csv')
                diseased.T.to_csv(f'data/{disease}/{group1}.T.csv')
                
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
                healthy_otu = healthy.drop(columns=list(metadata.columns))
                healthy_otu.to_csv(f'data/{disease}/{group0}_{transformation}.csv')
                diseased_otu = diseased.drop(columns=list(metadata.columns))
                diseased_otu.to_csv(f'data/{disease}/{group1}_{transformation}.csv')
                run_glasso('src/GLASSO.R', config_file)
                print('--- Finished Step 1 ---')

                print('--- Step 2. Obtain statistically significant pairs ---')
                data_sparse_healthy, data_sparse_diseased, keep_nodes_healthy, keep_nodes_diseased = get_sig_cor_pairs_glasso(
                    f'data/{disease}/glasso_{group0}_{transformation}.csv', 
                    f'data/{disease}/glasso_{group1}_{transformation}.csv', 
                    list(healthy_otu.columns)
                )
                print('--- Finished Step 2 ---')

                healthy = healthy[keep_nodes_healthy]
                diseased = diseased[keep_nodes_diseased]

                sparse_method = 'Graphical LASSO'
                sparse_method_fp = 'glasso'

            # make directory for graph
            if not os.path.exists(f'graphs/{disease}'): os.mkdir(f'graphs/{disease}')

            print('--- Step 3. Running PC with depth 2 starting with remaining edges ---') 
            print(f'--- {group0.upper()} ---')
            run_pc_depth2(healthy, data_sparse_healthy, 0.01, f'{disease}/{sparse_method_fp}_{group0}_{transformation}')
            print(f'--- {group1.upper()} ---')
            run_pc_depth2(diseased, data_sparse_diseased, 0.01, f'{disease}/{sparse_method_fp}_{group1}_{transformation}')
            print('--- Finished Step 3 ---')

            # print(f'--- Step 3. Converting {sparse_method} results to adjacency matrices ---')
            # print(f'--- {group0.upper()} ---')
            # healthy_adj = sparse_to_adj(data_sparse_healthy, list(healthy.columns))
            # print(f'--- {group1.upper()} ---')
            # diseased_adj = sparse_to_adj(data_sparse_diseased, list(diseased.columns))
            # print('--- Finished Step 3 ---')

            # # make directory for graph
            # if not os.path.exists(f'graphs/{disease}'): os.mkdir(f'graphs/{disease}')
                
            # print('--- Step 4. Graphing adjacency matrices ---')
            # graph_networkx(healthy_adj, list(healthy.columns), f'{disease}/{sparse_method_fp}_{group0}')
            # graph_networkx(diseased_adj, list(diseased.columns), f'{disease}/{sparse_method_fp}_{group1}')
            # print('--- Finished Step 4 ---')
            
            # print('--- Step 5. Running our algorithm ---')
            # print(f'--- {group0.upper()} ---')
            # healthy_ouralg = run_ouralg(healthy_adj, healthy, list(healthy.columns), fisherz)
            # print(f'--- {group1.upper()} ---')
            # diseased_ouralg = run_ouralg(diseased_adj, diseased, list(diseased.columns), fisherz)
            # print('--- Finished Step 5 ---')
    
            # print('--- Step 6. Graphing resulting adjacency matrices ---')
            # graph_networkx(healthy_ouralg, list(healthy.columns), f'{disease}/{sparse_method_fp}_{group0}_ouralg')
            # graph_networkx(diseased_ouralg, list(diseased.columns), f'{disease}/{sparse_method_fp}_{group1}_ouralg')
            # print('--- Finished Step 6 ---')

        check_1b = input('Would you like to generate graphs for the microbe-disease network? (Y/N): ')
        if check_1b.lower() == 'y':
            # make directory for graph
            if not os.path.exists(f'graphs/{disease}'): os.mkdir(f'graphs/{disease}')
                
            # 1b) Logistic LASSO + CD-NOD
            run_lasso('src/LASSO.R', config_file)
            data_loglasso = prune_lasso(merged, metadata, f'data/{disease}/lasso_covariates_{transformation}.txt')
            print(data_loglasso.shape)
            print(data_loglasso.columns)
            run_cdnod(data_loglasso, disease, f'{disease}/cdnod_{transformation}')

    
    ##### BIRDMAN #####
    if 'birdman' in targets:
        # otu_table must have otus as the index, samples on the columns and in biom Table format
        otu_table = pd.read_csv(otu_table_fp, index_col=0).T
        otu_table = Table(otu_table.to_numpy(), observation_ids=list(otu_table.index), sample_ids=list(otu_table.columns))
        
        metadata = pd.read_csv(metadata_fp, index_col=0)
        
        print('Starting BIRDMAn')
        run_birdman(otu_table, metadata, ' + '.join(list(metadata.columns)), f'{disease}/birdman')
        print('Finished BIRDMAn')

    
    ##### PREDICT VIA VAE #####
    if 'predict' in targets:
        print('prediction TBD')

    
    ##### CLEAN #####
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
    config_file = sys.argv[1]
    targets = sys.argv[2:]

    # check config_file is json
    if config_file.endswith('.json'):
        if 'all' in targets:
            main([config_file] + ['eda', 'graph', 'birdman', 'predict'])
        else:
            main([config_file] + targets)
    else:
        raise Exception('Please provide a json file! e.g. `python config.json all`')