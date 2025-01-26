'''
graph.py contains functions used to run causal discovery algorithms on the clean dataset
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import pydot

from causallearn.search.ConstraintBased.PC import pc
from causallearn.search.ConstraintBased.FCI import fci
from causallearn.search.ScoreBased.GES import ges
from causallearn.search.ConstraintBased.CDNOD import cdnod
from causallearn.utils.GraphUtils import GraphUtils
from causallearn.utils.cit import CIT
from scipy.stats import spearmanr

import itertools
import networkx as nx


def run_sparcc():
    '''
    Run SparCC on pairs of microbes for each cohort. Dataframe must be in the format, microbes (rows) x samples (columns). 
    TODO: figure out how to access the sparcc git repo from here.

    '''


def get_sig_cor_pairs(sparcc_fp, sparcc_pval_fp):
    '''
    Returns a dataframe of the pairs of genera that had a significant correlation (p <= 0.05). 

    :param: sparcc_fp: filepath to the sparcc correlations
    :param: sparcc_pval_fp: filepath to the p-values of the sparcc correlations

    :return: pandas DataFrame of the significant correlations
    '''
    sparcc = pd.read_csv(sparcc_fp)
    sparcc = sparcc.iloc[:, 1:]
    sparcc.columns = pcos_hc.columns

    sparcc_pval = pd.read_csv(sparcc_pval_fp)
    sparcc_pval = sparcc_pval.iloc[:, 1:]
    sparcc_pval.columns = pcos_hc.columns

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


def run_cdnod(data, fp):
    '''
    Runs the CD-NOD causal discovery algorithm on the clean dataset to obtain the microbe-disease network. 
    Also, prints the names of the microbes that are directly linked to the 'group' (disease status) node. 

    :param: data: clean dataset (must contain the 'group' node corresponding to disease status)
    :param: fp: filepath of resulting causal graph

    :return: the causallearn object, CausalGraph
    '''
    cg = cdnod(data.iloc[:, 1:].values, data[['group']].values) # 'group' must be the first column
    pyd = GraphUtils.to_pydot(cg.G, labels=list(data.columns[1:]) + [data.columns[0]])
    pyd.write_png(f"graphs/{fp}.png")

    adj_nodes = []
    for node in cg.G.get_adjacent_nodes(cg.G.get_node(f'X{data.shape[1]}')):
        adj_nodes.append(int(node.get_name().replace('X', '')))

    print('The following genera are directly linked to the \'group\' node: ', ' '.join(list(data.columns[adj_nodes])), '.')

    return cg.G

    