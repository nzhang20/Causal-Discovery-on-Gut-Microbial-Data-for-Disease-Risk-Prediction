'''
graph.py contains functions used to run causal discovery algorithms on the clean dataset
'''

import pandas as pd
import numpy as np
from numpy import ndarray
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import pydot
import itertools
import networkx as nx

from itertools import combinations
from typing import List
from tqdm.auto import tqdm
from causallearn.search.ConstraintBased.CDNOD import cdnod
from causallearn.graph.GraphClass import CausalGraph
from causallearn.utils.PCUtils import SkeletonDiscovery, UCSepset, Meek
from causallearn.utils.PCUtils.BackgroundKnowledge import BackgroundKnowledge
from causallearn.utils.PCUtils.BackgroundKnowledgeOrientUtils import orient_by_background_knowledge
from causallearn.utils.PCUtils.Helper import append_value
from causallearn.utils.GraphUtils import GraphUtils
from causallearn.utils.cit import CIT
from scipy.stats import spearmanr


def run_cdnod(data, disease, fp):
    '''
    Runs the CD-NOD causal discovery algorithm on the clean dataset to obtain the microbe-disease network. 
    Also, prints the names of the microbes that are directly linked to the cohort (disease status) node. 

    :param: data: pandas DataFrame of the LASSO-pruned features, must contain the three covariates in the first three columns
    :param: disease: disease of interest
    :param: fp: filepath of resulting causal graph

    :return: the causallearn object, CausalGraph
    '''
    # implement CD-NOD with two c_indx variables
    data_aug = data.to_numpy()
    indep_test = CIT(data_aug, 'fisherz')
    cg_1 = SkeletonDiscovery.skeleton_discovery(data_aug, 0.05, indep_test, stable=True)

    # for T2D, remove edges between 'Gender' and 'Ethnicity'
    if 'Gender' in data.columns:
        nodes = cg_1.G.get_nodes()
        bk = BackgroundKnowledge() \
            .add_forbidden_by_node(nodes[1], nodes[2]) \
            .add_forbidden_by_node(nodes[2], nodes[1])

        cg_1 = SkeletonDiscovery.skeleton_discovery(data_aug, 0.05, indep_test, stable=True, background_knowledge=bk)

    # c1 = first covariate
    c1_indx_id = 1
    for i in cg_1.G.get_adjacent_nodes(cg_1.G.nodes[c1_indx_id]):
        cg_1.G.add_directed_edge(cg_1.G.nodes[c1_indx_id], i)

    # c2 = second covariate
    c2_indx_id = 2
    for i in cg_1.G.get_adjacent_nodes(cg_1.G.nodes[c2_indx_id]):
        cg_1.G.add_directed_edge(cg_1.G.nodes[c2_indx_id], i)

    cg_2 = UCSepset.uc_sepset(cg_1, 2)
    cg = Meek.meek(cg_2)
    
    pyd = GraphUtils.to_pydot(cg.G, labels=list(data.columns))
    pyd.write_png(f"graphs/{fp}.png")

    # due to limited visibility on the graph, print the microbes that are directly linked to the cohort and covariate nodes
    adj_nodes = []
    for node in cg.G.get_adjacent_nodes(cg.G.get_node('X1')):
        adj_nodes.append(int(node.get_name().replace('X', '')))

    data_aug_col = np.array(['placeholder'] + list(data.columns))

    print(f'The following genera are directly linked to the \'{data.columns[0]}\' node: \n', '\n'.join(list(data_aug_col[adj_nodes])))

    adj_nodes = []
    for node in cg.G.get_adjacent_nodes(cg.G.get_node('X2')):
        adj_nodes.append(int(node.get_name().replace('X', '')))

    print(f'The following genera are directly linked to the \'{data.columns[1]}\' node: \n', '\n'.join(list(data_aug_col[adj_nodes])))

    adj_nodes = []
    for node in cg.G.get_adjacent_nodes(cg.G.get_node('X3')):
        adj_nodes.append(int(node.get_name().replace('X', '')))

    print(f'The following genera are directly linked to the \'{data.columns[2]}\' node: \n', '\n'.join(list(data_aug_col[adj_nodes])))

    return cg.G


def skeleton_discovery_depth2(data, data_sparse, alpha, indep_test, stable = True,
                             background_knowledge = None,
                             verbose = False,
                             show_progress = True,
                             node_names = None):
    
    assert type(data) == np.ndarray
    assert 0 < alpha < 1

    no_of_var = data.shape[1]
    cg = CausalGraph(no_of_var, node_names)
    cg.set_ind_test(indep_test)

    depth = -1
    
    pbar = tqdm(total=no_of_var) if show_progress else None
    # while cg.max_degree() - 1 > depth:
    while depth < 2:
        depth += 1
        edge_removal = []
        if show_progress:
            pbar.reset()
        for x in range(no_of_var):
            if show_progress:
                pbar.update()
            if show_progress:
                pbar.set_description(f'Depth={depth}, working on node {x}')
            Neigh_x = cg.neighbors(x)
            if len(Neigh_x) < depth - 1:
                continue
            for y in Neigh_x:
                knowledge_ban_edge = False
                sepsets = set()
                if background_knowledge is not None and (
                        background_knowledge.is_forbidden(cg.G.nodes[x], cg.G.nodes[y])
                        and background_knowledge.is_forbidden(cg.G.nodes[y], cg.G.nodes[x])):
                    knowledge_ban_edge = True
                if knowledge_ban_edge:
                    if not stable:
                        edge1 = cg.G.get_edge(cg.G.nodes[x], cg.G.nodes[y])
                        if edge1 is not None:
                            cg.G.remove_edge(edge1)
                        edge2 = cg.G.get_edge(cg.G.nodes[y], cg.G.nodes[x])
                        if edge2 is not None:
                            cg.G.remove_edge(edge2)
                        append_value(cg.sepset, x, y, ())
                        append_value(cg.sepset, y, x, ())
                        break
                    else:
                        edge_removal.append((x, y))  # after all conditioning sets at
                        edge_removal.append((y, x))  # depth l have been considered

                Neigh_x_noy = np.delete(Neigh_x, np.where(Neigh_x == y))
                for S in combinations(Neigh_x_noy, depth):
                    if not (
                        data_sparse[(data_sparse['genus_A'] == cg.G.nodes[x]) & (data_sparse['genus_B'] == cg.G.nodes[y])].shape[0] == 1
                    ) or not (
                        data_sparse[(data_sparse['genus_A'] == cg.G.nodes[y]) & (data_sparse['genus_B'] == cg.G.nodes[x])].shape[0] == 1):
                        p = cg.ci_test(x, y, S)
                        if p > alpha:
                            if verbose:
                                print('%d ind %d | %s with p-value %f\n' % (x, y, S, p))
                            if not stable:
                                edge1 = cg.G.get_edge(cg.G.nodes[x], cg.G.nodes[y])
                                if edge1 is not None:
                                    cg.G.remove_edge(edge1)
                                edge2 = cg.G.get_edge(cg.G.nodes[y], cg.G.nodes[x])
                                if edge2 is not None:
                                    cg.G.remove_edge(edge2)
                                append_value(cg.sepset, x, y, S)
                                append_value(cg.sepset, y, x, S)
                                break
                            else:
                                edge_removal.append((x, y))  # after all conditioning sets at
                                edge_removal.append((y, x))  # depth l have been considered
                                for s in S:
                                    sepsets.add(s)
                        else:
                            if verbose:
                                print('%d dep %d | %s with p-value %f\n' % (x, y, S, p))
                if (x, y) in edge_removal or not cg.G.get_edge(cg.G.nodes[x], cg.G.nodes[y]):
                    append_value(cg.sepset, x, y, tuple(sepsets))
                    append_value(cg.sepset, y, x, tuple(sepsets))

        if show_progress:
            pbar.refresh()

        for (x, y) in list(set(edge_removal)):
            edge1 = cg.G.get_edge(cg.G.nodes[x], cg.G.nodes[y])
            if edge1 is not None:
                cg.G.remove_edge(edge1)

    if show_progress:
        pbar.close()

    return cg
    

def run_pc_depth2(data, data_sparse, alpha, fp):
    '''
    Run PC with a max depth of 2. 

    :param: data: pandas DataFrame
    :param: alpha: significance level

    :return: the causallearn object, CausalGraph
    '''
    cg_1 = skeleton_discovery_depth2(data = data.to_numpy(),
                                     data_sparse = data_sparse,
                                     alpha = alpha,
                                     indep_test = CIT(data.to_numpy(), 'fisherz'),
                                     node_names = list(data.columns))
    cg_2 = UCSepset.uc_sepset(cg_1, 2)
    cg = Meek.meek(cg_2)
    
    pyd = GraphUtils.to_pydot(cg.G, labels=list(data.columns))
    pyd.write_png(f"graphs/{fp}.png")
    
    return cg.G
    

def fisherz(genus_A, genus_B, genus_C, data):
    '''
    Run the Fisher Z (conditional) independence test between genus_A and genus_B conditioned on genus_c.

    :param: genus_A: index of the first genus to test independence
    :param: genus_B: index of the second genus to test independence
    :param: genus_C: index/list of index of the conditioning set
    :param: data: data containing the genus values
    '''
    fisherz_obj = CIT(data.values, "fisherz")
    pValue = fisherz_obj(genus_A, genus_B, genus_C)
    return pValue


def run_ouralg(adj_matrix, data, genus_lst, indep_test):
    '''
    Run our algorithm on the data and save the plot of the resulting causal graph. The methodology can be found in the associated paper. Returns the resulting adjacency matrix.

    :param: adj_matrix: adjacency matrix of the correlation graph from the paper
    :param: data: data for the corresponding cohort
    :param: genus_lst: list of genus corresponding to the adjacency matrix as node labels
    :param: indep_test: choice of independence test

    :return: the updated adjacency matrix
    '''
    adj_df = pd.DataFrame(adj_matrix, columns=genus_lst, index=genus_lst)

    
    def create_tuples(genus_lst, genus_A, genus_B):
        '''
        Creates tuples of 2 combinations for all items in genus_lst except genus_a and genus_b.
        
        :param: genus_lst: list of genus names
        :param: genus_A: index of the first genus to exclude
        :param: genus_B: index of the second genus to exclude
    
        :return: list of tuples, where each tuple contains two distinct genus names from genus_lst, excluding genus_a and genus_b.
        '''
        filtered_genus_lst = [x for x in range(len(genus_lst)) if (x != genus_A) and (x != genus_B)]
        combinations = list(itertools.combinations(filtered_genus_lst, 2))
        
        return combinations
        
    
    removed_edges = {}

    print('Conditioning set size: 1')
    for i in range(len(genus_lst)):
        print(genus_lst[i])
        for j in range(len(genus_lst)):
            if adj_df.iloc[i, j] == 1.0:
                genus_A = i
                genus_B = j
                for k in range(len(genus_lst)):
                    if k != i and k != j:
                        genus_C = [k]
                        try: 
                            p_value = indep_test(genus_A, genus_B, genus_C, data)
                        except:
                            p_value = 1
    
                        if p_value < 0.05: # no multiple testing correction
                            adj_df.iloc[i, j] = 0
                            removed_edges[(i, j)] = [genus_C]
                            break

    print('Conditioning set size: 2')
    for i in range(len(genus_lst)):
        print(genus_lst[i])
        for j in range(len(genus_lst)):
            if adj_df.iloc[i, j] == 1.0:
                genus_A = i
                genus_B = j
                combinations = create_tuples(genus_lst, genus_A, genus_B)
    
                for k in combinations:
                    genus_C = k
                    try :
                        p_value = indep_test(genus_A, genus_B, list(genus_C), data)
                    except:
                        p_value = 1
    
                    if p_value < 0.05:
                        adj_df.iloc[i, j] = 0
                        removed_edges[(i, j)] = [genus_C]
                        break
                    else:
                        adj_df.iloc[i, j] = 2
                        
    return adj_df.values
    

def graph_networkx(adj_matrix, genus_lst, fp):
    '''
    Create a NetworkX graph of the adjacency matrix and corresponding genus as node labels. Currently uses the circular layout for optimal viewing.

    :param: adj_matrix: adjacency matrix of the genus to genus correlations
    :param: genus_lst: list of genus corresponding to the adjacency matrix as node labels
    '''
    if 1 in np.unique(adj_matrix):
        edge_val = 1
        edge_color = "gray"

    if 2 in np.unique(adj_matrix):
        edge_val = 2
        edge_color = "black"
        
    rows, cols = np.where(adj_matrix == edge_val)
    edges = zip(rows.tolist(), cols.tolist())
    edges = zip(rows.tolist(), cols.tolist())
    edge_list = list(edges)
    graph_labels = {i: genus_lst[i]for i in range(len(genus_lst))}
    
    plt.figure(figsize=(40,30))
    G = nx.Graph()
    G.add_nodes_from(range(len(genus_lst)))
    G.add_edges_from(edge_list)
    pos_spaced = nx.circular_layout(G)
    nx.draw_networkx(G, pos_spaced, node_size=500, font_size=6, labels=graph_labels, with_labels=True, edge_color=edge_color, node_color="#b9d1f0")
    plt.savefig(f"graphs/{fp}.png")


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current Axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom", fontsize=12)

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels, fontsize=12)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels, fontsize=12)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-60, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="white", linestyle='-', linewidth=1)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

    