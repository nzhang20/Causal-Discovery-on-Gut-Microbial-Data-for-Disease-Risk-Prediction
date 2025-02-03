'''
eda.py contains functions to generate graphs for EDA
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import itertools


def indiv_per_cohort(data, disease, cohort_col, cohort_labels):
    '''
    Graphs a bar plot of the number of individuals per cohort (healthy vs diseased).

    :param: data: pandas DataFrame containing the relevant information
    :param: disease: name of the disease of interest
    :param: cohort_col: column name corresponding to the cohort in data
    :param: cohort_labels: list of cohort labels, e.g. ['HC', 'PCOS'] or ['IS', 'IR']
    '''
    fig, ax = plt.subplots()
    cohort_size = data.groupby(cohort_col).count().reset_index().iloc[:, 0:2]
    other_col = cohort_size.columns[cohort_size.columns != cohort_col][0]
    bars = ax.bar(cohort_size[cohort_col], cohort_size[other_col], width=0.7)
    ax.set_title('Individuals in Each Cohort')
    ax.set_xlabel('Cohort')
    ax.set_ylabel('Number of Individuals')
    ax.set_xticks([0, 1], labels=cohort_labels)
    ax.bar_label(bars)
    fig.savefig(f'plots/{disease}/cohort_bar_chart.png')
    plt.close()


def indiv_per_cohort_cov(data, disease, cohort_col, covariate_col, cohort_labels, covariate_labels):
    '''
    Graphs a grouped bar plot of the number of individuals in each cohort for each category of the covariate column.

    :param: data: pandas DataFrame containing the relevant information
    :param: disease: name of the disease of interest
    :param: cohort_col: column name corresponding to the cohort in data
    :param: covariate_col: column name corresponding to the covariate of interest in data
    :param: cohort_labels: list of cohort labels, e.g. ['HC', 'PCOS'] or ['IS', 'IR']
    :param: covariate_labels: list of covariate labels, e.g. ['Male', 'Female']
    '''
    fig, ax = plt.subplots()
    data_pivot = data.pivot_table(index=covariate_col, columns=cohort_col, aggfunc='size')
    data_pivot.plot(kind='bar', rot=0, ax=ax)
    ax.set_title(f'Individuals per {covariate_col.capitalize()}', fontsize=20)
    ax.set_xlabel(covariate_col.capitalize(), fontsize=16)
    ax.set_ylabel('Number of Individuals', fontsize=16)
    ax.set_xticks(range(data_pivot.shape[0]), labels=covariate_labels)
    for container in ax.containers: 
        ax.bar_label(container, fontsize=12)
    ax.legend(title='Cohort', labels=cohort_labels)
    fig.set_figheight(6)
    fig.set_figwidth(np.log(data_pivot.shape[0]) * 10)
    fig.savefig(f'plots/{disease}/{covariate_col.lower()}_bar_chart.png')
    plt.close()


def check_linearity(data, disease):
    '''
    Graphs the scatter plot of every pair of genera in the dataset to assess linearity.

    :param: data: pandas DataFrame containing only the genera 
    :param: disease: name of the disease of interest
    '''
    pairs = list(itertools.combinations(data.columns, 2))

    for i in range(len(pairs)):
        fig, ax = plt.subplots()
        x = pairs[i][0]
        y = pairs[i][1]
        ax.scatter(x = data[x], y = data[y])
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_title(f'Scatter Plot of {x} and {y}')
        fig.savefig(f'plots/{disease}/linearity/{x}_vs_{y}.png')
        plt.close()


def check_normality(data, disease):
    '''
    Graphs the qqplots of every genera in the dataset to assess normality.

    :param: data: pandas DataFrame containing only the genera
    :param: disease: name of the disease of interest
    '''
    for i in range(len(data.columns)):
        fig, ax = plt.subplots()
        sm.qqplot(data.iloc[:, i], line='q', fit=True, ax=ax)
        ax.set_title(f'QQPlot of {data.columns[i]}')
        fig.savefig(f'plots/{disease}/normality/{data.columns[i]}.png')
        plt.close()


def plot_corr_heatmap(data, disease):
    '''
    Plots the correlation matrix heatmap for the genera in the dataset.

    :param: data: pandas DataFrame containing only the genera
    :param: disease: name of the disease of interest
    '''
    genera = list(data.columns)
    figsize = min(40, len(genera) / np.log(len(genera)))
    fig, ax = plt.subplots(figsize=(figsize, figsize))
    cax = ax.matshow(abs(data.corr(method='spearman')), cmap='Greens')
    fig.colorbar(cax, fraction=0.046, pad=0.04)
    
    xaxis = np.arange(len(genera))
    ax.set_xticks(xaxis, labels=genera)
    ax.set_yticks(xaxis, labels=genera)
    ax.set_title(f'Spearman Correlation Matrix for {disease}')
    
    plt.setp(ax.get_xticklabels(), rotation=-60, ha='right', rotation_mode='anchor')
    
    fig.tight_layout()
    fig.savefig(f'plots/{disease}/correlation_heatmap.png')
    plt.close()

    