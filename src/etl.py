'''
etl.py contains functions used to clean the raw dataframes and transforming data into adjacency matrices
'''

import pandas as pd
import numpy as np
import os


def clean_data_t2d(subject_fp, sample_fp, gut_fp):
    '''
    Returns the otu table and metadata for the T2D study.

    :param: subject_fp: filepath to the subject dataset
    :param: sample_fp: filepath to the sample info dataset
    :param: gut_fp: filepath to the gut abundance dataset

    :return: pandas DataFrame of the otu table, pandas DataFrame of the metadata
    '''
    subjects = pd.read_csv(subject_fp)
    samples = pd.read_csv(sample_fp)
    gut_16s = pd.read_csv(gut_fp, sep='\t')

    # filter out subjects who had an unknown "disease" (IR/IS) status
    subjects_IRIS_known = subjects[subjects['IRIS'] != 'Unknown'][['IRIS', 'Gender', 'Ethnicity', 'SubjectID']]

    # filter out samples that were not on a "Healthy" visit
    samples_healthy = samples[(samples['Gut_16S'] == 1) & (samples['CL4'] == 'Healthy')][['SubjectID', 'SampleID']]

    # only get genera, transform into percentages
    genera = []
    for col in gut_16s.columns:
        if 'genus_' in col:
            genera.append(col)
    gut_16s_genera = gut_16s[genera] * 100
    gut_16s_genera = pd.concat([gut_16s[['SampleID']], gut_16s_genera], axis=1)
    

    # merge three dataframes
    merged_df = pd.merge(gut_16s_genera, samples_healthy, on='SampleID', how='inner')
    merged_df = pd.merge(subjects_IRIS_known, merged_df, on='SubjectID', how='inner')

    # remove 'SubjectID' and set index to 'SampleID'
    merged_df = merged_df.drop(columns=['SubjectID'])
    merged_df = merged_df.set_index('SampleID')

    # convert categories to numbers 
    merged_df['IRIS'] = merged_df['IRIS'].map({'IS': 0, 'IR': 1, 'Unknown': 2})
    merged_df['Gender'] = merged_df['Gender'].map({'M': 0, 'F': 1})
    merged_df['Ethnicity'] = merged_df['Ethnicity'].map({'C': 0, 'A': 1, 'B': 2, 'H': 3, 'unknown': 4})

    otu_table = merged_df.loc[:, [x for x in merged_df.columns if 'genus_' in x]]
    metadata = merged_df.loc[:, ['IRIS', 'Gender', 'Ethnicity']]
    return otu_table, metadata


def clean_data_pcos(raw_fp):
    '''
    Returns the otu table and metadata for the PCOS study.

    :param: raw_fp: filepath to the raw dataset

    :return: pandas DataFrame of the otu table, pandas DataFrame of the metadata
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

    otu_table = pcosyang2024.loc[:, [x for x in pcosyang2024.columns if (x != 'group') & (x != 'region') & (x !='study_site')]]
    metadata = pcosyang2024.loc[:, ['group', 'region', 'study_site']]
    return otu_table, metadata


def filter_rare_otus(data, k=1):
    '''
    Filters the dataset for rare OTUs based on a chosen threshold. 

    :param: data: clean dataset, only containing abundance values
    :param: k: rare OTU abundance threshold (in percentage), default = 1%

    :return: a pandas DataFrame
    '''
    rare_OTUs = data.columns[(data <= k).sum() == data.shape[0]]
    filter_rare = data.drop(columns=rare_OTUs)
    return filter_rare