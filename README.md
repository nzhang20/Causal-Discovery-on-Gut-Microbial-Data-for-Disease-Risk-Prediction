# Causal Discovery on Gut Microbial Data for Disease Risk Prediction

Over the last few years, the scientific community has gained profound understanding of the importance of the gut microbiome to human health. However, most of these studies are purely associational, and as we all have learned "association is not causation". Thus, we are interested in tackling this causality question in the realm of gut microbiome research using causal discovery techniques such as constraint based causal graph search algorithms to understand 
1) the differences in microbial interactions between a diseased cohort and a healthy control cohort, and
2) the specific microbes that are directly linked to disease status.

Although we may claim to discover the causal structure of the gut microbiome for these two scenarios, conventionally, causality is established through randomized controlled trials. These can be costly and time-intensive, but we hope that our project can help aid in narrowing down the specific microbes to test for when conducting these wet-lab experiments. 

The difficult thing about doing research on the gut microbiome is that every healthy microbiome can look very different! The only (current) hallmark of a healthy microbiome is a diverse and balanced microbiome. To aid in precision medicine, we also build a prediction model to predict one's risk for disease given a certain microbiome. 

## File Descriptions

### Data
We have two data files corresponding to T2D and PCOS, respectively. In order to use BIRDMAn and be aligned with typical gut microbiome analysis where one has a table of the OTU counts, and a separate table for the metadata, each dataset is cleaned such that their OTU relative abundance values can be found in `data/{disease}/otu_table.csv` and the metadata file can be found in `data/{disease}/metadata.csv`, where `{disease}` is the disease of interest (T2D or PCOS).

(remove)
The PCOS dataset comes from the Yang 2024 systematic review, `data/raw/pcosyang2024.xlsx`. It contains normalized gut microbiome abundances for patients from 14 different PCOS studies conducted in China, Poland, Austria, and Russia. After data cleaning and feature pruning, the clean dataset can be found in `data/pcos/clean.csv`.

### Causal Discovery Algorithms & Graphs
To address research question 1), we develop our own causal discovery algorithm to reduce the number of conditional independence tests conducted, which can inflate our false discovery rate. 

To address research question 2), we implement Constraint-based causal Discovery from heterogeneous/NOnstationary Data (CD-NOD) where our domain index is the disease status. 

## Usage
Create and activate the conda envrionment with `conda env create -f environment.yml` and `conda activate capstone`. 

To replicate our results, run...
- `python run.py all`: runs the entire process from data cleaning to causal discovery algorithms and graphs
- `python run.py data`: recreates our clean dataset from the raw data
- `python run.py eda`: runs our EDA to produce plots assessing linearity, Gaussianity, and correlations between features
- `python run.py graph`: runs the causal discovery algorithms and generates results in the form of causal graphs
- `python run.py clean`: remove all generated files (e.g. the clean dataset, plots, and graphs)