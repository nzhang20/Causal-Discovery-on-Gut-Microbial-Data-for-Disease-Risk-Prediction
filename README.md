# Causal Discovery on Gut Microbial Data for Disease Risk Prediction

Over the last few years, the scientific community has gained profound understanding of the importance of the gut microbiome to human health. However, most of these studies are purely associational, and as we all have learned "association is not causation". Thus, we are interested in tackling this causality question in the realm of gut microbiome research using causal discovery techniques such as constraint based causal graph search algorithms to understand 
1) the differences in microbial interactions between a diseased cohort and a healthy control cohort, and
2) the specific microbes that are directly linked to disease status.

Although we may claim to discover the causal structure of the gut microbiome for these two scenarios, conventionally, causality is established through randomized controlled trials. These can be costly and time-intensive, but we hope that our project can help aid in narrowing down the specific microbes to test for when conducting these wet-lab experiments. 

The difficult thing about doing research on the gut microbiome is that every healthy microbiome can look very different! The only (current) hallmark of a healthy microbiome is a diverse and balanced microbiome. To aid in precision medicine, we also build a prediction model to predict one's risk for disease given a certain microbiome. 

## File Descriptions

### Data
We have two data files corresponding to T2D and PCOS, respectively. In order to use BIRDMAn and be aligned with typical gut microbiome analysis where one has a table of the OTU counts, and a separate table for the metadata, each dataset is cleaned such that their OTU relative abundance values can be found in `data/{disease}/otu_table.csv` and the metadata file can be found in `data/{disease}/metadata.csv`, where `{disease}` is the disease of interest (T2D or PCOS).

- `data/t2d` comes from the NIH HMP2 dataset, available at:
- `data/pcos` comes from a meta-analysis of 14 PCOS studies conducted in Asia and Europe, available at: 

### Causal Discovery Algorithms & Graphs
To address research question 1), we develop our own causal discovery algorithm to reduce the number of conditional independence tests conducted, which can inflate our false discovery rate. 

To address research question 2), we implement Constraint-based causal Discovery from heterogeneous/NOnstationary Data (CD-NOD) where our domain index is the disease status. 

## Usage
Create and activate the conda envrionment with `conda env create -f environment.yml` and `conda activate causal-gut-v2`. 

Before running the scripts, please confirm that your configuration file is correct. Please see `src/config-t2d.json` or `src/config-pcos.json` for examples on the T2D and PCOS dataset respectively. Below is a description of the configuration file:
- **otu_table_fp**: filepath to the otu table (filetype may be .csv or .biom)
- **metadata_fp**: filepath to the metadata table (filetype may be .csv or .tsv)
- **name**: name of the disease; this will be the name of the directory for outputs in `data`, `plots`, and `graphs`
- **sparse_microbe_microbe**: the sparse method to reduce microbe-microbe edges; may be one of "graphical lasso" or "sparcc"
- **disease_column**: column name in metadata corresponding to disease status
- **cohort_names**: names of the two disease status groups; the order in which names are presented will be mapped 0 and 1 respectively
- **covariates_map**: a dictionary of the column name of covariates and the names of the categories (order matters!)
- **rare_otu_threshold**: the percentage value to filter rare microbes
- **transformation**: transformation method; may be on of "norm" or "clr"


To replicate our results for t2d, run...
- `python run.py src/config-t2d.json all`: runs the entire process from data cleaning to causal discovery algorithms and graphs
- `python run.py src/config-t2d.json eda`: runs our EDA to produce plots assessing linearity, Gaussianity, and correlations between features
- `python run.py src/config-t2d.json graph`: runs the causal discovery algorithms and generates results in the form of causal graphs
- `python run.py src/config-t2d.json birdman`: runs BIRDMAn with the default negative binomial model 
- `python run.py src/config-t2d.json clean`: remove all generated files (e.g. the clean dataset, plots, and graphs) for t2d

To run a different configuration file, simply replace `src/config-t2d.json` with your configuration file. 

For causal inference results, you will need to consult the graphs produced and determine which microbes you want to include in a model of your choice. Our example analysis and results on T2D and PCOS can be found in `src/do-calculus`, where we use logistic regression to estimate causal effects, but there are many other causal inference estimators that can be used to answer this question. 