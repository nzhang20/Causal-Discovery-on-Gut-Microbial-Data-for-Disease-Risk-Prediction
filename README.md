# Causal Discovery on Gut Microbial Data for Disease Risk Prediction

**Mariana Paco Mendivil, Candus Shi, Nicole Zhang**

## Overview
Did you know the bacteria in our gut play a huge role in our health, including diseases like Type 2 Diabetes (T2D)? These gut bacteria, also known as microbes, hold immense potential for creating personalized treatments and aiding in the early detection of diseases to prevent their onset. 

However, most existing studies on the gut microbiome are purely associational or provide correlational interpretations. While statistically valid, these approaches fall short of explaining the underlying science or revealing the causal mechanisms behind why certain microbes play a more significant role. The high-dimensional nature of omics data adds to the complexity, with limited sample sizes and numerous variables per sample often complicating traditional multiple testing procedures and asymptotic assumptions.

Our project uses causal discovery and causal representation learning techniques to go beyond correlations. We aim to uncover the deeper cause-and-effect relationships within gut microbes and between microbes and diseases. By finding the most critical microbial connections, we hope to predict and prevent diseases more effectively, paving the way for groundbreaking discoveries in health research.

## Problem Statement and Motivation
While research has revealed associations between the gut microbiome and diseases like T2D and polycystic ovary syndrome (PCOS), these findings do not illuminate the causal mechanisms involved. Understanding these mechanisms is crucial for developing targeted therapies and improving diagnostic precision. Yet, analyzing gut microbiome data is inherently challenging due to its high dimensionality, variability across populations, and clustering caused by diverse study designs.

Our motivation lies in addressing these limitations by:
1. Utilizing causal discovery methods to analyze microbiome data rigorously.
2. Identifying microbes directly influencing disease states to enhance therapeutic strategies.
3. Facilitating precision medicine by uncovering meaningful causal pathways in microbiome-disease interactions.

By addressing these gaps, our work aims to provide actionable insights that can transform microbiome research and its applications in medicine.

## Objectives
Our project seeks to achieve the following objectives:

1. **Discover Causal Structures**: Apply advanced causal discovery algorithms to identify relationships between gut microbes and diseases such as T2D and PCOS.

2. **Compare Microbial Interactions**: Analyze differences in microbial interactions and causal pathways between healthy and diseased cohorts to detect disease-specific patterns.

3. **Identify Microbial Biomarkers**: Highlight specific microbes with direct causal links to disease states, enabling the development of targeted diagnostics and treatments.

4. **Develop Predictive Models**: Build a robust prediction model that uses gut microbiome composition to estimate disease risk, leveraging causal insights for enhanced precision.

By fulfilling these objectives, we aim to advance understanding of the gut microbiome’s role in metabolic health and lay the groundwork for innovative clinical applications.

---

Over the last few years, the scientific community has gained profound understanding of the importance of the gut microbiome to human health. However, most of these studies are purely associational, and as we all have learned "association is not causation." Thus, we are interested in tackling this causality question in the realm of gut microbiome research using causal discovery techniques such as constraint-based causal graph search algorithms to understand:

1. The differences in microbial interactions between a diseased cohort and a healthy control cohort.
2. The specific microbes that are directly linked to disease status.

Although we may claim to discover the causal structure of the gut microbiome for these two scenarios, conventionally, causality is established through randomized controlled trials. These can be costly and time-intensive, but we hope that our project can help aid in narrowing down the specific microbes to test for when conducting these wet-lab experiments.

The difficult thing about doing research on the gut microbiome is that every healthy microbiome can look very different! The only (current) hallmark of a healthy microbiome is a diverse and balanced microbiome. To aid in precision medicine, we also build a prediction model to predict one's risk for disease given a certain microbiome.

---

## File Descriptions

### Data
We have one data file from the Yang 2024 systematic review, `data/raw/pcosyang2024.xlsx`. It contains normalized gut microbiome abundances for patients from 14 different PCOS studies conducted in China, Poland, Austria, and Russia. After data cleaning and feature pruning, the clean dataset can be found in `data/clean.csv`.

### Causal Discovery Algorithms & Graphs
To address research question 1), we develop our own causal discovery algorithm to reduce the number of conditional independence tests conducted, which can inflate our false discovery rate.

To address research question 2), we implement Constraint-based causal Discovery from heterogeneous/Nonstationary Data (CD-NOD), where our domain index is the disease status.

### SparCC
To reduce the number of edges from a complete graph, we implement SparCC to find significant correlations between pairs of microbes. The `SparCC` folder comes directly from the https://github.com/dlegor/SparCC/tree/master repository, cloned on January 29th 2025. There are a few tweaks to allow for running the method from the top directory. 

---

## Usage

Create and activate the conda environment with the following commands:

```bash
conda env create -f environment.yml
conda activate capstone
```

To replicate our results, run:

- `python run.py all`: Runs the entire process from data cleaning to causal discovery algorithms and graphs.
- `python run.py data`: Recreates our clean dataset from the raw data.
- `python run.py eda`: Runs our EDA to produce plots assessing linearity, Gaussianity, and correlations between features.
- `python run.py graph`: Runs the causal discovery algorithms and generates results in the form of causal graphs.
- `python run.py clean`: Removes all generated files (e.g., the clean dataset, plots, and graphs).


Asked ChatGPT

but also add the other info in our  report pls
