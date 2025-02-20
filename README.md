# Causal Discovery on Gut Microbial Data for Disease Risk Prediction

## Authors
**Mariana Paco Mendivil**, **Candus Shi**, **Nicole Zhang**  
**Mentors**: Dr. Biwei Huang, Dr. Jelena Bradic  

## Background
### Association vs. Causation
Traditional statistical models identify associations, but they often fail to capture true causal relationships. Conducting randomized experiments is expensive, slow, and sometimes unethical. Causal discovery and inference techniques allow us to infer causal relationships from observational data.

### Gut Microbiome and Disease
The gut microbiome plays a crucial role in human health. However, due to its complexity and heterogeneous effects across populations, studying its causal relationships with diseases requires advanced computational methods.

### Causal Discovery in Gut Microbiome Research
Existing research has explored causal inference in the gut microbiome using algorithms like PC-Stable and do-calculus. More recent methods, such as CD-NOD, handle heterogeneous data, making them particularly useful for gut microbiome studies where data come from diverse sources.

## Research Questions
1. **Microbe-Microbe Interactions**: How do microbial interaction networks differ between healthy and diseased participants?
2. **Microbe-Disease Causal Relationships**: Which microbes have a causal effect on disease status?
3. **Prediction**: Can we predict disease status using causal representation learning?

## Data
We use gut microbiome datasets related to **Type 2 Diabetes (T2D)** and **Polycystic Ovary Syndrome (PCOS)**:

- **T2D Data**: NIH Human Microbiome Project (HMP2) dataset (filtered for healthy visits, 16S sequencing).  
  - **Samples**: 153 insulin-sensitive (IS) & 178 insulin-resistant (IR).

- **PCOS Data**: Aggregated from 14 clinical studies (16S sequencing).  
  - **Samples**: 435 healthy controls (HC) & 513 PCOS patients.

## Causal Discovery Methods
We use **causal discovery algorithms** alongside predictive modeling to analyze microbe-microbe and microbe-disease relationships. Our methodology includes:

1. **Feature Selection**: Reduce data dimensionality using sparse estimation and sure-screening.
   - Remove rare OTUs (operational taxonomic units) with <1% relative abundance.
   - Apply **SparCC** and **Graphical Lasso** to filter microbe-microbe interactions.
   - Use **Lasso Logistic Regression** to identify disease-relevant microbes.

2. **Causal Discovery Algorithms**:
   - **PC-Stable** (max depth = 2) for microbe-microbe networks.
   - **CD-NOD** to model heterogeneity in microbe-disease relationships.

3. **Causal Effect Estimation**:
   - Use **do-calculus** and logistic regression to estimate individual microbial effects.

## Key Findings
### Type 2 Diabetes (T2D)
From our **microbe-disease network**, we identified five genera with causal relationships to T2D:
- **Butyricimonas**, **Clostridium XIVb**, **Odoribacter**, **Unclassified Bacteria**, **Unclassified Firmicutes**

Using do-calculus, we estimated their causal effects on T2D. (See Table 1 in the original poster for statistical results.)

### Polycystic Ovary Syndrome (PCOS)
Our **microbe-disease network** for PCOS identified nine causal genera:
- **Alistipes, Blautia, Burkholderia, Desulfovibrio, Holdemanella, Knoellia, Prevotellaceae NK3B31 group, Ruminococcus, Ruminococcus gnavus group**

Causal effects were further validated using do-calculus (See Table 2 in the original poster).

## Conclusion & Future Work
- Answered key research questions on microbial interactions and disease risk.
- Developed a framework combining **causal discovery** and **predictive modeling** for microbiome analysis.
- Future work includes improving **causal representation learning** using Variational Autoencoders (VAE) and integrating **BIRDMAn** for Bayesian inference.

## References

## More Information
Visit our project GitHub:  
[https://github.com/nzhang20/Causal-Discovery-on-Gut-Microbial-Data-for-Disease-Risk-Prediction/](https://github.com/nzhang20/Causal-Discovery-on-Gut-Microbial-Data-for-Disease-Risk-Prediction/)
