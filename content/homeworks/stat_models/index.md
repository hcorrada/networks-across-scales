---
title: "Homework: Statistical Analysis of Network Data"
date: "2017-10-24"
---

**DUE**: Monday 11/4/2019, 11:59pm  
**Posted**: 10/21/2019  
**Last Update**: 10/21/2019  


We will use ecological network data from the paper "How Structured is the Entangled Bank?..." by Kefi et al. http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002527. This paper used Stochastic Block Models to analyze both trophic and non-trophic species interaction networks. Your goal in this homework is to (a) partially replicate some of the analyses in the paper, and reanalyze this data using other statistical methods for analysis. 


## Data 

Data for the trophic and non-trophic networks along with species metadata is available for download http://datadryad.org/resource/doi:10.5061/dryad.b4vg0. Download the adjacency matrix for the trophic network and species metadata.

## Stochastic Block Models

Fit a stochastic block model to the trophic network. 

In R, use the `blockmodels` package: https://cran.r-project.org/package=blockmodels 

In Python, you can use the `graph-tool` https://graph-tool.skewed.de/ library. More information about SBM and extensions found here: https://graph-tool.skewed.de/static/doc/demos/inference/inference.html#the-stochastic-block-model-sbm.  

### (Qual only)

Implement an inference procedure for the SBM model. You have two options:

- Variational EM algorithm: this is what `blockmodels` implements. The reference implementation is described in this paper: https://arxiv.org/pdf/1011.1813.pdf. The preprint for the `blockmodels` package has a smaller introduction: https://arxiv.org/pdf/1602.07587.pdf

- MCMC: this is what `graph-tool` implements. The reference implementation is described in this paper https://arxiv.org/pdf/1310.4378.pdf

## Non-probabilistic modularity methods

Use a non-probabilistic method based on modularity to find network communities on the trophic network. You can use either Girvan-Newman or Blondel et al. (Louvain method) https://arxiv.org/abs/0803.0476. Both of these are implemented in `igraph` and their corresponding R or python interfaces.

## Report 

Report on the resulting structure of the class membership probability distributions of the species using SBM. Comment on how your result relates to the results reported in the paper. Compare the resulting communities with that obtained based on modularity methods.

For both SBM and modularity communities, is there any correlation between class membership and vertex attributes reported for these species? One way to answer the latter is to perform a regression analysis modeling vertex attributes dependent on class membership, (for SBM, use the most likely class label for each species).

## Submission

As before, use Rmarkdown document or Jupyter notebook to prepare your submission including all code and discussion. Knit to pdf and submit to ELMS.
