---
title: "Homework: Statistical Analysis of Network Data"
date: "2017-10-24"
---

**DUE**: Monday 11/5/2018, 11:59pm  
**Posted**: 10/25/2018  
**Last Update**: 10/30/2018  


We will use ecological network data from the paper "How Structured is the Entangled Bank?..." by Kefi et al. http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002527. This paper used Stochastic Block Models to analyze both trophic and non-trophic species interaction networks. Your goal in this homework is to (a) partially replicate some of the analyses in the paper, and reanalyze this data using other statistical methods for analysis. 


## Data 

Data for the trophic and non-trophic networks along with species metadata is available for download http://datadryad.org/resource/doi:10.5061/dryad.b4vg0. Download the adjacency matrix for the trophic network and species metadata.

## Stochastic Block Models

Fit a stochastic block model to the trophic network. You can use the `mixer` R package as described in SAND Ch. 6. for this task. Alternatively, you can use any appropriate implementation of the EM algorithm and a suitable model selection criterion (AIC, BIC, or ICL). **UPDATE: The `mixer` package is not available for R 3.5, use package `blockmodels` instead.**

In Python, you can use the `graph-tool` https://graph-tool.skewed.de/. More information about SBM and extensions found here: https://graph-tool.skewed.de/static/doc/demos/inference/inference.html#the-stochastic-block-model-sbm.  

Report on the resulting structure of the class membership probability distributions of the species. Comment on how your result relates to the results reported in the paper. Is there any correlation between class membership and vertex attributes reported for these species? One way to answer the latter is to perform a regression analysis modeling vertex attributes dependent on class membership, taking the most likely class label for each species.

## Exponential Random Graph Models

Fit an ERGM to the trophic network. You should use the `ergm` R package as described in SAND Ch. 6 for this task. Fit a model that includes both network metrics (for example, number of edges and geometrically weighted edge-wise sharing pairs (gwesp)) and vertex attributes. Unfortunately, there is no easy-to-use package for ERGM in python. One option is to use R within python using the `rpy2` (or `rpy`) packages. 

- Do any of the vertex attributes have a significant effect in the log odds of edge presence in this network? 

- How does this result depend on your choice of network metric included in the model. To answer this question fit a model with a different set of network metrics as covariate (based on the degree distribution for example) and compare to the previous result.

For both of these make sure to include an analysis of goodness of fit as part of your discussion.

## Submission

Submit a pdf including code and answers to questions above to ELMS.
