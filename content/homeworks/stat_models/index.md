---
title: "Homework: Statistical Analysis of Network Data"
date: "2017-10-24"
---

**DUE**: Wednesday 11/1/2017, 11:59pm  
**Posted**: 10/24/2017  
**Last Update**: 10/24/2017  


We will use ecological network data from the paper "How Structured is the Entangled Bank?..." by Kefi et al. http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002527. This paper used Stochastic Block Models to analyze both trophic and non-trophic species interaction networks. Your goal in this homework is to (a) partially replicate some of the analyses in the paper, and reanalyze this data using other statistical methods for analysis. 


## Data 

Data for the trophic and non-trophic networks along with species metadata is available for download http://datadryad.org/resource/doi:10.5061/dryad.b4vg0. Download the adjacency matrix for the trophic network and species metadata.

## Stochastic Block Models

Fit a stochastic block model to the trophic network. You can use the `mixer` R package as described in SAND Ch. 6. for this task. Alternatively, you can use any appropriate implementation of the EM algorithm and a suitable model selection criterion (AIC, BIC, or ICL) to perform this task. 

Report on the resulting structure of the class membership probability distributions of the species. Comment on how your result relates to the results reported in the paper. Is there any correlation between class membership and metadata reported for these species? 

## Exponential Random Graph Models

Fit an ERGM to the trophic network. You should use the `ergm` R package as described in SAND Ch. 6 for this task. Fit a model that includes both network metrics and metadata covariates. 

- Do any of the metadata have a significant effect in the log odds of edge presence in this network? 

- How does this result depend on your choice of network metric included in the model. To answer this question fit a model with a different set of network metrics as covariate and compare to the previous result.

For both of these make sure to include an analysis of goodness of fit as part of your discussion.

## Submission

Submit a pdf including code and answers to questions above to ELMS.
