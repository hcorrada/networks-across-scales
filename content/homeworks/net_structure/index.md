---
title: "Homework: Network Structure Inference"
date: "2017-11-16"
---

**Posted**: 11/6/2018  
**Last Update**: 11/6/2018

## Data

Download timeseries data derived from the ADHD-200 resting fMRI samples (http://fcon_1000.projects.nitrc.org/indi/adhd200/).

[{{< baseurl >}}data/roi_timeseries.csv]({{< baseurl >}}/data/roi_timeseries.csv).

The data consists of standardized measurements over 176 time points of 39 regions of interest (ROI). The rows of the data correspond to the 39 regions, the first column contains a label indicating the name of the brain ROI, with the remaining columns corresponding to the fMRI timeseries.

## Correlation Network

Build a correlation network from the timeseries data as follows:

- Compute correlation (Pearson) \\(\hat{\rho}_{ij}\\) for every pair of ROIs
- Use Fisher's transform to compute a $z$-value from the estimated correlation

<div>
$$
z_{ij} = \frac{1}{2} \log{\frac{(1+\hat{\rho}_{ij})}{(1-\hat{\rho}_{ij})}}
$$
</div>

- Compute a $p$-value for the hypothesis test of no correlation between each pair of regions using a normal distribution with mean 0 and variance $1/(n-3)$ with $n$ the number of timepoints (see SAND CH. 7).

- Use a multiple testing correction procedure to derive a final list of network edges.

## Partial Correlation Network

Next, build a partial correlation network from the timeseries data. Use partial correlations between ROIs $i$ and $j$ as follows:

- For every ROI $k \neq i$ and $k \neq j$, compute partial correlation $\hat{\rho}_{ij|k}$ as the correlation between vectors $r_i$ and $r_j$, where $r_i$ are the residuals of linear model $x_i = b_0 + b_1 x_k$ where $x_i$ and $x_k$ are the timeseries for ROIs $i$ and $k$ respectively.

- Compute a $p$-value $p_{ij|k}$ for the hypothesis test of no correlation between $i$ and $j$ (now conditioned on $k$) using Fisher's transformation and a normal distribution with mean $0$ and variance $1/(n-1-3)$.

- Combine the $p$-values for $i$ and $j$ using the Wille rule

<div>
$$
p_{ij} = \max \{ p_{ij|k} | \forall k \neq i,j \}
$$
</div>

- Use a multiple testing correction procedure as above to derive the final list of network edges

## Graphical Lasso

Implement (or use a library) the graphical lasso procedure to estimate a sparse inverse covariance matrix with entries $\omega_{ij}$. 

- Recall that the sparse inverse covariance estimate depends on a penalty parameter $\lambda$. Build at least two networks with $\lambda$ values that return a very sparse and a very dense network. Additionally, choose a model selection criterion (see SAND Ch.7 pg. 128 for initial ideas), and use the network derived with the optimal $\lambda$ as defined by your chosen criterion.

- Note, you can derive partial correlation coefficients for $i$ and $j$ from the inverse covariance matrix as

You can find implementations of the graphical lasso in R package `glmnet` and in python module `sklearn.covariance`.

<div>
$$
\hat{\rho}_{ij|V\setminus\{i,j\}} = \frac{-\omega_{ij}}{\sqrt{\omega_{ii}\omega_{jj}}}
$$
</div>

## Discussion

1. Compare the networks derived from these procedures using relevant network structure statistics, e.g., vertex degree distribution. Support with appropriate figures.

2. Provide a qualitative comparison of the derived networks. E.g., are there brain ROIs that are consistently linked in all networks?

## Submission

Submit a writeup describing your methods (e.g., multiple testing correction procedure used, etc.) and discussion. Include code as an appendix, or within the writeup if using jupyter notebook, Rmarkdown or similar.



