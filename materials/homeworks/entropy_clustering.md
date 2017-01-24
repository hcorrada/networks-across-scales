---
title: Relative entropy and clustering
author: CMSC423 Fall 2015
date: September 29, 2015
geometry: margin=1in
fontfamily: utopia
fontsize: 11pt
papersize: letterpaper
---

**Question 1.** Consider profile $P$:

$$
\begin{matrix}
\mathrm{A:} & 0.0 & 0.0 & 0.0 \\
\mathrm{C:} & 0.5 & 0.5 & 0.5 \\
\mathrm{G:} & 0.2 & 0.3 & 0.5 \\
\mathrm{T:} & 0.3 & 0.2 & 0.0
\end{matrix}
$$

and background frequencies $b_A=0.3$, $b_C=0.3$, $b_G=0.2$ and $b_T=0.2$.

(a) What is the _entropy_ of profile $P$?  
(b) What is the _relative entropy_ of profile $P$ with respect to the background frequencies?

**Question 2.** In the soft k-means (EM) clustering algorithm:

(a) What are the parameters we want to estimate?  
(b) What does $\textrm{HiddenMatrix}_{ij}$ correspond to in this algorithm? E.g., it corresponds to the probability that _fill-in-the-blank_ $j$ is generated from _fill-in-the-blank_ $i$.    
(c) Given $\textrm{HiddenMatrix}_{ij}$ and data points $\textrm{Data}$, how is the $i$-th center calculated on the algorithm's M-step?
