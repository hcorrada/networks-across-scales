---
title: HW3 Question 1 Proof
author: CMSC 423 Fall 2014
date: Nov. 10 2014
options: -s
geometry: margin=1in
fontsize: 11pt
header-includes:
  - \newcommand{\md}{\mathrm{MaxDist}}
  - \newcommand{\dat}{\mathrm{Data}}
---

(A) First we show that the $\md$ quantity is non-increasing as more centers are added in farthest-first traversal. Let $X^{k-1}$ be the set of centers selected in the first $k-1$ iterations of farthest-first traversal, and $x_k$ be the center selected in the $k$-th step, then $\md(\dat,X) \leq \md(\dat,X^{k-1})$.

\begin{eqnarray}
\md(\dat,X) & = & \max_{u \in \dat} \min_{x \in X} d(u,x) \\
{} & = & \max_{u \in \dat} \min \left[ \min_{x \in X^{k-1}} d(u,x), d(u,x_k) \right] \\
{} & \leq & \max_{u \in \dat} \min_{x \in X^{k-1}} d(u,x) \\
{} & = & \md(\dat,X^{k-1})
\end{eqnarray}

From this it follows that $\md(\dat,X) \leq \md(\dat,X^t)$ for all t=\[1,\ldots,k-1\].

(B) Using (A), we show that $d(x_i,x_j) \geq \md(\dat,X)$ for all $x_i,x_j \in X$. 

Assume $i < j$, that is, $x_i$ was chosen as a center in an earlier iteration than $x_j$. Then,

$$
d(x_i,x_j) \geq \min_{x_t:t=\left[1,\ldots,j-1\right]} d(x_t, x_j) = \md(\dat,X^{j-1})\geq \md(\dat,X)
$$

where the last inequality follows from (A).

(C) Let $u \in \dat$ be such that $d(u,X)=\md(\dat,X)$, then by definition $d(u,x)\geq \md(\dat,X)$ for all $x\in X$.

(D) From (B) and (C), it follows that $X \cup \{u\}$ is a set of $k+1$ points in $\dat$ with distance between every pair of points greater than or equal to $\md(\dat,X)$.

(E) The optimal $k$ clustering must include two of the points defined in (D) in one of it's clusters. Therefore, $\md(\dat,X_{opt}) \geq \frac{\md(\dat,X)}{2}$. Which proves the result.




