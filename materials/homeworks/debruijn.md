---
title: Eulerian Path Exercise
author: CMSC423 
fontfamily: utopia
geometry: margin=1in
---

Name(s):   
UID(s):  

**Question 1.** Solve the string reconstruction problem for this set of eight 3-mers:

$$
\{\mathtt{AGT,AAA,ACT,AAC,CTT,GTA,TTT,TAA}\}
$$

(a) Construct the DeBruijn graph with 8 edges corresponding to these 3-mers (string overlap, Eulerian path approach)

(b) Find a Eulerian path (8 edges) which visits each edge exactly once. Does this path visit every vertex of the graph at most one time? 

(c) Write the reconstructed string corresponding to this Eulerian path.

**Question 2.** Consider nodes $a=$`ACCTG` and $b=$`CCTGT` in a DeBruijn graph $G$ built from the $k$-mer composition of string $s$.
Suppose graph $G$ contains 5 edges connecting node $a$ to node $b$. How many times does $k$-mer `ACCTGT` appear in string $s$. 
