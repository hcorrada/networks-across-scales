---
layout: default
title: Homework 3
---

**Due Monday Nov. 10 in class**

## Programming Questions ##

Submit your answers to Problems 9 and 10 in the
[Rosalind final submission page](http://rosalind.info/classes/156/).
**THESE ARE DUE ON FRIDAY NOV.  7**

## Clustering ##

**Question 1**
Prove the approximation bound for farthest-first traversal to the $$k$$-centers problem:

Let $$X$$ be a solution found by FarthestFirstTraversal and $$X_{opt}$$ be the optimal solution of the $$k$$-Center Clustering Problem. Prove that $$\mathrm{MaxDistance}(\mathrm{Data},X) ≤ 2 \times \mathrm{MaxDistance}(\mathrm{Data},X_{opt})$$.

## Assembly ##

For this exercise you will use your code for string reconstruction using DeBruijn graphs to analyze genomes from the Ebola virus. Download data for this exercise at this URL:

[http://cbcb.umd.edu/~hcorrada/CMSC423/downloads/hw3data.zip](http://cbcb.umd.edu/~hcorrada/CMSC423/downloads/hw3data.zip)

Unzip this file in your working directory, you should see the following four files:

- `ebola.fa`: This FASTA file contains assembled genomes from twenty different patients sampled during the 2014 ebola outbreak. These are a subset of the genomes analyzed as part of the sequencing effort described in this Science paper: [http://www.sciencemag.org/content/345/6202/1369.full](http://www.sciencemag.org/content/345/6202/1369.full). This effort was chronicled in this New Yorker piece: [http://www.newyorker.com/magazine/2014/10/27/ebola-wars](http://www.newyorker.com/magazine/2014/10/27/ebola-wars). All data related to this effort can be obtained from the NCBI [http://www.ncbi.nlm.nih.gov/bioproject/257197](http://www.ncbi.nlm.nih.gov/bioproject/257197).
- `mystery_paired_200.txt`: Simulated paired-end reads with $$d=200$$ from the same genome. Both files containing paired-end reads use the format in Rosalind:

`...AGACATCCGAACCATAGAGGATTC|CTCTATGTGCTGTGATG...`

**Question 2** Assemble the ebola genome from the given reads and see if it matches any of the twenty genomes in the `ebola.fa` file. Write down the sequence name from the `ebola.fa` file that matches you assembled genome (it starts with 'KM')

**Question 3** How hard is it to assemble this genome from these simulated 101 bp reads? Provide the  distribution of in-degrees and out-degrees in the DeBruijn graph for this set of reads. 

**Question 4** Generate all 101-mers from this genome using your 'String Composition' (Problem 38 in the pre-lecture exercises). How repetitive is this genome when using 101 bp reads?

**Question 5** What is the largest $$k$$ for which the $$k$$-mer composition of this genome contains at least 1 repeat.

## How to submit ##

Submit your answers to the five questions above in writing before
lecture on **November 10**. Include a listing of additional code you used to answer these questions (e.g., how did you read the `ebola.fa` file and checked if your assembled genome matches any of the twenty genomes there).



