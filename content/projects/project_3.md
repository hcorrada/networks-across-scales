---
date: 2017-02-07
title: "Project 3: Genome Assembly and String Alignment"
---

**Due: Tuesday April 18**  
**Posted: April 4, 2017**  
**Last Update: April 4, 2017**  



## Programming Questions ##

Submit your answers to Problems 6 and 7 in the
[Rosalind final submission page](http://rosalind.info/classes/401/).

### Code Grading (60 pts) ###

Same guidelines as [Project 1](projects/projects_1/).

## Assembly ##

For this exercise you will use your code for string reconstruction using DeBruijn graphs to analyze genomes from the Ebola virus. Download data for this exercise at this URL:

[http://cbcb.umd.edu/~hcorrada/CMSC423/downloads/hw3data.zip](http://cbcb.umd.edu/~hcorrada/CMSC423/downloads/hw3data.zip)

Unzip this file in your working directory, you should see the following four files:

- `ebola.fa`: This FASTA file contains assembled genomes from twenty different patients sampled during the 2014 ebola outbreak. These are a subset of the genomes analyzed as part of the sequencing effort described in this Science paper: [http://www.sciencemag.org/content/345/6202/1369.full](http://www.sciencemag.org/content/345/6202/1369.full). This effort was chronicled in this New Yorker piece: [http://www.newyorker.com/magazine/2014/10/27/ebola-wars](http://www.newyorker.com/magazine/2014/10/27/ebola-wars). All data related to this effort can be obtained from the NCBI [http://www.ncbi.nlm.nih.gov/bioproject/257197](http://www.ncbi.nlm.nih.gov/bioproject/257197).
- `mystery_paired_200.txt`: Simulated paired-end reads with \\(d=200\\) from the same genome. Both files containing paired-end reads use the format in Rosalind:

`...AGACATCCGAACCATAGAGGATTC|CTCTATGTGCTGTGATG...`

**Question 1** (15pts) Assemble the ebola genome from the given reads and see if it matches any of the twenty genomes in the `ebola.fa` file. Write down the sequence name from the `ebola.fa` file that matches you assembled genome (it starts with 'KM')

**Question 2** (10 pts) How hard is it to assemble this genome from these simulated 101 bp paired-end reads? Provide the  distribution of in-degrees and out-degrees in the DeBruijn graph for this set of reads.

**Question 3** (5 pts) Generate all 101-mers from this genome using your 'String Composition' (Problem 21 in the pre-lecture exercises). How repetitive is this genome when using 101 bp reads?

**Question 4** (10 pts) What is the largest \\(k\\) for which the \\(k\\)-mer composition of this genome contains at least 1 repeat.

## How to submit ##

Writeup answers to the four questions above. Submit to ELMS along with *all* code you used to answer these questions.
