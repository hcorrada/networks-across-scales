---
title: "Homework: Practice with Network Statistics"
date: "2017-09-02"
---


**DUE**: Friday 9/13/2017, 11:59pm  
**Posted**: 9/1/2017  
**Last Update**: 9/4/2017

## Data 

Obtain data for the year genetic interaction network at http://thecellmap.org/costanzo2016/. Download 
data in matrix format from http://thecellmap.org/costanzo2016/data_files/Raw%20genetic%20interaction%20datasets:%20Matrix%20format.zip, 
and after unzipping use file `Data File S2. Raw genetic interaction datasets: Matrix format/SGA_NxN_clustered.cdt` This is a `tab-separated value` file
with the following format:

- The first five rows contain metadata. Rows six to the end contain interaction data
- The first six columns also contain metadata. Columns seven to the end contain interaction datasets
- Each column corresponds to an array containing all double knockdowns for a specific gene, with each row corresponds to the second gene knocked-down (the query gene).
For example, column seven corresponds to the array with knockdowns of gene `PDR1`, and the entry in row 6 in column 7 corresponds to knocking down gene `RSC1` in that array.

For simplicity, let's make a few transformations to turn this into the adjacency matrix of an undirected graph. 

- Convert the numeric interaction score given in the matrix into edge indicators (0/1). Here is one way of doing it:
  - Set entries with absolute value greater than 0.2 as 1 (there are more principled ways to do this)
  - Ensure the diagonal of the matrix is 0
  - Set missing entries to 0

- Use only genes with entries in both rows and columns (that is, they were used both as queries and arrays). Use ORF ids to match rows and columns.
- Make the adjacency matrix symmetric. 

There is R code at the end of this page to perform all of these steps. Feel free to use this as you wish as part of your submission.

## Exercises

1. Compute and report the following quantities for the graph:  
  a) number of vertices,  
  b) number of edges,  
  c) average degree,  
  d) density
  
2. Compute the degree of each node and make the following two plots (these are similar to panels b. and c. of Image 2.4 of the Barabasi textbook):  
  a) a histogram of degrees,   
  b) a log-log plot of the degree distribution   
  

3. Compute the distance (length of shortest paths) between every pair of nodes. You have a few options here:

  - For each vertex in the network, use breadth-first search to compute distance to every other node
  - Use the Floyd-Warshall algorithm https://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
  - Section 10.3 of the Newman textbook discusses more options

4. Report:  
  a) the average distance  
  b) the network diameter  
  c) a plot of the distance distribution (simiar to 2.18 a in Barabasi textbook)
  
5. Is the network connected? How many components are there? Here you have a few options as well

  - Use BFS (Box 2.6 of Barabasi)
  - Use the eigenvalues of the graph Laplacian (Section 16.3.2 of Newman textbook)

6. Compute the clustering coefficient for each vertex in the network. Make a log-log plot of average clustering coefficient as a function of vertex degree (Image 2.18 b of Barabasi)

## What to turn in

Turn in a single `pdf` containing plots and answers to each exercise. Make sure to comment on how you processed data to get an adjacency matrix.
If you used the method described above state so, if you did something else, please describe what you did. Include all code used to answer this in your pdf. 
If using R, I recommend you use Rmarkdown to do this project. If using Python, I recommend you use a Jupyter notebook. 

You can work in self-organizing groups of at most 3 of your classmates. On ELMS, please have every member of the group submit the same pdf.


## R code to process data

```r
library(readr)
library(tidyverse)

# read data file
mat <- read_tsv("data/Data File S2. Raw genetic interaction datasets: Matrix format/SGA_NxN_clustered.cdt")

amat <- mat %>%
  # extract data rows and columns
  slice(-(1:5)) %>%
  select(starts_with("dma")) %>%
  
  # convert to numeric matrix
  type_convert() %>%
  as.matrix()
  
# turn into unweighted, undirected adjacency matrix  
amat <- 1 * (abs(amat) > 0.2)
  
# extract data about arrays (columns)
coldata <- mat %>%
  slice(1:5) %>%
  select(GID, starts_with("dma")) %>%
  slice(2) %>% select(-1) %>%
  gather(dma, orf, starts_with("dma"))
 
# extract ORF ids for queries (rows)
row_orf <- mat %>%
  slice(-(1:5)) %>%
  pull(ORF)

# match row and column ORF ids
m <- match(row_orf, coldata$orf)
rows_to_use <- !is.na(m)
cols_to_use <- m[rows_to_use]

# subset matrix into ORFs found in both rows and columns
amat <- amat[rows_to_use, cols_to_use]

# set diagonal and missing entries in matrix to 0
diag(amat) <- 0
amat[is.na(amat)] <- 0

# make the adjacency matrix symmetric
amat <- ceiling(0.5 * (amat + t(amat)))
```
