---
date: 2017-02-07
title: "Project 1: Finding origin of DNA replication in a bacterial genome"
---

**Due: Monday Feb. 20, 2017**  
**Posted: Feb 9, 2017**   
**Last Update: Feb 9, 2017** 

# Part 1: Programming Exercises #

Submit your answers to Problems 1-3 in the
[Rosalind final submission page](http://rosalind.info/classes/401/).

**NOTE: THESE PROGRAMS NEED TO BE SUBMITTED TO ROSALIND BY Monday 02/20 at 11:59PM**

**Question 1. (10pts)** For each of the 3 solutions submitted provide a
runtime analysis of your solution.
For full credit, include a short description of your algorithm and
show how you derive your answer.

**Code Grading (45 pts)**

Each of the 3 Rosalind exercises will be evaluated based on correctness, efficiency and style.

1. Correctness: if you pass the Rosalind correctness check, you get full credit. Otherwise, your grade is determined by
	* Program running without error when called as `python <script.py> <input_file>` where `<input_file>` is in the format used for that problem in Rosalind.
	* Program reading input in the required format.  
	* Program printing output in correct format as described in Rosalind.  
	* Program implementing an algorithm that addresses the problem (i.e., you get more points if the algorithm is correct but you have a bug, e.g., with indexing)  

2. Efficiency: points are awarded for providing efficient solutions. For example, a linear algorithm will score better than a quadratic algorithm. For full credit, implement the efficient algorithms described in textbook and discussed in class.

3. Is the code in your submission clean and easy to read, with non-obvious statements
properly commented? Are functions used appropriately for clarity and organization?

# Part 2: Skew Diagrams #

*Campylobacter jejuni* is a well-known bacterial pathogen, recently
 found
 [to be associated with childhood dysentery](http://genomebiology.com/2014/15/6/R76)
 in developing countries. Here, you will apply your new skills in a
 preliminary analysis of the genome of this pathogen.

**Step 1. (5 pts)** Use the `BioPython` library to download the *Campylobacter
  jejuni* genome from NCBI:

```
from Bio import Entrez, SeqIO
Entrez.email = "me@example.com"
campy_id = "AL111168.1"

# open an url handle for query
handle = Entrez.efetch(db="nucleotide", id=campy_id, rettype="gb", retmode="text")

# read query result records
record = SeqIO.read(handle, "genbank")
handle.close()

# write sequence to fasta file (so you don't have to request again)
SeqIO.write(record, "campy.fa", "fasta")
```

Variable `record` is a `SeqIO` object containing the *Campylobacter
jejuni* genome. Check
[Chapter 2 of the Biopython tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
again to see how these objects work.
Also see Chapter 9 of the tutorial for more info on how to obtain DNA
sequences from NCBI.

**Step 2. (15 pts)** Write a script to plot the skew diagram for
  *Campylobacter jejuni*. One easy way of plotting is using pylab and IPython. If you start IPython as follows:

```
ipython
```

and `skew` is the skew vector you computed for *Campylobacter jejuni*
you can plot with

```
%pylab
plot(skew)
```

**Question 2 (10 pts)** Does this skew diagram look like the ones you've seen
so far in book and discussion? What do you think may account for any
differences? Where do you think replication origin for
  *Campylobacter jejuni* is located?

Save your skew diagram `pdf` or other image file and include in your submission (see below).

# Part 3: DnaA boxes #

Use your code to find most frequent words with mismatches (at most d=2 mismatches) and reverse complements to find candidate DnaA binding sequences in the oriC candidate region you found above (250bp on each side of position minimizing skew).

**Question 3 (10 pts)** Construct a table that shows for each k=3,4,5,6,7,8,9:   
(a) the number of distinct frequent words, and   
(b) the number of times each frequent word occurs in the genomic region.  

**Question 4 (5 pts)** Based on your result, did you find a reasonable candidate DnaA sequence you would provide to a biologist to test. If so, write down the candidate sequence and explain why you chose this particular candidate sequence. If not, explain why not.

# How to submit #

Prepare a writeup answering the four questions above and save as `pdf`. Submit any additional code (besides your solutions to the Rosalind exercises) used for this part of the homework. I highly recommend working with a `Jupyter` notebook, exporting the notebook as `pdf` and submitting that.

Submit to ELMS by Monday 02/20 11:59pm here: https://myelms.umd.edu/courses/1218381/assignments/4362579
