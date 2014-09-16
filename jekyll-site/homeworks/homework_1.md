---
layout: default
title: Homework 1
---

**Due Monday Sept. 22 in class**

## Programming Questions ##

Submit your answers to Problems 1-4 in the [Rosalind final submission page](http://rosalind.info/classes/156/).

**Question 1.** For each of the 4 solutions posted, provide a runtime analysis of the solution you provided.
For full credit, include a description of how you derive your answer.

## Finding origin of replication ##

**Step 1.** Use the `BioPython` library to download the XX genome from NCBI:

{% highlight python %}
import Bio
{% endhighlight %}

**Step 2.** Use your code for the `Minimum Skew Problem` to find a candidate replication origin location in XX.

**Step 3.** Write a program to plot the skew diagram for genome XX. You can use matplotlib for this:

{% highlight python %}
import matplotlib
{% endhighlight %}

## Finding candidate DnaA binding sites ##

**Step 1.** Extract the 500bp DNA sequence from the XX genome centered around the point of replication origin you derived in the previous section.

**Step 2.** Using your code, find a candidate 9-mer that could be the DnaA binding site sequence in this region.

**Question 2.** How confident are you this is a good candidate binding site?

## How to submit ##

Prepare a one page `pdf` document answering the two questions above, and including the skew diagram for XX. Include
a listing of any additional code (besides your solutions to the Rosalind exercises) used for this part of the homework.
Print out and handin before class on **Monday Sept. 22**.


