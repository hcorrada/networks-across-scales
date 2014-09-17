---
layout: default
title: Homework 1
---

**Due Monday Sept. 22 in class**

## Programming Questions ##

Submit your answers to Problems 1-4 in the
[Rosalind final submission page](http://rosalind.info/classes/156/).
**THESE ARE DUE ON FRIDAY SEPT. 19**

**Question 1.** For each of the 4 probems you posted, provide a
runtime analysis of your solution.
For full credit, include a short description of your algorithm and
show how you derive your answer.

## Skew Diagrams ##

*Campylobacter jejuni* is a well-known bacterial pathogen, recently
 found
 [to be associated with childhood dysentery](http://genomebiology.com/2014/15/6/R76)
 in developing countries. Here, you will apply your new skills in a
 preliminary analysis of the genome of this pathogen.
 
**Step 1.** Use the `BioPython` library to download the *Campylobacter
  jejuni* genome from NCBI:

{% highlight python %}
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
{% endhighlight %}

Variable `record` is a `SeqIO` object containing the *Campylobacter
jejuni* genome. Check
[Chapter 2 of the Biopython tutorial](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
again to see how these objects work.
Also see Chapter 9 of the tutorial for moew info on how to obtain DNA
sequences from NCBI.

**Step 2.** Write a script to plot the skew diagram for
  *Campylobacter jejuni*. If  you start IPython as follows:

{% highlight bash %}
ipython --pylab=auto
{% endhighlight %}

and `skew` is the skew vector you computed for *Campylobacter jejuni*
you can plot with

{% highlight python %}
plot(skew)
{% endhighlight %}

**Question 2** Does this skew diagram look like the ones you've seen
so far in book and discussion? What do you think may account for any
differences? Where do you think is the replication origin for
  *Campylobacter jejuni*?
  
## How to submit ##

Prepare a one page `pdf` document answering the two questions above,
and including the skew diagram for *Campylobacter jejuni*. Include
a listing of any additional code (besides your solutions to the Rosalind exercises) used for this part of the homework.
Print out and handin before class on **Monday Sept. 22**.


