---
layout: default
title: Project 4
---

**Due Thrusday Dec. 11**

You will implement approximate string matching using a
seed-and-check strategy based on exact matching using the
Burrows-Wheeler transform. We are providing starter code
to get you going. Full instructions and code are here:

[https://gitlab.umiacs.umd.edu/hcorrada/cmsc423_project4](https://gitlab.umiacs.umd.edu/hcorrada/cmsc423_project4)

I highly, highly recommend that you use [`git`](https://git-scm.com/book/en/v1/Getting-Started)

To get code you will be using for this project using git:

{% highlight bash %}
git clone https://gitlab.umiacs.umd.edu/hcorrada/cmsc423_project4.git
{% endhighlight %}

This will create a `project4` directory with all code and data required for this project.

## Programming Questions ##

Check and submit your solution to Final Rosalind Problem 8
[Rosalind final submission page](http://rosalind.info/classes/233/).

Instructions on how to prepare your solution are given in
[project code repository](https://gitlab.umiacs.umd.edu/hcorrada/cmsc423_project4).

## SNP Finding ##

Use your code  to analyze sequences for the [neuraminidase Influenza gene (NA)](http://en.wikipedia.org/wiki/Influenza_neuraminidase) of
two H1N1 Human Influenza strains.

There are two FASTA files in directory `data` of the project repository. The `reference.fa` file contains a reference sequence you will use as target. The `reads.fa` file contains 50 bp reads you will use as queries. Use your approximate matching BWT solution to align reads to reference allowing up to 3 mismatches.

**Question 1**. A specific mutation (H275Y) in the NA gene confers resistance to [Oseltamivir](http://en.wikipedia.org/wiki/Neuraminidase_inhibitors), making
the drug [less effective](http://www.ncbi.nlm.nih.gov/pubmed/22837199). A note about nomenclature: the code H275Y
encodes a substitution in position 275 (1-based indexing), changing aminoacid H (Histidine) to Y (Tyrosine).

Do the reads correspond to a gene sequence with this mutation? How can you tell?

###Notes:

1. As we've done so far use the `Biopython` to read and import `fasta` files.
2. The following `pileup.py` script in the project repository can help you answer this question. See
When used as follows you can get a list of mismatches found in reads, and the number of times each kind of mismatch is observed in each position.

{% highlight python %}
from pileup import PileUp
from approximate_matcher import ApproximateMatcher

# initialize object
am = ApproximateMatcher(reference)
pileup = PileUp(reference)

for read in reads:
	# find matching positions for a given read
	# assumes positions is a list (even if only a single match is found)
	# with matching positions
	positions = am.get_matches(read, d)
	if len(positions) > 0:
		# add to pileup object
		pileup.insert(positions, read)

# prints out mismatching positions
# output is:
# (<position>, <reference_character>, [(<variant_character>,
# <num_times_aligned>)])
# argument filters mismatch by frequency in which variant character
# is observe, e.g., .01 means variant character has to be seen at least
# once for every 10 aligned nucleotides
p.print_mismatches(.01)
{% endhighlight %}

where positions are indices in reference where a match for `read` was found.

## How to submit ##

On ELMS you will submit two things

(1) A `diff` of your solution and the original code posted by us. Using git you can do this as follows
once you have committed all your code changes in directory `approximate_matcher`. Note that you do not need to change code outside of this directory to solve the Rosalind problem:

{% highlight bash %}
git diff origin/master approximate_matcher > project4_bwt_code.patch
{% endhighlight %}

Submit file `project4_bwt_code.patch` on ELMs

(2) An IPython notebook exported as pdf or html containing your answer to Question 1 above along with code used to answer it.
You have to use your approximate matching code from part I to match reads to the reference, otherwise, you are free to use Biopython as needed to answer this question. E.g., to do DNA->aminoacid translation.
