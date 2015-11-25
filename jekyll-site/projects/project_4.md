---
layout: default
title: Homework 5
---

**Due Friday Dec. 12**

## Programming Questions ##

Submit your answer to Problems 13 and 14 in the
[Rosalind final submission page](http://rosalind.info/classes/156/).

## SNP Finding (this is a bonus, extra 10% on final HW grade) ##

Use your code for multiple, approximate pattern matching (Rosalind final submission # 14) to analyze sequences for the [neuraminidase Influenza gene (NA)](http://en.wikipedia.org/wiki/Influenza_neuraminidase) of
two H1N1 Human Influenza strains. 

This `fasta` file contains a reference sequence you will use as target:

[http://cbcb.umd.edu/~hcorrada/CMSC423/downloads/reference.fa](http://cbcb.umd.edu/~hcorrada/CMSC423/downloads/reference.fa)

This `fasta` file contains 50 bp reads you will use as queries:

[http://cbcb.umd.edu/~hcorrada/CMSC423/downloads/reads.fa](http://cbcb.umd.edu/~hcorrada/CMSC423/downloads/reads.fa)

Use your approximate matching BWT solution to align reads to reference allowing up to 3 mismatches.

**Question 1**. A specific mutation (H275Y) in the NA gene confers resistance to [Oseltamivir](http://en.wikipedia.org/wiki/Neuraminidase_inhibitors), making
the drug [less effective](http://www.ncbi.nlm.nih.gov/pubmed/22837199). A note about nomenclature: the code H275Y
encodes a substitution in position 275, changing H (Histidine) to Y (Tyrosine).

Do the reads correspond to a sequence harboring this mutation? How can you tell?

###Notes:

1. As we've done so far use the `Biopython` to read and import `fasta` files.
2. The following script can help you answer this question:

[http://cbcb.umd.edu/~hcorrada/CMSC423/downloads/pileup.py](http://cbcb.umd.edu/~hcorrada/CMSC423/downloads/pileup.py)

When used as follows you can get a list of mismatches found in reads, and the number of times each kind of mismatch is observed in each position.

{% highlight python %}
from pileup import PileUp

# initialize object
p = PileUp(reference)
for read in reads:
	# find matching positions for a given read
	# assumes positions is a list (even if only a single match is found)
	# with matching positions
	positions = find_approximate_matches(read, bwt, d)
	
	# add to pileup object
	p.insert(positions, read)

# prints out mismatching positions
# output is: 
# (<position>, <reference_character>, [(<variant_character>, 
# <num_times_aligned>)])
# argument filters mismatch by frequency in which variant character
# is observe, e.g., .01 means variant character has to be seen at least 
# once for every algined base
p.print_mismatches(.01)
{% endhighlight %}

where positions are indices in text where a match for `read` was found.

## How to submit ##

Submit your answer to this question, along with any additional code you used to answer it by email
to **both** `hcorrada@cs.umd.edu` and `wikum@cs.umd.edu` by 5pm, December 12.


