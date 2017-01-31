## Key outcomes

- The idea that we perform _sequence analysis_ to find where a biochemical event (binding of DNA polymerase) occurs in a genome
- One aspect of _sequence analysis_ is analyzing the frequency at which substrings occur in a genome
- How we can setup _probability_ models to analyze frequency, so we can make statements about substrings that occur _surprisingly often_
- One way to setup a probability model is to count the number of strings of length $k$ that satisfy a property out of the set of all possible strings of length $k$
- The ability to reason about the time complexity of various ways of counting $k$-mers

## Things to discuss:

- DNA replication and the OriC region
- The FrequentWords implementation in page 9. What is the time complexity of this algorithm?
- The ways of speeding up frequent words (pages 39-44). Don't spend to much time on the 'PatternToNumber' function. That's just hashing. Time complexity analysis of each of these
- Pattern probabilities pg 52-56

## Group exercises

Make them get into groups of no more than 4, discuss each problem for 3 minutes or so, then go over it with them

1) Analyze the time complexity of Pg. 9 algorithm
2) Analyze the time complexity of FrequentWordsWithSorting
3) What is the probabilty of generating a palindromic (e.g., ATCGAAGCTA) DNA string of length $k$
