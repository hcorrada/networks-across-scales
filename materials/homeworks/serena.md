---
title: Win probabilities with dynamic programming
author: CMSC423 
fontfamily: utopia
geometry: margin=1in
---

Name(s):   
UID(s):  

## Serena Williams Win Probability

Could Serena Williams be a *more* dominant player if she played best-out-of-five-sets instead of best-out-of-three?

[Five-thirty-eight thinks so](http://fivethirtyeight.com/datalab/serena-williams-grand-slam-us-open-best-of-five-sets/)

How can we use dynamic programming to address this question?

Let's use this notation:

- $p$: probability that Serena wins a set over another player (say, Olympic gold medalist Monica Puig)   
- $q=(1-p)$: probability that the other player, Monica Puig, wins a set

Suppose Serena and Monica play a game until one player wins $n$ sets.

Denote $P(i,j)$ as the probability that Serena wins the game when:  
- Serena needs to win $i$ more sets, and  
- Monica needs to win $j$ more sets
 
Using this notation, $P(n,n)$ is the probability Serena wins the game. For example,
the probability Serena wins a best-out-of-five game is $P(3,3)$, and a best-of-three game is $P(2,2)$.

Your goal is to write a dynamic programming algorithm that computes this probability.

Step 1) Write a recursive solution to compute the probability:

a) What are base cases $P(0,j) \, \forall j=0 \dots n$ and $P(i,0) \, \forall i=0 \dots n$.  
b) Write an expression for $P(i,j)$ in terms of $p$, $q$, $P(i-1,j)$ and $P(i,j-1)$. Why don't we use $P(i-1,j-1)$ in this expression?

Step 2) Write pseudo code for a dynamic programming algorithm that computes $P(n,n)$.
