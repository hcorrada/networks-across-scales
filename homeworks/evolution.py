import sys
import numpy as np
import random
import re

version = int(sys.argv[1])

debug = False

ntimes = 10
nyears = 100
num_offspring = 5 
mutation_rate = 1. / 15
max_age = 10
max_population_size = 100

alphabet = list("ACGT")

first_string = 'TCGTACGGTATT'
exp = re.compile("T[CG]GT[ACGT]{4}T[AG][ACGT]T")
#exp = re.compile(first_string)

def reproduce(population, num_children, mutation_rate):
    new_population = []
    for age, string in population:
        new_population.append((age, string))

        for j in xrange(num_children):
            stringlist = list(string)
            for i in xrange(len(string)):
                if random.random() <= mutation_rate:
                    choice = random.choice(alphabet)
                    while choice == stringlist[i]:
                        choice = random.choice(alphabet)
                    stringlist[i] = choice
            new_population.append((0, "".join(stringlist)))
    return new_population

def remove_nonbinding(population):
    population = filter(lambda x: re.match(exp, x[1]) is not None, population)
    return population

def remove_elders(population, max_age):
    population = filter(lambda x: x[0] <= max_age, population)
    return population

def yearend(population, max_population_size):
    if len(population) > max_population_size:
        random.shuffle(population)
        population = population[:max_population_size]
    population = map(lambda x: (x[0] + 1, x[1]), population)
    return population

def calculate_entropy(population):
    k = len(population[0][1])
    t = len(population)
    
    res = 0.
    for i in xrange(k):
        counts = dict.fromkeys(alphabet, 0)
        for j in xrange(t):
            c = population[j][1][i]
            counts[c] += 1. / t
        res -= sum([p * (np.log2(p) if p > 0 else 0) for c, p in counts.iteritems()])
    return res

entropies = []
for i in xrange(ntimes):
	# age and sequence of primodial organism
    population = [(0, first_string)] 
    for j in xrange(nyears):
        if debug: print len(population),
		# five new individuals per
		# individual in the population, with random mutations at each
		# position with given probability
        population = reproduce(population, num_offspring, mutation_rate)
        if debug: print len(population),
	    # remove members of the population with non-binding sequence
		# (only in version 1 of the game)
        if version == 1:
            population = remove_nonbinding(population)
        if debug: print len(population),
	    # remove members of the population that are too old (10 years old)
        population = remove_elders(population, max_age)
        if debug: print len(population),
        # increase the age of each individual and keep at most 100
        #  individuals, choose randomly if populations is larger that 100
        population = yearend(population, max_population_size)
        if debug: print len(population)
    entropies.append(calculate_entropy(population))
print np.mean(entropies)
