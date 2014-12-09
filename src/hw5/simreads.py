from kmer_composition import kmer_composition
import numpy as np
import math

alphabet = list('ACGT')

def simulate_reads(text, k, coverage, error_rate):
    kmers = kmer_composition(k, text)
    nreads = (len(kmers) * coverage) / k

    reads = list(np.random.choice(np.array(kmers), size=nreads, replace=True))

    errors_per_read = math.floor(error_rate * k)
    nerrors = np.random.poisson(errors_per_read, len(reads))
    
    for i in xrange(len(reads)):
        if nerrors[i] == 0:
            continue
        
        read = list(reads[i])
        positions = np.random.choice(len(read), size=nerrors[i])
        for p in positions:
            read[p] = np.random.choice(alphabet)
        reads[i] = ''.join(read)
        
    return reads

