import itertools
import numpy as np
import re

# User defined for any length of file, or any number of mers
# This doesn't need to be changed, but could, in case we wanted
# to look at say the first 2000 sequences for 5-mers.
num_seqs = 1000
n_mer = 4
file = 'subset.fastq'

# 'merdict' will be used as a collection of unique k-mer keys with quality
# for each instance appended to its value as a list. For example,
# merdict = {'GATC' : [33, 31, 12, 25, 16], 'GATG' : [22, 11, 5, 49]}
merdict = {}
def compile_merdict(hits, mer_quality):
    for i in hits:
        if i in merdict:
            merdict[i].append(mer_quality)  # k-mer present? Add the new qscore
        else:
            merdict[i] = [mer_quality]      # k-mer not present? Add it to the dict
    return merdict

# Creates 4^k-mer combinations as a filter
all_mers = []
nuc = ['G','A','T','C']
for i in range(n_mer, n_mer+1):
    all_possible = itertools.product(*itertools.repeat(nuc, i))
    for j in all_possible:
        all_mers.append("".join(list(j)))

# Quality score list
p33 = {'!':0,'\"':1,'#':2,'$':3,'%':4,'&':5,'\'':6,'(':7,')':8,'*':9,'+':10,
       ',':11,'-':12,'.':13,'/':14,'0':15,'1':16,'2':17,'3':18,'4':19,'5':20,
       '6':21,'7':22,'8':23,'9':24,':':25,';':26,'<':27,'=':28,'>':29,'?':30,
       '@':31,'A':32,'B':33,'C':34,'D':35,'E':36,'F':37,'G':38,'H':39,'I':40,
       'J':41,'K':42}

# Straightforward, use the corresponding index to find the average value
def compute_qual(mer_quality):
    scores = []
    for i in mer_quality:
        scores.append(p33[i])
    return np.mean(scores)

# Read the file, replace Gs AND Ns with a wildcard '.' for regex matching.
with open(file, 'r+') as fq:
    n = 0
    while n < num_seqs*4:
        n += 4
        
        head = fq.readline().strip('\n')
        seq = fq.readline().strip('\n').replace('G','.')
        seq = seq.replace('N','.')
        sep = fq.readline().strip('\n')
        qual = fq.readline().strip('\n')

# Begin iterating through the k-mers, matching them against 'all_mers' list
        for i in range(len(seq)):
            if i+n_mer <= len(seq):
                mer = seq[i:i+n_mer]
                query = re.compile(mer)
                hits = list(filter(query.match, all_mers))  # Creates list of all hits
                qmer = qual[i:i+n_mer]                      # Quality string of mer
                mer_quality = compute_qual(qmer)
                compile_merdict(hits, mer_quality)

# Build a dictionary with unique k-mer keys and count occurences
library = {}
for k,v in merdict.items():
    library[k] = len(v)

# Determine max, along with Q1, median, and Q3
mostcommon = max(library.keys(), key=(lambda k: library[k]))
qmq = np.percentile(merdict[max(library)], [25,50,75])

# Report
with open('results.txt', 'w') as res:
    print(str(n_mer) + "-mer: " + mostcommon, file=res)
    print("Count: " + str(library[mostcommon]), file=res)
    print("Quality lower quartile: " + str(qmq[0]), file=res)
    print("Quality median: " + str(qmq[1]), file=res)
    print("Quality upper quartile: " + str(qmq[2]), file=res)
