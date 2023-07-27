def readFastq(filename):
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # skip name line
            seq = fh.readline().rstrip()  # read base sequence
            fh.readline()  # skip placeholder line
            qual = fh.readline().rstrip() # base quality line
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

def readGenome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header line with genome information
            if not line[0] == '>':
                genome += line.rstrip()
    return genome

def reverseComplement(s):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t

def naive(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

def naive_with_rc(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                match = False
                break
        if match and match not in occurrences:
            occurrences.append(i)  # all chars matched; record

    if p != reverseComplement(p):
        p = reverseComplement(p)
        for i in range(len(t) - len(p) + 1):  # loop over alignments
            match = True
            for j in range(len(p)):  # loop over characters
                if t[i+j] != p[j]:  # compare characters
                    match = False
                    break
            if match and match not in occurrences:
                occurrences.append(i)  # all chars matched; record
    return occurrences


import urllib.request
# urllib.request.urlretrieve("http://d396qusza40orc.cloudfront.net/ads1/data/phix.fa", "phix.fa")
#
phix_genome = readGenome('phix.fa')
occurrences = naive_with_rc('ATTA', phix_genome)

# print('offset of leftmost occurrence: %d' % min(occurrences))
#
# print('# occurrences: %d' % len(occurrences))

# urllib.request.urlretrieve("https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/lambda_virus.fa", "lambda_virus.fa")

lambdavirus = readGenome('lambda_virus.fa')
# occurrences = naive_with_rc('AGTCGA', lambdavirus)
# print(occurrences)
#
# print('offset of leftmost occurrence: %d' % min(occurrences))
#
# print('# occurrences: %d' % len(occurrences))

def naive_2mm(p, t):
    occurrences = []
    #check mismatches
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        counter = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
               counter += 1
               if counter == 3:
                   match = False
                   break
        if match:
            occurrences.append(i)  # all chars matched; record

    return occurrences


urllib.request.urlretrieve("https://d28rh4a8wq0iu5.cloudfront.net/ads1/data/ERR037900_1.first1000.fastq", "human.fastq")

read, qual = readFastq('human.fastq')

print(qual[0:5])



