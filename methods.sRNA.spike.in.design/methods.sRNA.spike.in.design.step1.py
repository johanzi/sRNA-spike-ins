#methods.sRNA.spike.in.design.step1.py

# Input = fasta file containing miRNAs (it can contain either T or U, it works for both)
# Output = random.fa (the Us are turned in Ts to allow proper mapping in bowtie)

import random, sys

# Provide as argument the fasta file containing the miRNAs
ref_fasta = sys.argv[1]

fasta =  open(ref_fasta,"r")
fasta = fasta.read()
fasta = fasta.split("\n")

list_seq = []

for seq in fasta:
    if seq[:1] != ">":
        seq = seq.strip()
        list_seq.append(seq)

#Remove empty line if there is
list_seq = [_f for _f in list_seq if _f]

# Create an empty dictionary
freqDict = {}

# Create a dictionary of 21 elements (0 to 20), each element has a key (position) and a value tuple of 4 values (frequency of each
# nucleotides
for i in range(0,21):
    freqDict[i] = (0,0,0,0)

# Retrieve only the first column of the input file (miRNA sequence) and assess the frequency of each
# nucleotide for each position and each miRNA
# Fill in the tuple for each element of the dictionary so that every position has 
# 4 frequency values
for seq in list_seq:
    for i in range(0,21):
        if len(seq) > i: 
            if seq[i] == 'A':
                freqDict[i] = (freqDict[i][0] + 1,freqDict[i][1],freqDict[i][2],freqDict[i][3])
            if seq[i] == 'G':
                freqDict[i] = (freqDict[i][0],freqDict[i][1] + 1,freqDict[i][2],freqDict[i][3])
            if seq[i] == 'C':
                freqDict[i] = (freqDict[i][0],freqDict[i][1],freqDict[i][2] + 1,freqDict[i][3])
            if seq[i] == 'T' or seq[i] == 'U':
                freqDict[i] = (freqDict[i][0],freqDict[i][1],freqDict[i][2],freqDict[i][3] + 1)


# Create a new dictionary that will contain the chain of each nucleotides based on their frequency
freqDict2 = {}

# Then add a number of CGUA equivalent to the frequencies values given in freqDict
# This will create a dictionary of 21 elements, each element containing x times each nucleotide 
for pos in freqDict.keys():
    baseList = []
    A_count = freqDict[pos][0]
    G_count = freqDict[pos][1]
    C_count = freqDict[pos][2]
    U_count = freqDict[pos][3]

    for i in range(0,A_count):
        baseList.append('A')
    for i in range(0,G_count):
        baseList.append('G')
    for i in range(0,C_count):
        baseList.append('C')
    for i in range(0,U_count):
        baseList.append('T')
    
    freqDict2[pos] = (baseList)

# Redirect stdout to random.fa
outFile = 'random.fa'
ofh = open(outFile,'w')

# Look for each elements in the dictionary freqDict2 and choose randomly one nucleotide
# in each element, and do that for the 21 elements so that a sequence of 21 nucleotides
# is built. The process is reiterated 1000 times and formated in fasta format
# Only the nucleotide from 5 to 17 are displayed 
for i in range(0,1000):
    randomSeq = ''
    for pos in freqDict2.keys():
        randomSeq = randomSeq + random.choice(freqDict2[pos])
    print >>ofh,'>' + str(i)
    print >>ofh,randomSeq[4:17]

