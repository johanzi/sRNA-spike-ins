# methods.sRNA.spike.in.design.step2.endogenous.miRNAs.py
#
# input: Fasta file containing miRNAs
# output: 1) folding structures calculated by RNAfold, 2) minimum free energies (MFEs) calculated by RNAfold

import sys, os

# Get first argument
inFile = sys.argv[1] 

# Execute RNAfold command on the input fasta file
cmd = 'cat ' + inFile + '| RNAfold -T 4 --noPS > ' + 'mature.miRNA.seqs_folded'
os.system(cmd)


# Retrieve mfes by parsing the output of RNAfold in mature.miRNA.seqs_folded
outFile = 'mature.miRNA.seqs_mfes'
ofh = open(outFile,'w')

inFile = 'mature.miRNA.seqs_folded'
ifh = open(inFile)
inLines = ifh.readlines()
ifh.close()

mfeList = []

for line in inLines:
    if line.startswith('.') or line.startswith('(') or line.startswith(')'):
        line = line.strip()
        mfe = float(line[-7:-1])
        mfeList.append(mfe)

for mfe in mfeList:
    print >>ofh,mfe

ofh.close()    
