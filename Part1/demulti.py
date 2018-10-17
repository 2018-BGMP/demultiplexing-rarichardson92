#!/usr/bin/env python3
#SBATCH --partition=short       ### Partition (like a queue in PBS)
#SBATCH --job-name=RRPS7      	### Job Name
#SBATCH --time=0-20:01:00       ### Wall clock 0ime limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=14     ### Number of tasks to be launched per Node
#SBATCH --mail-user=rarichardson92@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

# Don't forget to load modules in bash, easybuild, prl, python/3.6.0 before running code!
#Install matplotlib for import and generation of plot!!

import argparse

def getarguments():
	parser=argparse.ArgumentParser(description = "Creates a distribution of mean quality scores across base positions in a zipped fastq file.")
	parser.add_argument("-f", "--file", help = "Defines name and path of file to use in program. Required, must be a string.", required = True, type = str)
	parser.add_argument("-n", "--numreads", help = "Defines number of reads in the file. Required, must be an integer.", required = True, type = int)
	parser.add_argument("-o", "--out", help = "Defines name and path of output plot file to use in program. Required, must be a string.", required = True, type = str)
	return parser.parse_args()

args=getarguments()
numreads=str(args.numreads)
file=str(args.file)
out=str(args.out)

#Initializes user defined variables, table for table of corresponding genes and protiens and data for data file of all protiens

#file = "./testfile.2.fastq.gz"
#defines input file

#out = "outthing"
#defines outfile name

#numreads = 10
#numreads = 363,246,735
#set for number of reads in fastq file

import numpy as np
import gzip

def convert_phred(letter):
    """Converts a single character into a phred score. Will not tolerate letters outside of Allchar,
    which are the letters used for +33 quality scores."""
    Allchar="!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJK"
    if letter in Allchar:
        return Allchar.find(letter)
    else:
        raise ValueError("Invalid entry for +33 phred scores.")

def getlength(fastq):
    """Finds the length of the first sequence in a fastq file and returns length"""
    with gzip.open(fastq, "rt") as fh:
        LN=0
    #LN = line number
        for line in fh:
            line=line.strip('\n')
            LN += 1
            #LN = line number, adds 1 per line
            if LN == 2:
                return len(line)


seqlen = getlength(file)
#Returns length for first sequence, assumed all lines in file are at that length
all_qscores=np.zeros((seqlen), dtype=int)
#generates a np array thats as large as the expected file inputs for every quality score for each position in the sequence

with gzip.open(file, "rt") as fh:
    LN=0
    #LN = line number
    for line in fh:
        line=line.strip('\n')
        LN += 1
        #LN = line number, adds 1 per line
        if LN%4 == 0:
            #Specifies quality score line
            letternum=0
            #position of base in quality score line, used to determine which index it should be stored in for all_qscores
            for letter in line:
                all_qscores[letternum]+=convert_phred(letter)
                letternum += 1
#iterates through all quality scores in the file and adds each converted phred quality score to the appropriate base position in the all_qscores array
    if LN%1000000 ==0:
        print("Working on line"+str(LN))

mean_scores=all_qscores/(LN/4)
#Assigns statistical lists to appropriate statistics


import matplotlib
import matplotlib.pyplot as plt

plt.bar(range(seqlen), mean_scores)
plt.title('Average quality scores at each base position')
plt.xlabel('Base position')
plt.ylabel('Quality score')
plt.savefig(out+".qscore.distribution.png")
#Plots means of each position
