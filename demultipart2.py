#!/usr/bin/env python3
#SBATCH --partition=long       ### Partition (like a queue in PBS)
#SBATCH --job-name=RRPS7          ### Job Name
#SBATCH --time=1-20:01:00       ### Wall clock 0ime limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=28     ### Number of tasks to be launched per Node
#SBATCH --mail-user=rarichardson92@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL

# Don't forget to load modules in bash, easybuild, prl, python/3.6.0 before running code!

import argparse

def getarguments():
    parser=argparse.ArgumentParser(description = "Demultiplexes gzipped index and read outputs. requires fastq files and a file with tab seperated indexes. Outputs files based on read index. Requires carefully choosen indexes for single base error correction.")
    parser.add_argument("-r1", help = "Defines name and path of read1 file to use in program. Required, must be a string.", required = True, type = str)
    parser.add_argument("-r2", help = "Defines name and path of read2 file to use in program. Required, must be a string.", required = True, type = str)
    parser.add_argument("-i1", help = "Defines name and path of index1 file to use in program. Required, must be a string.", required = True, type = str)
    parser.add_argument("-i2", help = "Defines name and path of index2 file to use in program. Required, must be a string.", required = True, type = str)
    parser.add_argument("-a","--allindex", help = "Defines name and path of tab seperated file with all indexes to use in program. Required, must be a string.", required = True, type = str)
    parser.add_argument("-n", help = "Number of records in files (must all be equal, will only process up to that number). Required, must be an integer.", required = True, type = int)
    parser.add_argument("-qi", help = "Mean quality score cutoff for indexes. Optional. Defaults to 30, must be an integer.", required = False, type = int, default=30)
    parser.add_argument("-qr", help = "Mean quality score cutoff for reads. Optional. Defaults to 30, must be an integer.", required = False, type = int, default=30)
    parser.add_argument("-m", help = "Index mismatch tolerance. Optional. Defaults to 1, must be an integer.", required = False, type = int, default=1)
    return parser.parse_args()

args=getarguments()

index1=str(args.i1)
index2=str(args.i2)
read1=str(args.r1)
read2=str(args.r2)
allindex=str(args.allindex)
number=int(args.n)
indexq=int(args.qi)
readq=int(args.qr)
mis=int(args.m)

#Assumes input gzip files
#Initializes user defined variables

#set for number of records in fastq file

import gzip

def reverse(index):
    """Generates the complimentary sequence of input"""
    index=index.strip('\n')
    comp=""
    #holds reverse sequence
    for i in reversed(range(len(index))):
        if index[i] == "A":
            comp=comp+"T"
        elif index[i] == "T":
            comp=comp+"A"
        elif index[i] == "C":
            comp=comp+"G"
        elif index[i] == "G":
            comp=comp+"C"
        elif index[i] == "N":
            comp=comp+"N"
    return(comp)

def makedict(allindex):
    """Generates a dictionary with indexes as keys and the reverse compliments as values"""
    indexdict={}
    #holds dictionary
    with open(allindex, "r") as ai:
        indexstring=ai.readline()
        indexstring=indexstring.strip('\n')
        list=indexstring.split('\t')
        for index in list:
            indexdict[index]=reverse(index)
    return(indexdict)

def convert_phred(letter):
    """Converts a single character into a phred score. Will not tolerate letters outside of Allchar,
    which are the letters used for +33 quality scores."""
    Allchar="!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJK"
    if letter in Allchar:
        return Allchar.find(letter)
    else:
        raise ValueError("Invalid entry for +33 phred scores.")

def checkqual(line):
    """Checks mean quality of input line."""
    line=line.strip('\n')
    line_qscores=0
    #holds total quality scores for a line
    letternum=0
    #counts letters in a line
    for letter in line:
        line_qscores+=convert_phred(letter)
        letternum += 1
    return (line_qscores/letternum)

def errorcorrect(readindex,indexdict, mis):
    """Corrects errors of base pairs in index reads. Prevents correction where the mismatch is ambiguous to which index it should be corrected to. Returns forward index or 0 for noncorrection/no match. Assumes reverse compliments are more than mismatch base difference."""
    onecorrect=0
    passingdict={}
    Doubleflag=False
    least=mis+1
    #hold index of corrected index or 0 for multiples and no matches
    readindex=readindex.strip('\n')
    length=len(readindex)
    #determines length of index, assumes equal lengths for index input and index in dictionary
    for index in indexdict.keys():
        mismatch=0
        for i in range(length-1):
            if index[i]!=readindex[i]:
                mismatch+=1
        if mismatch<least:
            passingdict[index]=mismatch

        readindex=reverse(readindex)
        #Reverse compliments read index, repeats above process
        for index in indexdict.keys():
            mismatch=0
            for i in range(length-1):
                if index[i]!=readindex[i]:
                    mismatch+=1
            if mismatch<least:
                if index not in passingdict:
                    passingdict[index]=mismatch
                elif mismatch<passingdict[index]:
                    passingdict[index]=mismatch

        for item in passingdict:
            if passingdict[item]==least:
                Doubleflag=True
            if passingdict[item]<least:
                least=passingdict[item]
                Doubleflag=False
                onecorrect=item
        if Doubleflag==True:
            onecorrect=0
    return(onecorrect)


indexdict={}
indexdict=makedict(allindex)
#makes index:reversecompliment dictionary
indexhop={}
indexgood={}
for index in indexdict.keys():
    indexhop[index]=0
    indexgood[index]=0
#makes dictionaries for good and bad index tallies; index hopping categorized by index1
# index:###


RN=0
#record number
with gzip.open(read1, "rt") as r1:
    with gzip.open(read2, "rt") as r2:
        with gzip.open(index1, "rt") as i1:
            with gzip.open(index2, "rt") as i2:
                with gzip.open("unknown.r1.gz", "at") as u1:
                    with gzip.open("unknown.r2.gz", "at") as u2:
                        while(RN < number):
                            r1header=r1.readline()
                            r1seq=r1.readline()
                            r1plus=r1.readline()
                            r1qual=r1.readline()
                            r2header=r2.readline()
                            r2seq=r2.readline()
                            r2plus=r2.readline()
                            r2qual=r2.readline()
                            i1header=i1.readline()
                            i1seq=i1.readline()
                            i1plus=i1.readline()
                            i1qual=i1.readline()
                            i2header=i2.readline()
                            i2seq=i2.readline()
                            i2plus=i2.readline()
                            i2qual=i2.readline()
                            RN+=1
                            #opens all files and saves one record as variables, increments RN by 1

                            indexcheck=0
                            #Used as a variable to check indexes with dictionary
                            checknumi1=checkqual(i1qual)
                            checknumi2=checkqual(i2qual)
                            checknumr1=checkqual(r1qual)
                            checknumr2=checkqual(r2qual)
                            #Checks quality for all sequence lines

                            if checknumi1>=indexq and checknumi2>=indexq and checknumr1>=readq and checknumr2>=readq:
                                #Quality check - passing to continue
                                if i1seq.strip('\n') not in indexdict.keys() and i1seq.strip('\n') not in indexdict.values():
                                    #Check if NOT in dictionary
                                    indexcheck=errorcorrect(i1seq, indexdict, mis)
                                    #Attempts error correction, if successful, forward index assigned for comparison (otherwise indexcheck==0)
                                else:
                                    if i1seq.strip('\n') in indexdict.keys():
                                        indexcheck=i1seq.strip('\n')
                                    else:
                                        indexcheck=reverse(i1seq.strip('\n'))
                                        #For above if, else: if in dictionary keys or values, assigns forward index as indexcheck

                                if indexcheck!=0:
                                    #Where not no match (AKA match found)
                                    if indexcheck!=i2seq.strip('\n') and reverse(indexcheck)!=i2seq.strip('\n') and indexcheck!=errorcorrect(i2seq.strip('\n'), indexdict, mis):
                                        #check if saved index or reverse is not equal to index 2, as well as corrected index2
                                        u1.write(r1header)
                                        u1.write(r1seq)
                                        u1.write(r1plus)
                                        u1.write(r1qual)
                                        u2.write(r2header)
                                        u2.write(r2seq)
                                        u2.write(r2plus)
                                        u2.write(r2qual)
                                        #Writes to unknown files

                                        if i2seq.strip('\n') in indexdict or i2seq.strip('\n') in indexdict.values() or errorcorrect(i2seq.strip('\n'), indexdict, mis) in indexdict:
                                            indexhop[indexcheck]+=1
                                            #Discriminates high QC sequence errors from index hopping
                                    else:
                                        #index1 == index2
                                        with gzip.open(indexcheck+".r1.gz", "at") as out1:
                                            with gzip.open(indexcheck+".r2.gz", "at") as out2:
                                                out1.write(r1header)
                                                out1.write(r1seq)
                                                out1.write(r1plus)
                                                out1.write(r1qual)
                                                out2.write(r2header)
                                                out2.write(r2seq)
                                                out2.write(r2plus)
                                                out2.write(r2qual)
                                        indexgood[indexcheck]+=1
                                else:
                                #No index could be found that matches the read (indexcheck==0)
                                    u1.write(r1header)
                                    u1.write(r1seq)
                                    u1.write(r1plus)
                                    u1.write(r1qual)
                                    u2.write(r2header)
                                    u2.write(r2seq)
                                    u2.write(r2plus)
                                    u2.write(r2qual)
                                    #Writes to unknown files
                            else:
                                #Failed Quality check, discard lines
                                u1.write(r1header)
                                u1.write(r1seq)
                                u1.write(r1plus)
                                u1.write(r1qual)
                                u2.write(r2header)
                                u2.write(r2seq)
                                u2.write(r2plus)
                                u2.write(r2qual)
                                #Writes to unknown files

goodsum=0
hopsum=0

print("")
print("Parameters:")
print("Demultiplexing run on the following files:", read1, read2, index1, index2)
print("Number of records set to "+str(number))
print("Index mismatch tolerance set to "+str(mis))
#Can probably be changed to be manipulated later on
print("Index quality score mean cutoff set to "+str(indexq)+" and read quality score cutoff set to "+str(readq))
print("Index file used for demultiplexing:", allindex)
print("Full set of input indexes as follows:")
print(indexdict.keys())
print("")
print("Please refer arguments in the \"--help\" menu to adjust parameters.")
print("")
print("Number of good reads per index, percentage of all reads:")
for index in indexgood.keys():
    if indexgood[index] > 0:
        print(index, indexgood[index], str((indexgood[index]/number)*100)+"%")
        goodsum+=indexgood[index]
print("Index hopping per index (sorted by index 1), percentage of all reads:")
for index in indexhop.keys():
    if indexhop[index] > 0:
        print(index, indexhop[index], str((indexhop[index]/number)*100)+"%")
        hopsum+=indexhop[index]
print("")
print("Total number of read sets: "+str(number))
print("Total percentage of high quality read sets: "+str(((goodsum+hopsum)/number)*100)+"%")
print("Total number of high quality read sets without index hopping: "+str(goodsum))
print("Total number of high quality read sets with index hopping: "+str(hopsum))
print("Percentage of high quality read sets with index hopping: "+str((hopsum/(goodsum+hopsum))*100)+"%")
