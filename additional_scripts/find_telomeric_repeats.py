#!/usr/bin/env python

#Find telomeric repeats

import argparse
import sys, re
import glob
import string
from Bio import SeqIO
from Bio.Seq import Seq

def turn_seq(seq):
    """ Returns complementary sequence """
    ttable = str.maketrans("ACGTNacgtn","TGCANtgcan")
    seq = list(seq.translate(ttable))
    seq.reverse()
    seq = "".join(seq)
    return seq

def load_sequences(contigFile,delimiter):
    """ This function is needed to load a set of fasta sequences into memory. 
    Params
    --------
    
    contigFile --> contains the name of the fasta file you want to load into memory
    delimiter --> Leave empty in case you don't want to split the header, set the separator you want to use if you want to split the fasta header
    
    Return
    -------
    
    This function will return a dictionary which contained the headers as keys and the sequences as values.
    
    """
    seqs = {}
    name = ""
    s = []
    for line in open(contigFile):
        line = line.strip()
        if ">" in line:
            if name != "":
                seqs[name] = "".join(s)
            if delimiter == "":
                name = line.replace(">","")
            else:
                name = line.replace(">","").split(delimiter)[0]
            s = []
        else:
            s.append(line.upper())
    if len(s) != 0:
        seqs[name] = "".join(s)
    return seqs

parser = argparse.ArgumentParser(description="Script used to count number of telomeric repeats at the edges of contigs")
parser.add_argument("-g",dest="genomeFile",action="store",required=True,help="genomeFile")
parser.add_argument("-r",dest="repeat",action="store",required=True,help="Expected repeat")
parser.add_argument("-l",dest="length",action="store",type=int,default=1000,help="Length of the telomeric region to scan")
args = parser.parse_args()


repeat = args.repeat
repeat_rev = turn_seq(repeat)

l = args.length

seqs = load_sequences(args.genomeFile," ")

contigs = list(seqs.keys())
contigs.sort()
print("ContigName\tContig Length\tNR forward 3'\tNR reverse 3'\tNR forward 5'\tNR reverse 5'")
for code in contigs:
    s1 = seqs[code][:l]
    s2 = seqs[code][-l:]
    c1 = s1.count(repeat)
    c2 = s1.count(repeat_rev)
    c3 = s2.count(repeat)
    c4 = s2.count(repeat_rev)
    if c1 != 0 or c2 != 0 or c3 != 0 or c4 != 0:
        print(code+"\t"+str(len(seqs[code]))+"\t"+str(c1)+"\t"+str(c2)+"\t"+str(c3)+"\t"+str(c4))
