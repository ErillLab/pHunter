from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from Bio.SeqUtils import GC as gc
import numpy as np 
import matplotlib.pyplot as plt
import json
import csv

'''Reads a csv file and outputs a GE value dictionary
    Input - file_name: str - Name of the .csv file with the dinucleotide values
    Ouptut - GE_val: dictionary - key | dinucleotide value | Gibbs Free Energy
'''
def GibbsList(file_name):
    GE_val = {}
    with open(file_name) as gcsv:
        reader = csv.reader(gcsv)
        for items in reader: #iterates through the csv file
            GE_val[items[0]] = float(items[1]) #key = dinucleotide. value = Gibbs Free Energy
    return(GE_val)
# Define GibbsFE function
"""
GibbsFE: computes the Gibbs Free Energy of a DNA sequence
reads in:
- DNA sequence (BioPython Seq object)
- moving average window length (to compute Free Energy over)
Assumes the DNA sequence only contains DNA characters (A, C, G, T)
"""
def GibbsFE(mySeq, moving_average_wlen=None):
    mySeq = mySeq.upper()   #upper case sequence
    #Gibbs Free Sequence components defined in settings.json (manipulatable)
    LUTfe = GibbsList('../data/GibbsFE.csv')
    #convert sequence to vector
    myVec = []
    #return empty vector if sequence is less than 2 nucleotides
    mySeqlen = len(mySeq)
    if mySeqlen<2:
        return(myVec)
    #iterate over the sequence, grabbing two elements at time
    #with a +1 step size
    for p in range(mySeqlen-1):
        dinuc = mySeq[p]+mySeq[p+1]  #create dinucleotide
        myVec.append(LUTfe[dinuc])   #look up GE value in table
    
    #perform moving average, if requested (wlen!=None)
    if (moving_average_wlen):
        myVeclen = len(myVec)
        #return original vector if window size over length
        if myVeclen<moving_average_wlen:
            moving_average_wlen = myVeclen
        myMAvec = []
        for p in range(myVeclen-moving_average_wlen+1):
            #compute average
            myMAvec.append(sum(myVec[p:p+moving_average_wlen])/moving_average_wlen)
        #return moving average vector
        return(myMAvec)
    else:
        return(myVec)  

'''Determines the minimum and maximum values of the Gibbs Free Energy dictionary
    Input - file_name: str - Name of the .csv file with the dinucleotide values
    Output - max_GE,min_GE: int - Lowest values of the dictionary's values
'''
def Gibbs_minmax(file_name):
    LUTfe = GibbsList('../data/' + file_name)
    GE_val = []
    for values in LUTfe.values(): #creates a list of GE values
        GE_val.append(values)
    min_GE = GE_val[0] #Initializes variables
    max_GE = GE_val[0]
    for items in GE_val: #Iterates through list searching for minimum and maximum values
        if min_GE>items:
            min_GE = items
        if items>max_GE:
            max_GE = items
    return(max_GE,min_GE)
        
''' 
#create a sequence record/Test Code
mydna = SeqRecord(Seq('GTCAGGTGATTGACGAAGATGTCTATCCGATTCTGTCGCTGCAATCGTGCCTCGACAAGCGTGCGGCAAAAGGCGGCGTCTCACCGCAGCAGGTGGCGCAGGCGATTGCTTTTGCGCAGGCTCGGTTAGGGTAAGAACATTTATATGTATAAATTTGAGCCTGGCTTATCGCCGGGCTTTTTTATGGCAAAAAAAAGCGGATCCTGGAGATCCGCAAAAGTTCACGTTGGCTTTAGTTATTCGAGTTGAGAAACTCTCGAAACGGGCAGTGACTTCAAGGGTTAAAAGAGGTGCCGCTCCGTTTCTGTGAGCAATTATCAGTCAGAATGCTTGATAGGGAGCGCCGTTCATTGCTATTCTACCTATCGCCATGAACTATCGTGGCGATGGAGGATGGATAATGAATATTCGTGATCTTGAGTACCTGGTGGCATTGGCTGAACACCGCCATTTTCGGCGTGCGGCAGATTCCTGCCACGTTAGCCAGCCGACGCTTAGCGG'),\
                        id='NC_000913.3:4158090-4158590',\
                        description='PoxyR_Ecoli')    

#call GibbsFE function - Come Back
energies = GibbsFE(mydna,50)

#plot results
plt.plot(energies)
'''
