from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from Bio.SeqUtils import GC as gc
import numpy as np 
import matplotlib.pyplot as plt
import json
import csv

def GibbsFE_Dict(file_name):
    '''
    Reads a csv file and outputs a GE value dictionary
     The Gibbs Free Energy of a region within the potential promoter sequence is used as a metric to calculate a likelihood in pHunter. Upon comparing this likelihood with others, a log-likelihood may be formed, and scores of potential promoters can be determined.  
     This is the case since Gibbs Free Energy describes the ease at which a particular sequence unzips. 
     In order for Gibbs Free Energy of a region to be determined, sequences must be parsed for dinucleotides and trinucleotides, and the Gibbs Free Energy of these nucleotides must be outputted. 
     
     The GibbsFE_Dict function takes in a file name in the csv format and parses the file using the csv.reader() function to output a Gibbs Free Energy dictionary. 
     
     Input:
         file_name: Str - Name of the .csv file with the dinucleotide/trinucleotide combinations and Gibbs Free Energy values
     Ouptut:
         GE_val: Dict - Dictionary where the key is the dinucleotides/trinucleotides and the value is the Gibbs Free Energy value. 
    '''
    GE_val = {}
    
    #Sets directory for a .csv file in the data folder in the local directory
    file_name = '../data/' + file_name
    
    #Opens and parses the .csv file
    with open(file_name) as gcsv:
        reader = csv.reader(gcsv)
        for items in reader: #iterates through the csv file
            GE_val[items[0]] = float(items[1]) #key = dinucleotide. value = Gibbs Free Energy
    return GE_val


def GibbsFE(mySeq, GibbsFEdict):
    '''
    Computes the Gibbs Free Energy of a DNA sequence
     The Gibbs Free Energy of a region within the potential promoter sequence is used as a metric to calculate a likelihood in pHunter. Upon comparing this likelihood with others, a log-likelihood may be formed, and scores of potential promoters can be determined.  
     This is the case since Gibbs Free Energy describes the ease at which a particular sequence unzips. 
     In order for Gibbs Free Energy of a region to be determined, sequences must be parsed for dinucleotides and trinucleotides, and the Gibbs Free Energy of these nucleotides must be outputted. 
     
     The GibbsFE function takes in a sequence to be parsed, breaks down the sequence into its dinucleotide and trinucleotide composition, references its GibbsFEdict dictionary to associate a particular dinucleotide or trinucleotide with a Gibbs Free Energy value, and constructs a list of these energy values in the order they appear in the sequence.
     
     Inputs:
         mySeq: Seq - Sequence to be iterated on
         GibbsFEdict: Dict - Dictionary where the keys are the particular dinucleotide/trinucleotides and the values are the Gibbs Free Energy values. 
     
     Returns:
         myVec: List - List containing Gibbs Free Energy values of dinucleotide/trinucleotides of sequence

    '''
    mySeq = mySeq.upper()   #upper case sequence
    
    #convert sequence to vector
    myVec = []
    
    #Sets variable to length of GibbsFEdict key (dinucleotide or trinucleotide set)
    nuclen = len(list(GibbsFEdict.keys())[0])
    
    #return empty vector if sequence is less than 3 nucleotides
    mySeqlen = len(mySeq)
    if mySeqlen<3:
        return(myVec)
    
    #iterate over the sequence, grabbing nuclen amount of elements at a time
    #with a +1 step size
    for p in range(mySeqlen-1):
        
        multinuc = mySeq[p:p+nuclen] #create multinucleotide (length of multinucleotide of sequence = nuclen [Length of dinucleotide/trinucleotide key in GibbsFEdict])
        myVec.append(GibbsFEdict[multinuc])   #look up GE value in dictionary using dinucleotide or trinucleotide key
        
    return myVec 


def GibbsFE_minmax(GibbsFEdict):
    '''
    Determines the minimum and maximum values of the Gibbs Free Energy dictionary
     The Gibbs Free Energy of a region within the potential promoter sequence is used as a metric to calculate a likelihood in pHunter. Upon comparing this likelihood with others, a log-likelihood may be formed, and scores of potential promoters can be determined.  
     This is the case since Gibbs Free Energy describes the ease at which a particular sequence unzips. 
     In order for the Gibbs Free Energy vector (created by GibbsFE) to become a likelihood, the minimum and maximum values must be determined for the bin edge array (np object). This bin_edge array, containing the ranges of each Gibbs Free Energy bin, allows for the creation of the array containing frequencies of Gibbs Free Energy values within bin_edge ranges(Gibbs Likelihood)
     
     The GibbsFE_minmax function takes in a Gibbs Free Energy dictionary (produced by GibbsFE_Dict), creates a list containing all of the dictionary's values, and iterates through the list to return the minimum and maximum values within the Gibbs Free Energy dictionary. 
     
     Input:
         file_name: Str - Name of the .csv file with the dinucleotide values
     Output:
         max_GE,min_GE: Int - Lowest values of the dictionary's values
    '''
    GibbsFEdict #Initializes Gibbs Free Energy dictionary
    GE_val = [] #Initializes Gibbs Free Energy values list
    
    #Iterates through the Gibbs Free Energy dictionary values to append the values to a list
    for values in GibbsFEdict.values(): #creates a list of GE values
        GE_val.append(values)
        
    #Initializes min and max variables as a random reference
    min_GE = GE_val[0] 
    max_GE = GE_val[0]
    
    #Iterates through Gibbs Free Energy values list searching for minimum and maximum values
    for items in GE_val: 
        if min_GE>items:
            min_GE = items
        if items>max_GE:
            max_GE = items
            
    return max_GE,min_GE
        

