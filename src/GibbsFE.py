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
    The Gibbs Free Energy of a region within the potential promoter sequence is used as a metric to calculate a likelihood in pHunter. Upon comparing this likelihood with others, a log-likelihood may be formed, and scores of potential promoters can be determined.  
    This is the case since Gibbs Free Energy describes the ease at which a particular sequence unzips. 
    In order for Gibbs Free Energy of a region to be determined, sequences must be parsed for dinucleotides and trinucleotides, and the Gibbs Free Energy of these nucleotides must be outputted. 
    
    The GibbsFE_Dict function takes in a file name in the csv format and parses the file using the csv.reader() function to output a Gibbs Free Energy dictionary. 
    
    Input - file_name: str - Name of the .csv file with the dinucleotide/trinucleotide combinations and Gibbs Free Energy values
    Ouptut - GE_val: dictionary - Dictionary where the key is the dinucleotides/trinucleotides and the value is the Gibbs Free Energy value. 
'''
def GibbsFE_Dict(file_name):
    GE_val = {}
    
    #Sets directory for a .csv file in the data folder in the local directory
    file_name = '../data/' + file_name
    
    #Opens and parses the .csv file
    with open(file_name) as gcsv:
        reader = csv.reader(gcsv)
        for items in reader: #iterates through the csv file
            GE_val[items[0]] = float(items[1]) #key = dinucleotide. value = Gibbs Free Energy
    return(GE_val)


'''Computes the Gibbs Free Energy of a DNA sequence
    The Gibbs Free Energy of a region within the potential promoter sequence is used as a metric to calculate a likelihood in pHunter. Upon comparing this likelihood with others, a log-likelihood may be formed, and scores of potential promoters can be determined.  
    This is the case since Gibbs Free Energy describes the ease at which a particular sequence unzips. 
    In order for Gibbs Free Energy of a region to be determined, sequences must be parsed for dinucleotides and trinucleotides, and the Gibbs Free Energy of these nucleotides must be outputted. 
    
    The GibbsFE function takes in a sequence to be parsed, breaks down the sequence into its dinucleotide and trinucleotide composition, references its GibbsFEdict parameter to associate a particular dinucleotide or trinucleotide with a Gibbs Free Energy value, and constructs a list of these energy values in the order they appear in the sequence.
    The function also possesses the ability to calculate a moving average depending on the inputted moving_average_wlen value, where an average of the vector values that were originally returned is taken to smooth the returned array. The increment at which the moving average is taken will depend on the inputted parameter. 
    
    Inputs:
        mySeq: Seq Object - Sequence to be parsed
        GibbsFEdict: Dictionary - Dictionary where the keys are the particular dinucleotide/trinucleotides and the values are the Gibbs Free Energy values. 
        moving_average_wlen: int - Integer describing the number of Gibbs Free Energy values that will be averaged to create a Moving Average vector.
    
    Returns:
        myVec: List - List containing Gibbs Free Energy values of dinucleotide/trinucleotides of sequence
        myMAvec: List - List containing moving average of Gibbs Free Energy values of sequence. 
'''
def GibbsFE(mySeq, GibbsFEdict, moving_average_wlen=None):
    mySeq = mySeq.upper()   #upper case sequence
    
    #GibbsFEdict defined in csv file named in settings.json
    LUTfe = GibbsFEdict
    
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
        myVec.append(LUTfe[dinuc])   #look up GE value in dictionary using dinuc key
    
    #perform moving average, if requested (moving_average_wlen!=None)
    if (moving_average_wlen):
        myVeclen = len(myVec)
       
        #Sets moving_average_wlen to myVec length if the moving average is larger than the length of the myVec list. 
        if moving_average_wlen>myVeclen:
            moving_average_wlen = myVeclen
            
        myMAvec = []
        
        #Iterating through myVec list, grabbing multiple Gibbs Free Energy values to calculate an average. Appends average and continues iterating to create the moving average. 
        for p in range(myVeclen-moving_average_wlen+1):
            #compute average
            myMAvec.append(sum(myVec[p:p+moving_average_wlen])/moving_average_wlen)
        #return moving average vector
        return(myMAvec)
    else:
        return(myVec)  

'''Determines the minimum and maximum values of the Gibbs Free Energy dictionary
    The Gibbs Free Energy of a region within the potential promoter sequence is used as a metric to calculate a likelihood in pHunter. Upon comparing this likelihood with others, a log-likelihood may be formed, and scores of potential promoters can be determined.  
    This is the case since Gibbs Free Energy describes the ease at which a particular sequence unzips. 
    In order for the Gibbs Free Energy vector (created by GibbsFE) to become a likelihood, the minimum and maximum values must be determined for the bin edge array (np object). This bin_edge array, containing the ranges of each Gibbs Free Energy bin, allows for the creation of the array containing frequencies of Gibbs Free Energy values within bin_edge ranges(Gibbs Likelihood)
    
    The GibbsFE_minmax function takes in a Gibbs Free Energy dictionary (produced by GibbsFE_Dict), creates a list containing all of the dictionary's values, and iterates through the list to return the minimum and maximum values within the Gibbs Free Energy dictionary. 
    
    Input - file_name: str - Name of the .csv file with the dinucleotide values
    Output - max_GE,min_GE: int - Lowest values of the dictionary's values
'''
def GibbsFE_minmax(GibbsFEdict):
    LUTfe = GibbsFEdict #Initializes Gibbs Free Energy dictionary
    GE_val = [] #Initializes Gibbs Free Energy values list
    #Iterates through the Gibbs Free Energy dictionary values to append the values to a list
    for values in LUTfe.values(): #creates a list of GE values
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
            
    return(max_GE,min_GE)
        

