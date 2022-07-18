from Bio import SeqIO
from Bio.SeqUtils import GC as gc
from Bio import motifs
from GibbsFE import GibbsFE, GibbsFE_minmax, GibbsFE_Dict
import json
import numpy as np
import matplotlib.pyplot as plt
import csv
import math
import warnings

def GibbsFE_likelihood(Seq_Object, sampling_size, GibbsFEdict, GibbsFE_n_bins, GibbsFE_pc):
    '''
    Determines the Gibbs Free Energy likelihood-estimate for a given reference set of sequences
     Gibbs Free Energy contributions to RNA polymerase binding are computed in pHunter as log-likelihood ratios. 
     This requires that likelihoods be computed for “true” (reference) promoters and for some reference (null) background.
     
     The GibbsFE_likelihood function takes in a Seq Object representing a positive/reference set and uses the GibbsFE function to translate segments of the sequence into Gibbs Free Energy vectors.
     Each segment's Gibbs Free Energy vector is averaged and appended to a list of averages of segments (ener_set). ener_set is then used to construct an array of frequencies, thereby generating an estimate of the likelihood for that set.
     
     Inputs:
         Seq_Object: Seq Object - Seq Object of sequence or sequences to generate a Gibbs Free Energy vector from
         sampling_size: Int - The window size by which the Gibbs Free Energy vector is averaged across. Valued as the amount of base pairs around the -10 sequence that get unzipped initially and applied to both the positive and negative set for fair comparison. 
         GibbsFEdict: Dict - containing keys (dinucleotides/trinucleotides) and values (Gibbs Free Energy values)
         GibbsFE_n_bins: Int - Number of bins created for the hist_ener and hist_freq arrays
         GibbsFE_pc: Int - Pseudocount value added to both the sum of the hist_ener array and each value of hist_ener. 
     Returns:
         hist_freq: nparray - Array containing the frequencies of segments in each bin
         bin_edge: nparray - Array containing the ranges of each bin
    '''
    ener_set = []
    hist_freq = []
      
    #Iterates through seq object and uses GibbsFE function to determine Gibbs Free Energy of each sequence in the seq object using a sliding window (sample_size) within each sequence. Appends the average for each Gibbs Free Energy vector to ener_set
    for aseq in Seq_Object:
        
        #Raises a warning if a sequence is less than the sampling size 
        if len(aseq.seq)<sampling_size:   
            warnings.warn(aseq.id + " - sequence is less than the sampling size. Sequence will be discarded for the creation of the positive or negative set")
            continue
        
        #Iterates through sequence, taking segments (sampling_size base pairs long) and converting these segments into Gibbs Free Energy values. Averages the Gibbs Free Energy values within each segment and appends to ener_set list
        for p in range(0,len(aseq.seq)-sampling_size+1): 
            es = GibbsFE(aseq.seq[p:p+sampling_size],GibbsFEdict)   #Calls GibbsFE function and gets dinucleotide combos within the step size window
            ener_set.append(sum(es)/len(es))   #stores average energy for the sequence
              
    max_GE,min_GE = GibbsFE_minmax(GibbsFEdict)
    
    #hist_ener = nparray containing counts of segments from (ener_set) in each bin identified in bin_edge
    #bin_edge = nparray containing ranges of each bin used in hist_ener
    hist_ener,bin_edge = np.histogram(ener_set, range=(min_GE , max_GE), density=False, bins=GibbsFE_n_bins)
    
    #Creates a likelihood estimate by dividing each hist_ener array value by the sum of hist_ener array
    hist_freq = (hist_ener+GibbsFE_pc)/(hist_ener+GibbsFE_pc).sum()
    return hist_freq, bin_edge


def pHunt(mySeq, lmot, lthrs, rmot, rthrs, minD, maxD, unzip_dist, posdist, negdist, binedges,GibbsFEDict):
    '''
    Predicts putative promoters in a given sequence, following the PromoterHunter approach.
     The final score of a putative promoter is determined by adding the log-likelihoods of the left motif, right motif, and Gibbs Free Energy. 
     This requires each putative promoter to be evaluated for the -35 motif and -10 motif, with restrictions placed on the possible locations of the -10 sequence based on the possible spacer values. 
     In addition, Gibbs Free Energy must be applied to the moving region and converted into a log likelihood. 

     The pHunt function takes in a sequence to iterate through and checks each potential location with the left motif for a value above the inputted threshold. 
     Sequences found to satisfy this requirement are tested with the right motif, where the right -10 sequence could be any distance from the minimum spacer value to the maximum spacer value away from the -35 sequence. 
     If this criteria is also satisfied, the sequence is reported as a promoter and motif scores are summed. The Gibbs Free Energy of a specified range (lerg and rerg) are used to create conditional probabilities for the promoter set and the background set to calculate a log likelihood ratio. The sum of the motif scores and Gibbs Free Energy log-likehood ratio are summed together for a final score. 
     
    Inputs:
     mySeq: Seq - DNA sequence to be iterated across
     lmot: Motif Object - -35 motif 
     lthrs: Int - Score threshold for -35 motif
     rmot: Motif Object - -10 motif 
     rthrs: Int - Score threshold for -10 motif
     minD: Int - Minimum spacer value
     maxD: Int - Maximum spacer value
     unzip_dist: Int - Distance around the start of the -10 sequence where unzipping begins and the Gibbs Free Energy Likelihood is evaluated across 
     posdist: np Array - Contains Gibbs Free Energy likelihood of promoter set
     negdist: np Array - Contains Gibbs Free Energy likelihood of negative set
     binedges: np Array - Contains bin ranges of positive and negative arrays 
     GibbsFEDict: Dictionary - Dictionary containing keys (dinucleotides/trinucleotides) and values (Gibbs Free Energy values)
    
    Returns:
     hmsp: List - List of promoter sequences with -10 and -35 sequences that surpass the motif thresholds. List composed of dictionaries with information on each identified promoter. 
    '''
    #GC = gc(mySeq)/100
    
    #score sequence with left motif
    lscrs=lmot.pssm.calculate(mySeq).tolist()
    
    #score sequence with left motif
    rscrs=rmot.pssm.calculate(mySeq).tolist()
    
    
    #compute the Gibbs Free Energy vector on sequence
    fes = GibbsFE(mySeq, GibbsFEDict)
    
    #list of spacer-adequate high-scoring motif pairs
    hsmp = []
    
    mySeqlen = len(mySeq)
    #go through sequence, first with left PSSM. Iterate through sequence with consideration for length of minimum spacer and both motifs. 
    for pl in range(mySeqlen-(lmot.pssm.length + minD + rmot.pssm.length)):
        
        if lscrs[pl] > lthrs:
            remrange = pl+lmot.pssm.length+maxD+1 if (pl+lmot.pssm.length+maxD+1 < mySeqlen - rmot.pssm.length) else mySeqlen - rmot.pssm.length #else statement for potential promoters close to the end of the sequence
            
            #go through sequence, now with right PSSM. Considers placement of right motif depending on any possible permutation off spacer length
            for pr in range(pl+lmot.pssm.length+minD,remrange):
                
                #if score above threshold
                if rscrs[pr] > rthrs:
                    
                    #Initializes information dictionary
                    element = {'lpos' : pl, 'lseq' : str(mySeq[pl:pl+lmot.pssm.length]),
                               'lscr' : lscrs[pl],
                               'rpos' : pr, 'rseq' : str(mySeq[pr:pr+rmot.pssm.length]),
                               'rscr' : rscrs[pr],
                               'spcr' : str(mySeq[pl+lmot.pssm.length:pr].lower()),
                              }
                    #compute overall (normalized) PSSM contribution
                    element['mscr'] = element['lscr'] + element['rscr']

                    #add Gibbs Free Energy component
                    #Initializes lrange and rrange as unzip_dist from the start of the -10 sequence
                    #Prevents lrange or rrange from being out of bounds (extending beyond the sequence)
                    lrange = pr-(unzip_dist//2) if pr-(unzip_dist//2)>0 else 0
                    rrange = pr+(unzip_dist//2) if pr+(unzip_dist//2)<mySeqlen else mySeqlen
                    
                    #Averages values from GibbsFE vector
                    element['escr'] = sum(fes[lrange:rrange])/(rrange-lrange)
                    
                    #Gets the frequency of the first bin with a higher Gibbs FE value than 'escr' for both the positive and negative set
                    p_pos = posdist[np.where(binedges>element['escr'])[0][0] - 1]
                    p_neg = negdist[np.where(binedges>element['escr'])[0][0] - 1]
                    
                    #Creates the log likelihood for Gibbs Free Energy
                    element['ellr'] = np.log2(p_pos / p_neg)
                    
                    #compute final score (global norm motif score + 2*norm energy score)
                    element['Fscr'] = element['mscr'] + element['ellr']
                    
                    #add high-scoring motif pair with all info to the list of HSMPs
                    hsmp.append(element)
                  
    #return list of spacer-adequate high-scoring motif pairs
    return hsmp

def go():
    '''
    Transforms the code into a program by:
     Getting information from settings.json,
     Initalizing parameters,
     Running the pHunt function, 
     Creating output file 
    '''
    #Stores settings.json in a local variable (json_file)
    with open('settings.json') as settings:
        json_file = json.load(settings) 
        
    #Parses both fasta files. File names stored within the json file and searches for the files within the data folder in the local directory.
    l_insta = SeqIO.parse("../data/" + json_file["left_motif"],'fasta')
    r_insta = SeqIO.parse("../data/" + json_file["right_motif"],'fasta')
    
    #Creates the -35 motif and assigns a psuedocount and threshold based on values in json_file (settings.json)
    lmot = motifs.create([l.seq for l in l_insta],'ACGT') 
    lmot._pseudocounts = json_file["left_motif_pseudocounts"]
    lthrs = json_file["left_motif_threshprec"]
    
    #Creates the -10 motif and assigns a psuedocount and threshold based on values in json_file (settings.json)
    rmot = motifs.create([r.seq for r in r_insta],'ACGT') #Creates the -10 motif
    rmot._psuedocounts = json_file["right_motif_pseudocounts"]
    rthrs = json_file["right_motif_threshprec"]

    #Accesses minimum and maximum spacer values from json_file (settings.json) 
    minD = json_file["min_spacer_value"]
    maxD = json_file["max_spacer_value"]
    
    #Distance around start of -10 sequence where unzipping begins 
    unzip_dist = json_file["sampling_size"]
    
    #Initializes the GibbsFE Dictionary from csv file in settings.json
    GibbsFEdict = GibbsFE_Dict(json_file["GibbsFE_file"])
    
    #Parses both the positive and background set. Creates a Seq Object to be passed into the GibbsFE_likelihood function
    seq_objectp = SeqIO.parse("../data/" + json_file["promoter_set"],'fasta')
    seq_objectb = SeqIO.parse("../data/" + json_file["background_set"],'fasta')
    
    #If the negative set is a whole genome (len(list{seq_objectb) == 1), the genome will be divided into multiple segments. Each segment will be then be ran into the GibbsFE_Likelihood function
    #If the negative set is a collection of intergenic sequences (len(list{seq_objectb) > 1), the set will be left as is
    if len(list(seq_objectb)) == 1:  
        
        #Copies the positive and background set Seq Objects to determine scan_wind and to iterate upon
        seq_objects = SeqIO.parse("../data/" + json_file["promoter_set"],'fasta')
        seq_objecti = SeqIO.read("../data/" + json_file["background_set"],'fasta')
        
        #Prepares seq_objectb to be rewritten as a list of segments
        seq_objectb = []
        
        #Scan_wind is the length of the segments that the genome will be divided into. It is calculated as the average length of the promoter sequences rounded to the nearest integer.
        w = [len(seq) for seq in seq_objects]
        scan_wind= round(sum(w)/len(w))
        
        #Divides the genome set into segments and appends these segments to the seq_objectb list
        for q in range(0,len(seq_objecti.seq)-scan_wind + 1,json_file["step_size_genome"]):
            seq_objectb.append(seq_objecti[q:q+scan_wind])
               
    
    #Creates the positive and negative sets by running the GibbsFE_likelihood() function on the promoter and non-promoter Seq Objects (parsed from the fasta files), using the GibbsFEdict, bin amount, psuedocount, and step size from the json_file.
    posdist,bin_edges = GibbsFE_likelihood(seq_objectp, unzip_dist, GibbsFEdict, json_file["GibbsFE_n_bins"], json_file["GibbsFE_pseudocounts"]) 
    negdist,bin_edges = GibbsFE_likelihood(seq_objectb, unzip_dist, GibbsFEdict, json_file["GibbsFE_n_bins"], json_file["GibbsFE_pseudocounts"]) 
    
    #plt.bar(range(json_file["GibbsFE_n_bins"]), negdist, alpha=0.5) #for graphing the GE frequencies of both sets
    #plt.bar(range(json_file["GibbsFE_n_bins"]), posdist, alpha=0.5)
        
    #Parses the input_sequences file stored in the json file, runs the pHunt function on each resulting sequence from parsing (using the variables initialized above), and returns a csv file with promoter information   
    mySeq = SeqIO.parse("../data/" + json_file["input_sequences"],'fasta')
    
    
    with open(json_file["output_information"], 'w',newline = '') as f:
        writer = csv.writer(f)
        
        #Writes headers
            #ID: String - ID of Input Sequence for which the promoter is generated for 
            #Spacer_length: Integer - Length of the spacer between the -35 and -10 
            #Range: String - Beginning of -35 sequence to end of -10 sequence
            #Left_Motif_Score: Float - Score of -35 sequence using -35 motif pssm
            #Right_Motif_Score: Float - Score of -10 sequence using -10 motif pssm
            #Average_GibbsFE: Float - Gibbs Free Energy of region of ##https://www.google.com/search?q=where+does+unzipping+start+for+transcription&client=opera&hs=mrD&sxsrf=ALiCzsaARycheqIz7p4mDX2kcxZ0T4HL8w%3A1657053195252&ei=C6DEYt_iDt2hptQP5Mi8uA8&ved=0ahUKEwiflM6YzOL4AhXdkIkEHWQkD_cQ4dUDCA0&uact=5&oq=where+does+unzipping+start+for+transcription&gs_lcp=Cgdnd3Mtd2l6EAMyBQghEKABMgUIIRCgAToHCAAQRxCwAzoFCCEQqwI6CAghEB4QFhAdSgQIQRgASgQIRhgAULwHWM0XYOkYaAFwAXgAgAFiiAG_CpIBAjE4mAEAoAEByAEIwAEB&sclient=gws-wiz
            #GibbsFE_Log-Likelihood: Numpy Float - Log-Likelihood of promoter and background sets given the Average_GibbsFE of the unzipping region
            #Final_Score: Numpy float - Sum of Left_Motif_Score, Right_Motif_Score, and GibbsFE_Log-Likelihood
            #Promoter_Sequence: String - Nucleotides of -35 sequence, spacer, and -10 sequence. #https://www.addgene.org/mol-bio-reference/promoters/
        writer.writerow(['ID','Spacer_length','Range','Left_Motif_Score','Right_Motif_Score','Average_GibbsFE','GibbsFE_Log-Likelihood','Final_Score', 'Promoter_Sequence'])
        
        #Iteration to write information on each promoter in csv file
        for items in mySeq:
            
            #Storing high-scoring motif pairs in mypromoters (dictionary information) for each sequence in the mySeq Seq Object
            mypromoters = pHunt(items.seq,lmot,lthrs,rmot,rthrs,minD,maxD,unzip_dist, posdist, negdist, bin_edges, GibbsFEdict)
            
            #sorting results by global score
            mypromoters = sorted(mypromoters, key = lambda k: k['Fscr'], reverse=True)
            
            #Iterating through mypromoters to write information about each reported promoter on output file
            for promoter in mypromoters:
                writer.writerow([
                    items.id,
                    promoter['rpos'] - promoter['lpos'] - 6,
                    str(promoter['lpos']) + '-' + str(promoter['rpos'] + len(promoter['rseq'])),
                    promoter['lscr'],
                    promoter['rscr'], 
                    promoter['escr'],
                    promoter['ellr'],
                    promoter['Fscr'],
                    promoter['lseq']+promoter['spcr']+promoter['rseq']
                    ])
    f.close()

#Creates output csv file upon running the Python file
if __name__ == "__main__":
    go()


