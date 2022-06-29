from Bio import SeqIO
from Bio.SeqUtils import GC as gc
from Bio import motifs
from GibbsFE import GibbsFE
from GibbsFE import Gibbs_minmax
import json
import numpy as np
import matplotlib.pyplot as plt

#parameters ener_sets = (wsize, promoterfield)

'''Determines the GE values of a given sequence
    Inputs:
        file_name: str - File name of the fasta sequence to be parsed
        step_size: intr - Sliding window used to evaluate the inputted sequence
    Returns:
        hist_freq: nparray - Array containing the frequencies of segments in each bin
        bin_edge: nparray - Array containing the ranges of each bin
'''
def ener_sets(file_name,step_size):
    ener_set = []
    hist_freq = []
    n_bins = 30
    pc = 1
    #Proceeds if there are multiple instances in the fasta file
    seq_object = SeqIO.parse(file_name,'fasta') #initializes the seq object
    for aseq in seq_object:
        for p in range(0,len(aseq)-step_size+1,step_size): #sets window
            es = GibbsFE(aseq[p:p+step_size], 0)   #Calls GibbsFE function and gets dinucleotide combos within the window
            ener_set.append(sum(es)/len(es))   #stores  average energy for the sequence
    max_GE,min_GE = Gibbs_minmax('GibbsFE.csv')
    hist_ener,bin_edge = np.histogram(ener_set, range=(min_GE , max_GE), density=False, bins=n_bins)
    #hist_ener = nparray containing counts of segments from (ener_set) in each bin
    hist_freq = (hist_ener+pc)/(hist_ener+pc).sum()
    return hist_freq, bin_edge



"""Predicts putative promoters in a given sequence, following the PromoterHunter approach.
   Inputs:
   - mySeq - DNA sequence [Seq, not SeqRecord object]
   - lmot - -35 motif object
   - lthrs - score threshold for -35 motif
   - rmot - -10 motif object
   - rthrs - score threshold for -10 motif
   - [minD, maxD] - range of spacer lengths
   - wsize - size of moving average window for Gibbs free energy
   - [lerg,rerg] - range surrounding -10 on which to compute energy score
   
   Returns:
   - List of High-scoring Motif Pairs (HSMP), with all the score and sequence information
     for the constituent motifs and energy
"""    
def pHunt(mySeq, lmot, lthrs, rmot, rthrs, minD, maxD, wsize, posdist, negdist, binedges):
    lerg = -wsize//2   #left window from -10 for energy
    rerg = +wsize//2    #right window from -10 for energy
    GC = gc(mySeq)/100
    #score sequence with left motif
    lscrs=lmot.pssm.calculate(mySeq).tolist()
    #score sequence with left motif
    rscrs=rmot.pssm.calculate(mySeq).tolist()
    #compute the Gibbs Free Energy vector
    fes = GibbsFE(mySeq, None)
    
    #list of spacer-adequate high-scoring motif pairs
    hsmp = []
    
    mySeqlen = len(mySeq)
    #go through sequence, first with left PSSM    
    for pl in range(mySeqlen-(lmot.pssm.length + minD + rmot.pssm.length)):
        #if score above threshold (Makes sure not to consider promoter seqs with no space for right sequence and adjusts the range of rightmotif start values if it no longer can extend with maxD value)
        if lscrs[pl] > lthrs:
            remrange = pl+lmot.pssm.length+maxD+1 if (pl+lmot.pssm.length+maxD+1 < mySeqlen - rmot.pssm.length) else mySeqlen - rmot.pssm.length
            #go through sequence, now with right PSSM, up to spacer
            for pr in range(pl+lmot.pssm.length+minD,remrange):
                #if score above threshold
                if rscrs[pr] > rthrs:
                    element = {'lpos' : pl, 'lseq' : str(mySeq[pl:pl+lmot.pssm.length]),
                               'lscr' : lscrs[pl],
                               'rpos' : pr, 'rseq' : str(mySeq[pr:pr+rmot.pssm.length]),
                               'rscr' : rscrs[pr],
                               'spcr' : str(mySeq[pl+lmot.pssm.length:pr].lower()),
                              }
                    #compute overall (normalized) PSSM contribution
                    element['mscr'] = element['lscr'] + element['rscr']

                    #add Gibbs Free Energy component
                    #as average betweeen coordinates (plus lerg, rerg margins)
                    lrange = pr+lerg if pr+lerg>0 else 0
                    rrange = pr+rerg if pr+rerg<mySeqlen-rmot.pssm.length else mySeqlen-rmot.pssm.length #check this
                    element['escr'] = sum(fes[lrange:rrange])/(rrange-lrange)
                    
                    p_pos = posdist[np.where(binedges>element['escr'])[0][0] - 1] #comeback
                    p_neg = negdist[np.where(binedges>element['escr'])[0][0] - 1]
                    
                    element['ellr'] = np.log2(p_pos / p_neg)
                    
                    #compute final score (global norm motif score + 2*norm energy score)
                    element['Fscr'] = element['mscr'] + element['ellr']
                    
                    #add high-scoring motif pair with all info to the list of HSMPs
                    hsmp.append(element)
                  
    #return list of spacer-adequate high-scoring motif pairs
    return(hsmp)
'''Initalizes parameters, runs PHunt function, and creates output file 
'''
def go():
    n_bins = 30
    #local variable json_variable
    with open('settings.json') as settings:
        json_file = json.load(settings) #Stores settings.json in a local variable
    l_insta = SeqIO.parse("../data/" + json_file["left_motif"],'fasta')
    r_insta = SeqIO.parse("../data/" + json_file["right_motif"],'fasta')
    mySeq = SeqIO.read("../data/" + json_file["input_sequences"],'fasta')
    lmot = motifs.create([l.seq for l in l_insta],'ACGT') #Creates the -35 motif
    lmot._pseudocounts = 0.05
    lthrs = lmot.pssm.distribution(precision=10**3).threshold_patser()
    rmot = motifs.create([r.seq for r in r_insta],'ACGT') #Creates the -10 motif
    rmot._psuedocounts = 0.05
    rthrs = rmot.pssm.distribution(precision=10**3).threshold_patser()
    minD = json_file["min_spacer_value"]
    maxD = json_file["max_spacer_value"]
    posdist,bin_edges = ener_sets("../data/" + json_file["promoter_set"],json_file["step_size_pos"]) #positive set GE frequencies
    negdist,bin_edges = ener_sets("../data/" + json_file["background_set"],json_file["step_size_neg"]) #negative set GE frequencies stored
    #plt.bar(range(n_bins), negdist, alpha=0.5) #for graphing the GE frequencies of both sets
    #plt.bar(range(n_bins), posdist, alpha=0.5)
    mypromoters = pHunt(mySeq.seq,lmot,lthrs,rmot,rthrs,minD,maxD,101, posdist, negdist, bin_edges)
    #sorting results by global score
    mypromoters = sorted(mypromoters, key = lambda k: k['Fscr'], reverse=True)
    '''for promoter in mypromoters: #reporting all promoters
        print('Range:', promoter['lpos'],'-',promoter['rpos'], '|', promoter['rpos']-promoter['lpos']-6, \
              '- L:',round(promoter['lscr'],2),'R:',round(promoter['rscr'],2),\
              'FE:',round(promoter['escr'],2), 'FELLR:',round(promoter['ellr'],2), 'T:',round(promoter['Fscr'],2),\
             '- S:',promoter['lseq']+promoter['spcr']+promoter['rseq'])'''
    f = open("Output_Information.txt", "w") #Output Text file
    for promoter in mypromoters:
        f.write(('Range: ' + str(promoter['lpos']) +'-' + str(promoter['rpos']) + ' | ' + 'Spacer_length:' + str(promoter['rpos']-promoter['lpos']-6)+ ' Left_Motif_Score: '))
        f.write(str(round(promoter['lscr'],2)) + ' Right_Motif_Score: ' + str(round(promoter['rscr'],2)) + ' Gibbs_Energy: ' + str(round(promoter['escr'],2)) + ' Log-Likelihood_GibbsFE: ')
        f.write(str(round(promoter['ellr'],2)) + str( ' Score:' ) + str(round(promoter['Fscr'],2)) + ' Promoter_Sequence: ' + str(promoter['lseq']+promoter['spcr']+promoter['rseq']))
        f.write('\n')
    f.close()

#Creates Output_Information.txt upon running the Python file
if __name__ == "__main__":
    go()


