from Bio import SeqIO
from Bio.SeqUtils import GC as gc
from Bio import motifs
from GibbsFE import GibbsFE, GibbsFE_minmax, GibbsFE_Dict
import json
import numpy as np
import csv
import warnings
import statistics as st
from time import ctime

def Det_Nuc_Freq(sequence):
    '''
    Calculates GC content of a inputted sequence and determines frequency
    of each base pair based on GC content
    
    The original Promoter Hunter uses GC content to determine the null hypothesis
    when computed the PSSM LLR score. 
    By converting GC content into a frequency of each base pair, the original 
    Promoter Hunter uses the frequency of a base pair as the null instead 
    of a 0.25 frequency. 
    
    This function calculates the GC content of a sequence and determines the AT content
    as 1 - GC content. It then refers to the GC frequency to calculate the frequency
    of Guanine and Cytosine (GC frequency/2) and AT frequency to calculate the frequency
    of Adenine and Thymine (AT frequency/2)
    
    Note: Either the GC content null hypothesis or the 0.25 frequency null hypothesis
    can be used to calculate scores in this program
    
    Parameters
    ----------
    sequence : str - sequence to calculate the GC content of 

    Returns
    -------
    A_cont, T_cont, G_cont, C_cont: float - Frequencies of each base pair, used for background of PSSM scores. 

    '''
    #Calculates GC frequency
    GC_content_freq = gc(sequence) * (1/100)
        
    #Assigns frequency of nucleotides identically to original Promoter
    A_cont = (1-GC_content_freq)/2
    T_cont = (1-GC_content_freq)/2
    G_cont = (GC_content_freq)/2
    C_cont = (GC_content_freq)/2
        
    return {"A":A_cont,"C":C_cont,"G":G_cont,"T":T_cont}
    

def GibbsFE_likelihood(Seq_Object, sampling_size, GibbsFEdict, GibbsFE_n_bins, GibbsFE_pc):
    '''
    Determines the Gibbs Free Energy likelihood-estimate for a given reference set of sequences
     Gibbs Free Energy contributions to RNA polymerase binding are computed in pHunt as log-likelihood ratios. 
     This requires that likelihoods be computed for “true” (reference) promoters and for some reference (null) background.
     
     The GibbsFE_likelihood function takes in a Seq Object representing a positive/reference set and uses the GibbsFE function to translate segments of the sequence into averaged Gibbs Free Energy vectors.
     This moving average vector, avg_ener_sets, is passed into the np.histogram() package to make an array of counts. Converting this array into a frequency array, a likelihood estimate is created for the inputted set.
     
     Inputs:
         Seq_Object: List - List of sequence or sequences to generate a Gibbs Free Energy vector from
         sampling_size: Int - The window size by which the Gibbs Free Energy vector is averaged across. Valued as the amount of base pairs around the -10 sequence that get unzipped initially and applied to both the positive and negative set for fair comparison. Passed into GibbsFE as the moving_avg_wlen 
         GibbsFEdict: Dict - containing keys (dinucleotides/trinucleotides) and values (Gibbs Free Energy values)
         GibbsFE_n_bins: Int - Number of bins created for the hist_ener and hist_freq arrays
         GibbsFE_pc: Int - Pseudocount value added to both the sum of the hist_ener array and each value of hist_ener. 
     Returns:
         hist_freq: nparray - Array containing the frequencies of segments in each bin
         bin_edge: nparray - Array containing the ranges of each bin
    '''
    hist_freq = []
    avg_ener_sets = []
      
    #Iterates through list and uses GibbsFE function to determine Gibbs Free Energy of each sequence in the seq object using a sliding window (sample_size) within each sequence. Appends the average for each Gibbs Free Energy vector to ener_set
    for aseq in Seq_Object:
        
        #Raises a warning if any sequence is less than the sampling size 
        if len(aseq.seq)<sampling_size:   
            warnings.warn(aseq.id + " - sequence is less than the sampling size. Sequence will be discarded for the creation of the positive or negative set") 
            continue
        
        #Passes each sequence into GibbsFE function, along with a window size so that the function may iterate through each sequence, convert the sequences into GibbsFE values depending on their dinucleotide/trinucleotide composition, and take segments of this ener_array (sampling_size base pairs long) to calculate an average. With a step size of one, a moving average vector is compiled (list of the average energy values of the sequence). 
        #Extends an ongoing list with averages of each sequence computed from the GibbsFE function
        
        avg_ener_sets.extend(GibbsFE(aseq.seq,GibbsFEdict,sampling_size))
    
    #Determines the maximum and minimum Gibbs Free Energy values present in the inputted dictionary. These values will be used in the np.histogram() function. 
    max_GE,min_GE = GibbsFE_minmax(GibbsFEdict)
    
    #hist_ener = nparray containing counts of segments from (ener_set) in each bin identified in bin_edge
    #bin_edge = nparray containing ranges of each bin used in hist_ener
    hist_ener,bin_edge = np.histogram(avg_ener_sets, range=(min_GE , max_GE), density=False, bins=GibbsFE_n_bins)
    
    #Creates a likelihood estimate by dividing each hist_ener array value by the sum of hist_ener array
    hist_freq = (hist_ener+GibbsFE_pc)/(hist_ener+GibbsFE_pc).sum()

    return hist_freq, bin_edge


def pHunt(mySeq, lmot, lthrs, rmot, rthrs, GC_null_hypo, GC_pseudocnt, minD, maxD,use_GibbsFE, mode, lerg = None, rerg = None, wsize = None, GibbsFEDict = None, posdist = None, negdist = None, binedges = None):
    '''
    Predicts putative promoter regions (-35, spacer, -10) in a given sequence, following the PromoterHunter approach.
     The final score of a putative promoter is determined by summing the contributions of the left motif, right motif, and Gibbs Free Energy.
     In the original PromoterHunter, these contributions were normalized and added together with arbtirary coefficients. 
     In the revised PromoterHunter, these contributions are summed together as log-likelihood ratios.
     
     the pHunt function takes in a sequence to iterate through and scores each potential location for the -35 sequence, using the left motif. 
     Locations with scores above the left threshold are used to determine potential locations of the -10 sequence, given the minimum and maximum spacer values inputted by the user.
     These potential locations for the -10 sequence are scored by the right motif.
     If a -10 sequence possesses a score above the right threshold, the -35, spacer region, and -10 are reported as a hit, and Gibbs Free Energy is calculated.
     The Gibbs Free Energy of a specified range (lerg and rerg) are used to create conditional probabilities for the promoter set and the background set to calculate a log likelihood ratio. The sum of the motif scores and Gibbs Free Energy log-likehood ratio are summed together for a final score. 
     
     Note: As configurable by the user, this function can either compute final score as a sum of normalized values or of LLRs.
     
     Further description of how GibbsFE is turned into a log-likelihood ratio:
     PSSM scores computed on potential locations for the -35 and -10 motifs are already log-likelihood ratios.
     Gibbs Free Energy is converted into an LLR value by passing promoter and non-promoter distributions into the function.
     In order to convert Gibbs Free Energy into a log-likelihood ratio, likelihood estimates for promoter and non-promoter regions
     are passed into the function, the average Gibbs Free Energy of a region around the -10 motif, as dictated by the user,
     is taken (using a double averaging method), and the average value is used to determine the frequency of said value 
     in both the promoter and non-promoter distributions. Dividing these frequencies and taking the logarithm 
     successfully converts Gibbs Free Energy into an LLR value. 
     
    Inputs:
     mySeq: Seq - DNA sequence to be iterated across to find hits
     lmot: Motif Object - -35 motif 
     lthrs: Int - Score threshold for -35 motif
         Note: If the user selects, left threshold can be calculated as the Phisite promoter hunter does (median + stdev)
     rmot: Motif Object - -10 motif 
     rthrs: Int - Score threshold for -10 motif
         Note: If the user selects, right threshold can be calculated as the Phisite promoter hunter does (median + stdev)
     GC_null_hypo: Boolean - Determines whether or not GC content of the sequence is used in the background for PSSM scores
     GC_pseudocnt: Boolean - Determines whether or not GC content of the sequence is used to determine pseudocount. 
     minD: Int - Minimum spacer value
     maxD: Int - Maximum spacer value
     use_GibbsFE: Boolean - True to incoporate GibbsFE into promoter hunter. 
     mode: str - LLR/norm | Used to dictate the mode of summing the contributions of the left/right motif and GibbsFE
     lerg: Int - Distance to the left of the start of the -10 sequence where unzipping begins and the Gibbs Free Energy Likelihood is evaluated across 
     rerg: Int - Distance to the right of the -10 sequence where unzipping begins and the Gibbs Free Energy Likelihood is evaluated across 
     wsize: Int - Window used to average the Gibbs Free Energy vector
     GibbsFEDict: Dictionary - Dictionary containing keys (dinucleotides/trinucleotides) and values (Gibbs Free Energy values)
         Used to convert nucleotides into a Gibbs Free Energy vector
     posdist: np Array - Contains Gibbs Free Energy likelihood estimate of promoter set
     negdist: np Array - Contains Gibbs Free Energy likelihood estimate of negative set
     binedges: np Array - Contains bin ranges of positive and negative arrays 
    
    Returns:
     hmsp: List - List of promoter sequences with -10 and -35 sequences that surpass the motif thresholds. List composed of dictionaries with information on each identified promoter. 
    '''
    #GC = gc(mySeq)/100
    
    #Used to remove GibbsFE contribution from final score
    disregard_GibbsFE = False
    
    
    #Calculates GC content if necessary for background or pseudocount
    if GC_null_hypo == True or GC_pseudocnt == True:
        
        #Initializes dictionary with frequencies of base pairs
        GC_cont_dict = Det_Nuc_Freq(mySeq)
        
        #Sets the background based on GC content
        if GC_null_hypo == True:
            
            #Sets the background of the right motif
            rmot.background = GC_cont_dict
            
            #Sets the background of the left motif
            lmot.background = GC_cont_dict
        
        #Sets the pseudocount value based on GC content
        if GC_pseudocnt == True:
            
            #Sets the pseudocount of the right motif
            rmot._pseudocounts = GC_cont_dict
            
            #Sets the pseudocount of the left motif
            lmot._pseudocounts = GC_cont_dict
        

        
    #score sequence with left motif
    lscrs=lmot.pssm.calculate(mySeq).tolist()
    
    #score sequence with right motif
    rscrs=rmot.pssm.calculate(mySeq).tolist()
    
    #Changes and Computes left threshold used in original pHunt if set by user
    if lthrs == 'original':
        lmed = st.median(lscrs)
        lstd = st.stdev(lscrs)
        lthrs = lmed + lstd
           
    #Changes and Computes right threshold used in original pHunt if set by user
    if rthrs == 'original':
        rmed = st.median(rscrs)
        rstd = st.stdev(rscrs)
        rthrs = rmed + rstd
    
    #Completes the following if using GibbsFE
    if use_GibbsFE == True:
        
        #Computes shift variable as half of window size
        #Used to move the range by which the second average is taken left
        shift = int(wsize/2)
        
        #compute the moving average Gibbs Free Energy vector on the sequence
        fes = GibbsFE(mySeq, GibbsFEDict,wsize)
    
    #list of spacer-adequate high-scoring motif pairs
    hsmp = []
    
    mySeqlen = len(mySeq)
    
    #go through sequence, first with left PSSM. Iterate through sequence with consideration for length of minimum spacer and both motifs. 
    for pl in range(mySeqlen - (lmot.pssm.length + minD + rmot.pssm.length)):
        
        if lscrs[pl] > lthrs:
            remrange = pl+lmot.pssm.length+maxD+1 if (pl+lmot.pssm.length+maxD+1 < mySeqlen - rmot.pssm.length) else mySeqlen - rmot.pssm.length #else statement for potential promoters close to the end of the sequence
            
            #go through sequence, now with right PSSM. Considers placement of right motif depending on any possible permutation off spacer length
            for pr in range(pl+lmot.pssm.length+minD,remrange):
                
                if rscrs[pr] > rthrs:
                    #Initializes information dictionary
                    element = {'lpos' : pl, 'lseq' : str(mySeq[pl:pl+lmot.pssm.length]),
                               'lscr' : lscrs[pl],
                               'rpos' : pr, 'rseq' : str(mySeq[pr:pr+rmot.pssm.length]),
                               'rscr' : rscrs[pr],
                               'spcr' : str(mySeq[pl+lmot.pssm.length:pr].lower()),
                              }
                    
                    #compute overall PSSM contribution
                    element['mscr'] = element['lscr'] + element['rscr']
                    
                    if use_GibbsFE == True:
                        #Initalizes the ranges as positions from which an average of a moving average Gibbs Free Energy vector will be calculated
                        #Shifts the range to the left by half the wsize so as to account for the first averaging
                        #Prevents lrange or range from being out of bounds (extending beyond the averaged vector)
                        lrange = pr - lerg - shift + len(rmot) if pr - lerg - shift + len(rmot) > 0 else 0
                        rrange = pr + rerg - shift + len(rmot) if pr + rerg - shift + len(rmot) < len(fes) - 1 else len(fes) - 1
                    
                        #Issue when both lrange and rrange are out of range of the averaged vector
                        #Prepares to raise a warning and disregard GibbsFE contribution in these cases
                        if lrange > len(fes) - 1 or rrange < 0:
                            warnings.warn('Lrange value is higher than last index of the averaged Gibbs Free Energy vector' + 
                                          ' or Rrange value is lower than first index of the vector. ' +
                                          'Use larger lerg/rerg values or a smaller window size. '
                                          + 'Disregarding Gibbs Free Energy contribution for this inputted sequence')
                            
                            disregard_GibbsFE = True
                    
                    #LLR specific computations
                    if mode == 'LLR':
                        
                        #Averages values from GibbsFE vector if lrange and rrange are in range and the program is using GibbsFE
                        if disregard_GibbsFE == False and use_GibbsFE == True:
                            element['escr'] = sum(fes[lrange:rrange])/(rrange-lrange)
                            p_pos = posdist[np.where(binedges>element['escr'])[0][0] - 1]
                            p_neg = negdist[np.where(binedges>element['escr'])[0][0] - 1]
        
        
                            #Creates the log likelihood for Gibbs Free Energy
                            element['ellr'] = np.log2(p_pos / p_neg)
                        
                        else:
                            #Sets GibbsFE contribution to 0
                            element['escr'] = 0
                            element['ellr'] = 0
                            
                        #compute final score (sum of LLRs)
                        element['Fscr'] = element['mscr'] + element['ellr']
                    
                    #Norm specific computations
                    else:
                        #compute normalized PSSM contribution
                        element['Nlscr'] = (element['lscr'] - lmot.pssm.min) / (lmot.pssm.max - lmot.pssm.min)
                        element['Nrscr'] = (element['rscr'] - rmot.pssm.min) / (rmot.pssm.max - rmot.pssm.min)
                        element['Nmscr'] = element['Nlscr'] + element['Nrscr']
        
                        if disregard_GibbsFE == False and use_GibbsFE == True:
                            #Averages values from GibbsFE vector
                            element['escr'] = sum(fes[lrange:rrange])/(rrange-lrange)
                            
                            #normalizes GibbsFE score
                            element['Nescr'] = (element['escr'] - 9.2) / (3 - 9.2)
                        
                        else:
                            #Sets GibbsFE contribution to 0
                            element['escr'] = 0
                            element['Nescr'] = 0
                        
                        #compute final score (global norm motif score + 2*norm energy score)
                        element['Fscr'] = element['Nmscr'] + (2 * element['Nescr'])
                        

                    #add high-scoring motif pair with all info to the list of HSMPs
                    hsmp.append(element)
                  
    #return list of spacer-adequate high-scoring motif pairs
    return hsmp

def write_file(writer,hit_dict,mySeq,seq_counter,mode_fscr):
    '''
    Writes a row for each hit, listing key/value information from their dictionary. 
    
    After finding suspected hits on a promoter sequence and storing these hits
    and their information within a list of dictionaries, an output file is generated
    listing several of these calculated metrics, including the PSSM scores, GibbsFE score, and 
    final score. 
    
    This function acts to write a row within the csv file by accessing information
    from the dictionary that is passed into the function. By running this function
    on the list of hits generated, each hit's information can be written into the csv file.
    Depending on the final score mode, different metrics are recorded on the csv file.
    
    writer: csvwriter - Permits writing onto the csv file
    hit_dict: dict - Contains keys and values listing the metrics of a hit
    mySeq: list - List of promoter sequences used to record the sequence that a hit comes from
    seq_counter: Int - Integer used to keep track of the input sequence
    mode_fscr: str - LLR or norm used to determine which metrics to write onto the csv
    '''
    
    if mode_fscr == 'LLR':
        writer.writerow([
            mySeq[seq_counter].id,
            hit_dict['rpos'] - hit_dict['lpos'] - len(hit_dict['lseq']),
            str(hit_dict['lpos']) + '-' + str(hit_dict['rpos'] + len(hit_dict['rseq'])),
            hit_dict['lscr'],
            hit_dict['rscr'], 
            hit_dict['escr'],
            hit_dict['ellr'],
            hit_dict['Fscr'],
            hit_dict['lseq']+hit_dict['spcr']+hit_dict['rseq'],
            ])
    else:
        writer.writerow([
            mySeq[seq_counter].id,
            hit_dict['rpos'] - hit_dict['lpos'] - len(hit_dict['lseq']),
            str(hit_dict['lpos']) + '-' + str(hit_dict['rpos'] + len(hit_dict['rseq'])),
            hit_dict['lscr'],
            hit_dict['rscr'], 
            hit_dict['escr'],
            hit_dict['Fscr'],
            hit_dict['lseq']+hit_dict['spcr']+hit_dict['rseq'],
            ])

def remove_dupl(hit_list):
    '''
    Removes duplicate hits in a list of hits for a promoter segment
    
    When promoter hunting, duplicate hits, defined as two hits with the same -35 or -10 position, appear.
    Only one of these locations can be a true promoter, so the duplicate hit with a lower final score is removed.
    
      
    This function acts to address these duplicate hits by iterating through a 
    list of hits for a promoter segment sorted by final score in descending order,
    checking if the right position is present within an initialized list of right positions,
    and either labeling the hit as a non-duplicate if it is not within the list and 
    adding the right position to the list or labeling the hit as a duplicate if it is 
    found within the list. The list of hits is then iterated through again
    in order to check for duplicates in the left position. 
    Duplicates are removed from the list and the list of dictionaries is returned
    list 

    Parameters
    ----------
    hit_list : list - List of hits for a sequence

    Returns
    -------
    mypromoters_edited : list - List of hits for a sequence without duplicates

    '''
    #Initializes left and right position arrays holding original positions
    first_position_array = []
    second_position_array = []
    
    #List of hits without duplicates
    mypromoters_edited = []
    
    #Iterates through list, beginning with right position
    for promoters in hit_list:
        promoters['non-Duplicate'] = True
        
        #Labels duplicate
        if promoters['rpos'] in second_position_array:
            promoters['non-Duplicate'] = False
        
        else:
            second_position_array.append(promoters['rpos'])
            
        #Repeats for left position
        if promoters['lpos'] in first_position_array:
            promoters['non-Duplicate'] = False

        else:
            first_position_array.append(promoters['lpos'])
    
    #Adds non-duplicates to list and returns
    for promoters in hit_list:
        if promoters['non-Duplicate']:
            mypromoters_edited.append(promoters)
            
    return mypromoters_edited

def file_format(file,extension):
    '''
    Helper function used to make sure a file has a necessary extension/is a particular format
    
    The function takes in an extension string, counts the length of this string, and
    checks the end of the file is the same extension string (proving necessary format). 
    
    Parameters
    ------------
    file: str - Name of file to be checked
    extension: str - Type of extension (csv, fas, fasta)
    
    Returns
    ------------
    Boolean: True if it has the proper extension and False if it does not
    '''
    
    extension = '.' + extension
    
    #Checks length of extension to be used to check the end of the file
    extension_len = len(extension)
    
    #Checks the end of the file string to see if the extension is the correct one
    if file[(-1 * extension_len):] == extension:
        return True
    else:
        return False
    
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
        
    #Parses -35 and -10 files for motifs or sequences. 
    #Creates and stores motifs in variables
    if '.fas' == json_file["left_motif"][-4:] or '.fasta' == json_file["left_motif"][-6:]:
        l_insta = SeqIO.parse("../data/" + json_file["left_motif"],'fasta')
        
        #Creates the -35 motif
        lmot = motifs.create([l.seq for l in l_insta],'ACGT')  
        
    elif '.jaspar' == json_file["left_motif"][-7:]:
        #Creates the -35 motif
        lmot = motifs.read(open("../data/" + json_file["left_motif"]),'jaspar') 
    
    #Non-valid motif file format
    else:
        raise(Exception('Non-valid motif file format. Please use a fasta file or jaspar file for the -35 motif'))
     

    if '.fas' == json_file["right_motif"][-4:] or '.fasta' == json_file["right_motif"][-6:]:   
        r_insta = SeqIO.parse("../data/" + json_file["right_motif"],'fasta')
        
        #Creates the -10 motif
        rmot = motifs.create([r.seq for r in r_insta],'ACGT')
        
    elif '.jaspar' == json_file["right_motif"][-7:]: 
        #Creates the -10 motif
        rmot = motifs.read(open("../data/" + json_file["right_motif"]),'jaspar') 
    
    #Non-valid motif file format
    else:
        raise(Exception('Non-valid motif file format. Please use a fasta file or jaspar file for the -10 motif'))
    
    #Initializes variable to determine whether pseudocount is dependent on sequence base pair frequency
    GC_pseudocnt = json_file["use_GCcont_pseudocnt"]
    
    #Assigns individual pseudocounts if the above statement is false
    if GC_pseudocnt == False:
        
        #Assigns left motif psuedocount based on values in json
        lmot._pseudocounts = json_file["non-GCcont_pseudocount"]["left_motif_pseudocounts"]
        
        #Assigns right motif psuedocount based on values in json
        rmot._pseudocounts = json_file["non-GCcont_pseudocount"]["right_motif_pseudocounts"]
    
    #Assigns null hypothesis for PSSM score to variable (GC content based or random)
    null_hypo = json_file["use_GCcont_background"]
    
    #Assigns threshold type (patser or original) from json
    threshold_type = json_file["motif_threshold"]
    
    #Checks to see whether the user inputs an acceptable threshold option
    #Calculates rmot and lmot patser thresholds if the user inputs patser
    if threshold_type == 'patser':
        
        rthrs = rmot.pssm.distribution(precision=10**3).threshold_patser()
        
        lthrs = lmot.pssm.distribution(precision=10**3).threshold_patser()
    
    #Assigns rthrs and lthrs to string value
    #Original phisite promoter hunter thresholds are dependent on the sequence
    elif threshold_type == 'original':
        
        rthrs = 'original'
        
        lthrs = 'original'
    
    #Raises an exception if the threshold_type is not original or patser
    else:
        raise Exception("You may only input patser or original for the \"motif threshold\" in the json file")
    
        
    #Accesses minimum and maximum spacer values from json_file (settings.json) 
    minD = json_file["min_spacer_value"]
    maxD = json_file["max_spacer_value"]
    
    #Reads whether the user wants to use Gibbs Free Energy in the hunt for promoters
    use_GibbsFE = json_file["use_GibbsFE"]
    

    #Prepares Gibbs Free Energy related parameters if the user chooses to use GibbsFE in the hunt for promoters
    if use_GibbsFE == True:
        
        #Distance around start of -10 sequence where unzipping begins. 
        #Necessary for the second averaging of Gibbs Free Energy. 
        #Lerg: Amount of base pairs to the left of the -10 sequence
        #Rerg: Amount of base pairs to the right of the -10 sequence 
        lerg = json_file["GibbsFE_related_parameters"]["lerg"]
        rerg = json_file["GibbsFE_related_parameters"]["rerg"]
    
        #Total amount of base pairs initially unzipped. 
        #Used in pos and neg sets as window over which average Gibbs Free Energy is computed
        unzip_dist = lerg + rerg
       
        #Only permits even window sizes for pHunt function
        if json_file["GibbsFE_related_parameters"]["GibbsFE_windowsize"] % 2 == 0:
            #Initializes the window size used to average Gibbs Free Energy of the inputted sequence
            wsize = json_file["GibbsFE_related_parameters"]["GibbsFE_windowsize"]
            
        else:
            raise Exception('Window size for averaging Gibbs FE vector must be even')
        
        #Initializes the GibbsFE Dictionary from csv file in settings.json
        #Contains GibbsFE values for dinucleotides/trinucleotides
        GibbsFEdict = GibbsFE_Dict(json_file["GibbsFE_related_parameters"]["GibbsFE_file"])
      
    #Intializes the user inputted Final Score computting technique 
    #Computes final score using an LLR summation or a sum of normalized values
    mode_fscr = json_file["mode_fscr"]
    
    #Parses the input_sequences file stored in the json file
    #Raises a warning if the file is not a fasta file
    if '.fas' == json_file["input_sequences"][-4:] or '.fasta' == json_file["input_sequences"][-6:]:
        mySeq = list(SeqIO.parse("../data/" + json_file["input_sequences"],'fasta'))
    
    else:
        raise Exception('You must use a fasta file (.fas or .fasta) for "input_sequences" in the json_file')

    #Initializes variable with output file name  
    #Raises an error if the output file's name is not a csv file
    if json_file["output_information"][-4:] == '.csv':
        output_file = json_file["output_information"] 
        
    else:
        raise Exception('Output file name, inputted in "output_information" in the json file, must be a .csv file')
        
    #Prepares LLR specific parameters and runs pHunt
    if mode_fscr == "LLR":
        
        #Only prepares the positive likelihood estimate if the program is using GibbsFE
        if use_GibbsFE == True:
            #Parses the positive set and generates likelihood array if the input is a fasta file
            if '.fas' == json_file["LLR_specific_parameters"]["promoter_set"][-4:] or '.fasta' == json_file["LLR_specific_parameters"]["promoter_set"][-6:]:   
                #Reads file
                seq_objectp = list(SeqIO.parse("../data/" + json_file["LLR_specific_parameters"]["promoter_set"],'fasta'))
                
                #Creates likelihood array
                #Uses bin amount and pseudocount values as specified in json file
                posdist,bin_edges = GibbsFE_likelihood(seq_objectp, unzip_dist, GibbsFEdict, json_file["LLR_specific_parameters"]["GibbsFE_n_bins"], json_file["LLR_specific_parameters"]["GibbsFE_pseudocounts"])
            
            else:
                #Raises an error if the file extension is invalid
                raise Exception('Please use a csv file or fas/fasta file for promoter set')
                
        
        #Prepares the background set only if the user is using GibbsFE
        if use_GibbsFE == True:
            
            #Parses the background set and generates likelihood array if the input is a fasta file
            if '.fas' == json_file["LLR_specific_parameters"]["background_set"][-4:] or '.fasta' == json_file["LLR_specific_parameters"]["background_set"][-6:]:   
                
                #Reads file
                seq_objectb = list(SeqIO.parse("../data/" + json_file["LLR_specific_parameters"]["background_set"],'fasta'))
           
                #If the negative set is a whole genome, 
                #the genome will be divided into multiple segments. 
                #Each segment will be then be ran into the GibbsFE_Likelihood function
                if len(seq_objectb) == 1:  
                    
                    #Scan_wind is the length of the segments that the genome will be divided into. 
                    #Calculated as the average length of the promoter sequences 
                    w = [len(seq) for seq in seq_objectp]
                    scan_wind= round(sum(w)/len(w))
                    
                    #Divides the genome set into segments and appends these segments to list
                    for q in range(0,len(seq_objectb[0].seq)-scan_wind + 1,json_file["LLR_specific_parameters"]["step_size_genome"]):
                        seq_objectb.append(seq_objectb[0][q:q+scan_wind])
                    
                    #Removes the whole genome entry
                    seq_objectb.pop(0)
            
                #Generates likelihood array
                negdist,bin_edges = GibbsFE_likelihood(seq_objectb, unzip_dist, GibbsFEdict, json_file["LLR_specific_parameters"]["GibbsFE_n_bins"], json_file["LLR_specific_parameters"]["GibbsFE_pseudocounts"]) 
                
            #Reads from file if user has stored lists in csv file
            elif '.csv' == json_file["LLR_specific_parameters"]["background_set"][-4:]:
                
                #Prepares to get negdist and bin_edges array from file
                with open('../data/' + json_file["LLR_specific_parameters"]["background_set"],newline = '') as b:
                    listreaderb = csv.reader(b)
                    
                    #Initializes lists for negdist and bin_edges
                    listreaderb = list(listreaderb)
                    negdist = listreaderb[0]
                    bin_edges = listreaderb[1]
                    
                    #Iterates through both lists to convert strings to floats
                    for indn in range(len(negdist)):
                        negdist[indn] = float(negdist[indn])
                        
                    for indb in range(len(bin_edges)):
                        bin_edges[indb] = float(bin_edges[indb])
                    
                    #Converts lists in numpy arrays
                    negdist = np.array(negdist)
                    bin_edges = np.array(bin_edges)
            
            #Raises Error if file extension is incorrect
            else:
                raise Exception('Please use a csv file or fas/fasta file for the background set')
            
        #Conducts promoter analysis
        seq_hits = []
        
        for sequence in mySeq:
            
            #If else statement removing parameters if the program is PSSM only
            if use_GibbsFE == True:
                #Stores hits for each promoter as a list of dictionaries
                hit = pHunt(sequence.seq, lmot, lthrs, rmot, rthrs, null_hypo, GC_pseudocnt, minD, maxD, use_GibbsFE, 'LLR', lerg, rerg, wsize, GibbsFEdict, posdist, negdist, bin_edges)
            
            else:
                hit = pHunt(sequence.seq, lmot, lthrs, rmot, rthrs, null_hypo, GC_pseudocnt, minD, maxD, use_GibbsFE, 'LLR')
            
            #Sorts list by score in descending order
            hit_sorted = sorted(hit, key = lambda k: k['Fscr'], reverse=True)
            
            #Removes duplicates
            hit_sorted = remove_dupl(hit_sorted)
            
            #List of hits for each sequence inputted            
            seq_hits.append(hit_sorted)
            
    #Generates hits using original pHunt techniques (normalized)
    elif mode_fscr == 'norm':
        
        seq_hits = []
        
        #Iterates through input_sequences to promoter hunter on each sequence inputted. 
        for sequence in mySeq:
            
            #If else statement removing variable inputs if the program is PSSM only
            if use_GibbsFE == True:
                #Stores hits for each promoter as a list of dictionaries
                hit = pHunt(sequence.seq, lmot, lthrs, rmot, rthrs, null_hypo, GC_pseudocnt, minD, maxD, use_GibbsFE, 'norm', lerg, rerg, wsize, GibbsFEdict)
            
            else:
                hit = pHunt(sequence.seq, lmot, lthrs, rmot, rthrs, null_hypo, GC_pseudocnt, minD, maxD, use_GibbsFE, 'norm')
                
            #Sorts hits based on final score
            hit_sorted = sorted(hit, key = lambda k: k['Fscr'], reverse=True)
            
            #Removes duplicates of sorted hit list
            hit_sorted = remove_dupl(hit_sorted)
            
            
            #Adds to list of promoter hits
            seq_hits.append(hit_sorted)
         
    #Raises an error if the user does not input LLR or norm
    else:
        raise Exception("You may only input LLR or norm for \"mode_fscr\" in the json file")

    #Writes output file        
    with open('../data/' + output_file,'w',newline = '') as f:
        
        #Initializes csv writer
        writer = csv.writer(f)
        
        #Initializes variables to record the sequence a hit corresponds to 
        mySeq = list(mySeq)
        seq_counter = 0
        
        #Initializes variable to determine what type of hits are written in the csv file
        output_type = json_file["output_type"]
        
        #Writes input name, mode, and date/time
        writer.writerow(['output file: ' + output_file,'Mode Final Score: ' + mode_fscr,'Time of running: ' + ctime()])
        
        #LLR specific recorded metrics
        if mode_fscr == 'LLR':
            
            #Initializes titles for column
            writer.writerow(['ID','Spacer_length','Range','Left_Motif_Score','Right_Motif_Score','Average_GibbsFE','GibbsFE_Log-Likelihood','Final_Score', 'Promoter_Sequence'])
            
            #Iterates through list of hits for all sequences to get the list of hits for each sequence
            for prom_hit in seq_hits:
                

                #Skips a sequence if it does not have hits
                if len(prom_hit) > 0:
                    
                    #Writes only the top hits for each sequence
                    if output_type == 'top hits':
                        write_file(writer,prom_hit[0],mySeq,seq_counter,'LLR')
                    
                    #Writes all hits for each sequence
                    elif output_type == 'all hits':
                    
                        #Iterates through hits to write each hit
                        for hit in prom_hit:
                            write_file(writer,hit,mySeq,seq_counter,'LLR')
                    
                    else:
                        raise Exception('Output type, inputted in "output_type" in the json file, may only be all hits or top hits')
                seq_counter += 1
        
        #Writing for normalized mode
        else:
            
            #norm specific recorded metrics
            writer.writerow(['ID','Spacer_length','Range','Left_Motif_Score','Right_Motif_Score','Average_GibbsFE','Final_Score', 'Promoter_Sequence'])
            
            #Iterates through list of hits for all sequences to get the list of hits for each sequence
            for prom_hit in seq_hits:
                
                #Skips a sequence if it does not have hits
                if len(prom_hit) > 0:
                    
                    #Writes only the top hits for each sequence
                    if output_type == 'top hits':
                        write_file(writer,prom_hit[0],mySeq,seq_counter,'norm')
                    
                    #Writes all hits for each sequence
                    elif output_type == 'all hits':
                
                        #Iterates through hits to write each hit
                        for hit in prom_hit:

                            write_file(writer,hit,mySeq,seq_counter,'norm')
    
                    else:
                        raise Exception('Output type, inputted in "output_type" in the json file, may only be all hits or top hits')
                
                seq_counter += 1
            
    f.close()

#Creates output csv file upon running the Python file
if __name__ == "__main__":
    go()
    



