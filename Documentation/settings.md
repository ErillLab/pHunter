# Settings.json File overview

The `settings.json` file is json, meaning variables are stored using a java script dictionary and accessed through the json library in python. The file allows for user input in the program. The text below seeks to clarify the meaning and use of the keys within the json dictionary.

## mode_fscr
- Type: `string`
- Allowed values:
	- 'LLR' or 'norm'
- Default value: 'LLR'
- Rationale/format:
	- While the PhiSite Promoter Hunter computes the Final Score of a putative promoter site as the summation of the normalized values of two PSSMs and the surrounding sequences' Gibbs Free Energy contribution, **pHunt** alows the user to sum PSSM and Gibbs Free Energy contributions together non-arbitrarily by converting the Gibbs Free Energy metric into a Log-Likelihood Ratio
	- By setting this parameter to LLR, the program will convert Gibbs Free Energy into a log-likelihood ratio by taking in a promoter and background set of sequences, constructing a positive likelihood estimate and a negative likelihood estimate, and using the average Gibbs Free Energy around the -10 sequence of a hit to convert these distributions into frequencies. By taking the logarithm of the ratio of these frequencies, the PSSM scores and Gibbs Free Energy can be summed together as LLRs.
	- By setting this parameter to norm, the program will forego the code to convert Gibbs Free Energy to a log-likelihood ratio and instead normalize the average Gibbs Free Energy value around the -10 sequence of a hit. By normalizing the PSSM scores as well and using a coefficient of 2 for the normalized Gibbs Free Energy value, final score is computed as the sum of normalized values

## LLR_specific_parameters
The following parameters are only required if the user enters 'LLR' for `mode_fscr`

### Positive and Negative Set
- Several parameters that control the creation of the positive and negative likelihood-estimate arrays

#### promoter_set
- Type: `string`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas, .fasta)
- Default value: `PromEC_seqs_filtered.fas`
- Rationale/Format:
	- The promoter set acts as the alternative hypothesis (positive set) in this model. The positive set, which is converted into a likelihood-estimate, is compared to the background set through a log-likelihood ratio to determine the likelihood by which a sequence of a particular average Gibbs Free Energy value satisfies the alternative hypothesis over the null.
	- Promoter sequences are held within a fasta file for creation of the positive set. The inputted file is parsed, sequences are evaluated for dinucleotides/trinucleotides, and they are translated into a Gibbs Free Energy vector based on Gibbs Free Energy values stored in  `GibbsFE_file`. Information from this Gibbs Free Energy vector is used to construct a positive set likelihood-estimate.

#### background_set
- Type: `string`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas, .fasta) OR Any csv file containing a likelihood estimate
- Default value: `NC_000913.3.fasta`
- Rationale/Format:
	- The background set acts as the null hypothesis (negative set) in this model. The negative set, which is converted into a likelihood-estimate, is compared to the promoter set through a log-likelihood ratio to determine the likelihood by which a sequence of a particular average Gibbs Free Energy value satisfies the alternative hypothesis over the null.
	- Genome sequences are held within a fasta file for creation of the negative set. The inputted file is parsed, evaluated for dinucleotides/trinucleotides, and translated into a Gibbs Free Energy vector based on Gibbs Free Energy values of the nucleotides that compose the genome (Gibbs Free Energy values of dinucleotides/trinucleotides determined by `GibbsFE_file`). Information from this Gibbs Free Energy vector is used to construct a negative set likelihood-estimate.
	- The background set can either be the entire genome of a particular species - for the use case of finding any potential promoters in the genome - or intergenic sequences of a particular species - for the use case of finding the potential promoters that are most likely true positives.
	- A csv file containing negative set likelihood estimate can also be used to cut down on computation time. The file must contain two Arrays, one for bin counts and one for bin edges. For a reference on format, please check out `negativeset.csv`

### GibbsFE_pseudocounts
- Type: `Integer`
- Allowed values:
	- Any integer above 0
- Default value: 1
- Rationale/Format:
	- When converting the Gibbs Free Energy vector into an array of counts with bin edges, several bins may possess a count of 0, an incompatible number for reference in the log-likelihood ratio. To prevent the creation of a frequency array with values equal to 0 and to bridge the gap between relative frequency and probability, a pseudocount is added to each value in the counts array.
	- When converting the counts array into a frequency array, GibbsFE_pseudocounts is added to each element in the counts array.

### GibbsFE_n_bins
- Type: `Integer`
- Allowed values:
	- Any integer above 1
- Default value: `30`
- Rationale/Format:
	- In order to convert the Gibbs Free Energy vector into a likelihood-estimate, binning is necessary to create an array of counts. The GibbsFE_n_bins controls the amount of bins created for the likelihood-estimate array for both the positive and negative sets.
	- The np.histogram package uses GibbsFE_n_bins, the Gibbs Free Energy vector, and the minimum and maximum Gibbs Free Energy values to create an array of counts and an array of bin edges.

## use_GibbsFE
- Type: `Boolean`
- Allowed values:
	- true or false
- Default value: true
- Rationale/format:
	- In the original PhiSite Promoter Hunter, the user could decide whether or not to use Gibbs Free Energy when hunting for promoters. This functionality has been continuned in **phunt**.
	- By selecting true, the program will incoporate Gibbs Free Energy in when computing Final Score as normal
	- By selecting false, the program will forego code necessary to compute the Gibbs Free Energy contribution of a hit. When the user selects 'LLR' for `mode_fscr`, the code will still forego the computation necessary to convert Gibbs Free Energy into a Log-Likelihood Ratio



## GibbsFE_related_parameters
The following parameters are only required if the user inputs true for `use_GibbsFE`

### GibbsFE_file
- Type: `string`
- Allowed values:
	- Any file name that follows csv naming conventions (.csv)
- Default value: `GibbsFE.csv`
- Rationale/Format:
	- The Gibbs Free Energy dictionary created from the GibbsFE_file acts as a source for Gibbs Free Energy values of particular dinucleotide/trinucleotides combinations. When creating a Gibbs Free Energy vector for a sequence, the program iterates through the inputted sequence and grabs dinucleotides/trinucleotides. By accessing the Gibbs Free Energy dictionary created from the GibbsFE_file, a vector of Gibbs Free Energy values is constructed, assiting the creation of a likelihood-estimate in the case of the positive/negative sets.
	- The Gibbs Free Energy dictionary is held in a headerless csv file. Each line within the csv file has a dinucleotide/trinucleotide (Before the comma) and a Gibbs Free Energy value (After the comma). The program uses the csv.reader package to parse the file and convert the results into a dictionary.
		- `AA,4.18`
		- `AC,6.02`
		- `AG,5.36`
		- `AT,3.68`

### GibbsFE_windowsize
- Type: `Integer`
- Allowed values:
	- Any even integer above 1
- Default value: 50
- Rationale/format:
	- Similar to the PhiSite Promoter Hunter, **pHunt** computes Gibbs Free Energy contribution of the area around the -10 sequence of a hit as the average across a range of the moving average array (a.k.a a double average). In order to take a double average, the moving average array must first be generated. `GibbsFE_windowsize` dictates the size of the window used to compute the moving average arrays
	- `GibbsFE_windowsize` must be even in order for `shift` to be calculated (1/2 of GibbsFE_windowsize). For more information about shift, please consult `user_manual.md`

### lerg
- Type: `Integer`
- Allowed values:
	- Any integer
- Default value: 12
- Rationale/Format:
	- To begin transcription, a segment of base pairs around the -10 sequence is unzipped at one time to assist the unzipping of additional sequences. This sequence that is initially unzipped is AT rich (comparatively easier to unzip). lerg represents the amount of base pairs to the left of the start of the -10 sequence where unzipping likely begins and is used to compute lrange, the left range of the moving average array where the Gibbs Free Energy contribution is computed from.


### rerg
- Type: `Integer`
- Allowed values:
	- Any integer
- Default value: 13
- Rationale/Format:
	- To begin transcription, a segment of base pairs around the -10 sequence is unzipped at one time to assist the unzipping of additional sequences. This sequence that is initially unzipped is AT rich (comparatively easier to unzip). rerg represents the amount of base pairs to the right of the start of the -10 sequence where unzipping likely begins and is used to compute rrange, the right range of the moving average array where the Gibbs Free Energy contribution is computed from.




**Parameters involved in the search for potential parameters**
### input_sequences
- Type: `string`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas, .fasta)
- Default value: `Input_Sequences.fas`
- Rationale/Format:
	- The sequences within the input_sequences file are iterated on to predict putative promoters. Using a strict spacer range (from `min_spacer_value` to `max_spacer_value`) between the consensus sequence, potential locations where the -35 and -10 pssms return a score above the threshold inputted by the user (`motif_threshold`) are reported as putative promoters in the output file. Among these putative promoters, average Gibbs Free Energy of a distance (`sampling_size`) around the start of the -10 sequence is used as to create a conditional distribution when comparing the likelihoods of the positive and negative sets within the Log-Likelihood Ratio.
	- The input_sequences file contains a Seq Object with either a single sequences or multiple sequences. This object is iterated on so that each sequence within the file may be checked for putative promoters. Iterating on each sequence, locations with -35 sequences that score above the threshold set by the user are found, and potential -10 sequences are searched for any distance from the `min_spacer_value` to the `max_spacer_value` base pairs away from the -35 sequence. If the -10 sequence scores above the threshold set, it is marked as a promoter.

### use_GCcont_background
- Type: `Boolean`
- Allowed values:
	- true or false
- Default value: false
- Rationale/Format:
	- When computing the PSSM score for a sequence, `Biopython` uses a null hypothesis that assumes that each base pair's frequency is 0.25 whereas `PhiSite` uses a null hypothesis based on the frequency of each base pair in the input sequence. When set to true, pHunt uses Biopython's null hypothesis. When set to false, pHunt uses PhiSite's null hypothesis.

### left_motif
- Type: `String`
- Allowed values:
	- Any file name that follows fasta or jaspar naming conventions (.fas,.fasta,.jaspar)
- Default value: `Eco-35.fas`
- Rationale/Format:
	- Fasta file containing instances of the left motif used to search for -35 sequences in each sequence in the `input_sequences` file.

### right_motif
- Type: `String`
- Allowed values:
	- Any file name that follows fasta or jaspar naming conventions (.fas,.fasta,.jaspar)
- Default value: `Eco-10.fas`
- Rationale/Format:
	- Fasta file containing instances of the right motif used to search for -10 sequences in each sequence in the `input_sequences` file.

### use_GCcont_pseudocnt
- Type: `Boolean`
- Allowed values:
	- *true* or *false*
- Default value: *false*
- Rationale/Format:
	- Pseudocount values are applied to the alternative hypothesis of the PSSM log-likelihood ratio to prevent the score from going to negative infinity (accounts for the difference between probability and frequency). The PhiSite promoter hunter uses the base pair frequency of the input sequence as a pseudocount value for the PSSM
	- The default value is false since the numerator of the PSSM log-likelihood ratio describes the binding affinity of the transcription factor, which is not beholden to the GC content of the sequence surrounding a binding site.

### non-GCcont_pseudocount
The following two parameters are only required if the input to `use_GCcont_pseudocnt` was false

#### left_motif_pseudocounts
- Type: `Float`
- Allowed values:
	- Any floating value above 0.0
- Default value: 0.25
- Rationale/Format:
	- Pseudocount value applied to the Position specific-scoring matrix of the left motif prevent any base at any position in the matrix from yielding a negative infinity value.
	- Float value assigned to ._pseudocounts propery of the left motif object.

#### right_motif_pseudocounts
- Type: `Float`
- Allowed values:
	- Any float value above 0.0
- Default value: 0.05
- Rationale/Format:
	- Pseudocount value applied to the Position specific-scoring matrix of the right motif prevent any base at any position in the matrix from yielding a negative infinity value.
	- Float value assigned to ._pseudocounts propery of the right motif object.

### motif_threshold
- Type: `string`
- Allowed values:
	- 'original' or 'patser'
- Default value: 'patser'
- Rationale/Format:
	- In order to determine whether or not a location is a transcription factor binding site, the pssm must score a sequence, and a score threshold must be set for putative binding sites. If the user sets `motif_threshold` to 'patser', this threshold is set as the patser threshold computed on the left and right motif using `Biopython`. If the user sets `motif_threshold` to 'original', this threshold is set as the median of PSSM scores + the standard deviation for each sequence (How PhiSite computes their threshold)

### min_spacer_value
- Type: `Integer`
- Allowed values:
	- Any integer value above 0
- Default value: 15
- Rationale/Format:
	- Minimum spacer value between the -35 and -10 consensus sequences. Used in pHunt function when evaluating sequences in `input_sequences` to test possible locations of the -10 sequence considering a -35 sequence has already been identified.
	- Integer value setting the closest the -10 sequence could be to the -35 sequence. Potential locations for start of -10 sequence are between `min_spacer_values` to `max_spacer_values` base pairs away from the -35 sequence.

### max_spacer_value
- Type: `Integer`
- Allowed values:
	- Any integer value equal to or greater than min_spacer_value
- Default value: 18
- Rationale/Format:
	- Maximum spacer value between the the -35 and -10 consensus sequences. Used in pHunt function when evaluating sequences in `input_sequences` to test possible locations of the -10 sequence considering a -35 sequence has already been identified.
	- Integer value setting the furthest the -10 sequence could be to the -35 sequence. Potential locations for start of -10 sequence are between `min_spacer_values` to `max_spacer_values` base pairs away from the -35 sequence.



## Information Output
- Parameter involved in the communication of findings of putative promoters.
### output_information
- Type: `string`
- Allowed values:
	- Any file name that follows csv naming conventions (.csv)
- Default value: `Output_Information.csv`
- Rationale:
	- Output_Information is the csv file containing information about each promoter hit. After a potential promoter sequence passes both left and right motif thresholds, a dictionary is initialized for the sequence, containing the -35 sequence, spacer nucleotides, -10 sequence, PSSM scores, Gibbs Free Energy information, and Final score. Dictionary is accessed and written in the output file by the csv.writer package.
	- Csv file written by the csv.writer package. Contains the headers specified in `user_manual.md` depending on the input for `mode_fscr`

### output_type
- Type: `string`
- Allowed values:
	- 'all hits' or 'top hits'
- Default value: 'all hits'
- Rationale:
	- Once a hit passes the left and right threshold, it is appended to a list of dictionaries with keys and values describing the nature of this hit. If the user selects 'top hits', the output file will only display the highest scoring hit for each input sequence. If the user selects 'all hits', the output file will display every hit passing the left and right threshold for each input sequence.
