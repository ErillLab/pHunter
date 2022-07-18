# Settings.json File overview

The Settings file of the program is a .json, meaning variables are stored using a java script dictionary and accessed through the json library in python. The file allows for user input in the program. The text below seeks to clarify the meaning and use of the keys within the json dictionary.

## Positive and Negative Set
- Several parameters that control the creation of the positive and negative likelihood-estimate arrays
### promoter_set
- Type: `string`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas, .fasta)
- Default value: `PromEC_seqs_filtered.fas`
- Rationale/Format:
	- The promoter set acts as the alternative hypothesis (positive set) in this model. The positive set, which is stored as a likelihood-estimate, is compared to the background set through a log-likelihood to determine the likelihood by which a sequence of a particular average Gibbs Free Energy value satisfies the alternative hypothesis over the null.
	- Promoter sequences are held within a fasta file for creation of the positive set. The inputted file is parsed, evaluated for dinucleotides/trinucleotides, and translated into a Gibbs Free Energy vector based on Gibbs Free Energy values of the nucleotides that compose the promoter sequences (Gibbs Free Energy values of dinucleotides/trinucleotides determined by `GibbsFE_file`). Information from this Gibbs Free Energy vector is used to construct a positive set likelihood-estimate.

### background_set
- Type: `string`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas, .fasta)
- Default value: `NC_000913.3.fasta`
- Rationale/Format:
	- The background set acts as the null hypothesis (negative set) in this model. The negative set, which is stored as a likelihood-estimate, is compared to the promoter set through a log-likelihood to determine the likelihood by which a sequence of a particular average Gibbs Free Energy value satisfies the alternative hypothesis over the null.
	- Genome sequences are held within a fasta file for creation of the negative set. The inputted file is parsed, evaluated for dinucleotides/trinucleotides, and translated into a Gibbs Free Energy vector based on Gibbs Free Energy values of the nucleotides that compose the genome (Gibbs Free Energy values of dinucleotides/trinucleotides determined by `GibbsFE_file`). Information from this Gibbs Free Energy vector is used to construct a negative set likelihood-estimate.
	- The background set can either be the entire genome of a particular species - for the use case of finding any potential promoters in the genome - or intergenic sequences of a particular species - for the use case of finding the potential promoters that are most likely true positives.
		- Whole genome fasta

### step_size_genome
- Type: `Integer`
- Allowed values:
	- Any positive Integer value. Recommended to be the average length of the promoters in the positive set or larger.
- Default value: 101
- Rationale/Format:
	- The negative and positive conditional distributions (derived from the set given a certain Gibbs Free Energy value) are compared to each other in the Log-likelihood ratio. In order to be compared to each other properly, both sets need to follow a similar format. In the case that the background set is the entire genome of a particular species, the background set must be edited to be properly compared to the positive set.
	- If the background set is an entire genome (json_file["background_set"] has one stored sequence), this sequence is iterated on, dividing the sequence into a list of segments (size = average length of the promoters in the positive set). This list of segments is then passed to GibbsFE_likelihood, creating a negative likelihood-estimate that is comparable to the positive likelihood-estimate.
	- step_size_genome = step size of the segment windows. If step_size_genome = average length of the promoters in the positive set, all the segments are contiguous.

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

### GibbsFE_n_bins
- Type: `Integer`
- Allowed values:
	- Any integer above 1
- Default value: `30`
- Rationale/Format:
	- In order to convert the Gibbs Free Energy vector into a likelihood-estimate, binning is necessary to create an array of counts. The GibbsFE_n_bins controls the amount of bins created for the likelihood-estimate array for both the positive and negative sets.
	- The np.histogram package uses GibbsFE_n_bins, the Gibbs Free Energy vector, and the minimum and maximum Gibbs Free Energy values to create an array of counts and an array of bin edges.

### GibbsFE_pseudocounts
- Type: `Integer`
- Allowed values:
	- Any integer above 0
- Default value: 1
- Rationale/Format:
	- When converting the Gibbs Free Energy vector into an array of counts with bin edges, several bins may possess a count of 0, an incompatible number with a log-likelihood ratio. To prevent the creation of a frequency array with values equal to 0 and to bridge the gap between relative frequency and probability, a pseudocount is added to each value in the counts array.
	- When converting the counts array into a frequency array, GibbsFE_pseudocounts is added to each element in the counts array.

## Putative Promoter Search
- Parameters involved in the search for potential parameters
### input_sequences
- Type: `string`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas, .fasta)
- Default value: `Input_Sequences.fas`
- Rationale/Format:
	- The sequences within the input_sequences file are iterated on to predict putative promoters. Using a strict spacer range (from `min_spacer_value` to `max_spacer_value`) between the consensus sequence, potential locations where the -35 and -10 pssms return a score above the threshold inputted by the user(`left_motif_threshpres` and `right_motif_threshpres`) are reported as putative promoters in the output file. Among these putative promoters, average Gibbs Free Energy of a distance (`sampling_size`) around the start of the -10 sequence is used as to create a conditional distribution when comparing the likelihoods of the positive and negative sets within the Log-Likelihood Ratio.
	- The input_sequences file contains a Seq Object with either a single sequences or multiple sequences. This object is iterated on so that each sequence within the file may be checked for putative promoters. Iterating on each sequence, locations with -35 sequences that score above the threshold set by the user are found, and potential -10 sequences are searched for any distance from the `min_spacer_value` to the `max_spacer_value` base pairs away from the -35 sequence. If the -10 sequence scores above the threshold set, it is marked as a promoter.

### left_motif
- Type: `String`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas,.fasta)
- Default value: `Eco-35.fas`
- Rationale/Format:
	- Fasta file containing instances of the left motif used to search for -35 sequences in each sequence in the `input_sequences` file.

### left_motif_threshprec
- Type: `float`
- Allowed values:
	- Any value above 0
- Default value: 0.03828859126248574 (Calculated using the patser threshold within Biopython)
- Rationale/Format:
	- In order to determine whether or not a location is a -35 sequence, the -35 pssm must score a sequence, and a score threshold must be set for putative binding sites. The threshold set by `left_motif_threshprec` is used by the Program to determine what qualifies as a -35 sequence.

### left_motif_pseudocounts
- Type: `Float`
- Allowed values:
	- Any float value above 0.0
- Default value: 0.05
- Rationale/Format:
	- Pseudocount value applied to the Position specific-scoring matrix of the left motif to prevent any base at any position in the matrix from yielding a negative infinity value.
	- Float value assigned to ._pseudocounts propery of the left motif object.

### right_motif
- Type: `String`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas,.fasta)
- Default value: `Eco-10.fas`
- Rationale/Format:
	- Fasta file containing instances of the right motif used to search for -10 sequences in each sequence in the `input_sequences` file.

### right_motif_threshprec
- Type: `float`
- Allowed values:
	- Any value above 0
- Default value: 1.9342901670836667 (Calculated using the patser threshold within Biopython)
- Rationale/Format:
	- In order to determine whether or not a location is a -10 sequence, the -10 pssm must score a sequence, and a score threshold must be set for putative binding sites. The threshold set by `right_motif_threshprec` is used by the program to determine what qualifies as a -10 sequence.

### right_motif_pseudocounts
- Type: `Float`
- Allowed values:
	- Any float value above 0.0
- Default value: 0.05
- Rationale/Format:
	- Pseudocount value applied to the Position specific-scoring matrix of the right motif to prevent any base at any position in the matrix from yielding a negative infinity value.
	- Float value assigned to ._pseudocounts propery of the right motif object.

### min_spacer_value
- Type: `Integer`
- Allowed values:
	- Any integer value above 0
- Default value: 15
- Rationale/Format:
	- Minimum spacer value between the -35 and -10 consensus sequences. Used in pHunt function when evaluating sequences in `input_sequences` to test possible locations of -10 sequence considering a -35 sequence has already been identified.
	- Integer value setting the closest the -10 sequence could be to the -35 sequence. Potential locations for start of -10 sequence are between `min_spacer_values` to `max_spacer_values` base pairs away from the -35 sequence.

### max_spacer_value
- Type: `Integer`
- Allowed values:
	- Any integer value equal to or greater than min_spacer_value
- Default value: 18
- Rationale/Format:
	- Maximum spacer value between the -35 and -10 consensus sequences. Used in pHunt function when evaluating sequences in `input_sequences` to test possible locations of -10 sequence considering a -35 sequence has already been identified.
	- Integer value setting the furthest the -10 sequence could be to the -35 sequence. Potential locations for start of -10 sequence are between `min_spacer_values` to `max_spacer_values` base pairs away from the -35 sequence.

### sampling_size
- Type: `Integer`
- Allowed values:
	- Any integer above 2
- Default value: 25
- Rationale/Format:
	- To begin transcription, a segment of base pairs around the -10 sequence are unzipped at one time to assist the unzipping of additional sequences. This sequence that is initially unzipped is AT rich (Low in Gibbs Free Energy meaning easy to unzip). sampling_size represents the length of this sequence within putative promoters.
	- sampling_size is the area around the start of the -10 sequence where unzipping initiates. The value is used to calculate a left range(-10_start_location - `unzip_dist` // 2 ) and right range (-10_start_location + `unzip_dist` // 2) from the start of the -10 sequence. An average of a Gibbs Free Energy vector of the sequence from these ranges is used to create a conditional distribution for the positive and negative sets within the Gibbs Free Energy Log Likelihood Ratio.
	- sampling_size is also used as the window size for the positive and negative sets within the `GibbsFE_likelihood` function. The positive and negative sequences are converted into Gibbs Free Energy values and are averaged over the window. Using this strategy, the positive and negative arrays are directly comparable to the value obtained when averaging the Gibbs Free Energy vector around the -10 sequence.

## Information Output
- Parameter involved in the communication of findings of putative promoters.
### output_information
- Type: `string`
- Allowed values:
	- Any file name that follows csv naming conventions (.csv)
- Default value: `Output_Information.csv`
- Rationale:
	- Output_Information is the csv file containing information about each promoter hit. After a potential promoter sequence passes both left and right motif thresholds, a dictionary is initialized for the sequence, containing the -35 sequence, spacer nucleotides, -10 sequence, PSSM scores, Gibbs Free Energy Log Likelihood Ratio, and Final score. Dictionary is accessed and written in the output file by the csv.writer package.
	- Csv file written by the csv.writer package. Contains the following headers: ID, Spacer_length, Range, Left_Motif_score, Right_Motif_Score, Average_GibbsFE, GibbsFE_Log-Likelihood, Final_Score, and Promoter_Sequence. Sorts promoters by Final_Score.
		- `ID`: ID of Input Sequence for which the promoter is generated for
		- `Spacer_length`: Length of the spacer between the -35 and -10
		- `Range`: Beginning of -35 sequence to end of -10 sequence
		- `Left_Motif_score`: Score of -35 sequence by -35 pssm
		- `Right_Motif_Score`: Score of -10 sequence by -10 pssm
		- `Average_GibbsFE`: Gibbs Free Energy of region around the beginning of the -10 sequence. Used to create a conditional distribution for the GibbsFE_Log-Likelihood Ratio
		- `GibbsFE_Log-Likelihood`: Log-Likelihood of promoter and background sets given the Average_GibbsFE of the initially unzipped region
		- `Final_Score`: Sum of the Left_Motif_score, Right_Motif_Score, and GibbsFE_Log-Likelihood ratio.
		- `Promoter_Sequence`: Base pairs of -35 sequence, spacer, and -10 sequence
