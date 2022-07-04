# Settings.json File overview

The Settings file of the program is a .json, meaning variables are stored using a java script dictionary and accessed through the json library in python. The file allows for user input in the program. The text below seeks to clarify the meaning and use of the keys within the json dictionary.

### promoter_set
- Type: `string`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas, .fasta)
- Default value: `PromEC_seqs_filtered.fas`
- Rationale/Format:
	- Promoter sequences are held within a fasta file for creation of the positive set. This positive set acts as a likelihood of Gibbs Free Energy and is compared to the background set. Altering the string allows the user to customize which fasta file has sequences that are used to create what is considered as the positive set. 

### background_set
- Type: `string`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas, .fasta)
- Default value: `NC_000913.3.fasta`
- Rationale/Format:
	- Background sequences are held within a fasta file for creation of the background set. This background set acts as a likelihood of Gibbs Free Energy and is compared to the promoter set. Altering the string allows the user to customize which fasta file has sequences that are used to create what is considered as the background set. 


### GibbsFE_file
- Type: `string`
- Allowed values:
	- Any file name that follows csv naming conventions (.csv)
- Default value: `GibbsFE.csv`
- Rationale/Format:
	- The Gibbs Free Energy dictionary is held in a headerless csv file. Entering the file name of the wanted csv file easily allows for the user to substitute dinucleotide Gibbs Free Energy values of their choosing.

### input_sequences
- Type: `string`
- Allowed values:
	- Any file name that follows fasta naming conventions (.fas, .fasta)
- Default value: `Input_Sequences.fas`
- Rationale

### output_information
- Type: `string`
- Allowed values:
	- Any file name that follows csv naming conventions (.csv)
- Default value: `Output_Information.csv`
- Rationale

### step_size_pos
- Type: `Integer`
- Allowed values:
	- Any integer from 1 to the size of the sequences in the promoter set
- Default value: `101`
- Rationale

### step_size_neg
- Type: `Integer`
- Allowed values:
	- Any integer from 1 to the size of the sequences in the background set
- Default value: `101`
- Rationale