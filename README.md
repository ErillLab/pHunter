# pHunt: Searching for Promoters in Bacterial Genomes

Based on the [PhiSite Promoter Hunter](http://www.phisite.org/main/index.php?nav=tools&nav_sel=hunter)
that use both *PSSMs* and *Gibbs Free Energy* to find putative promoter sites.  

*PSSMs* are used to find the -35 and -10 sequences where the RNA polymerase holoenzyme
binds in order to translate a gene. *Gibbs Free Energy* is used to find sections in the
genome that are thermodynamically predisposed to unzipping (sites where the holoenzyme likely binds to in order to initiate transcription).

## About

The general process by which the program finds promoters is as follows:

1. Generates a motif object for the -35 and -10 sequences using the `Biopython` Library

2. Scans a sequence inputted by the user for promoter hunting

3. Computes -35 and -10 PSSM scores across the inputted sequence, applying strict
spacer limits between potential -35 and -10 sequences.

    * Spacer = # of Base Pairs between the putative -35 and -10 sequence.

4. Proceeds with potential promoters with -35 and -10 scores that are above the left and right thresholds.

    - Only hits that pass the threshold barrier will be reported on the output file

5. Computes the *Gibbs Free Energy* of the area around the -10 sequence

6. Calculates a Final Score by summing the contributions of the *PSSMs* and *Gibbs Free Energy*

7. Generates a file of hits for an input sequence, sorted by Final Score

The original PhiSite Promoter Hunter sums the contributions of the PSSM scores and Gibbs Free Energy
by normalizing each value and applying arbitrary coefficients. Although the program still has functionality to permit the user to compute Final Score as the sum of normalized values, computing the Final Score as the summation of Log-likelihood Ratios is suggested. The program provides this functionality by converting Gibbs Free Energy into a Log-Likelihood Ratio using Free Energy distributions. This method of computing final scores is non-arbitrary and as shown in `benchmarking.md`, superior in performance.

To read more about the steps above, please take a look at `user_manual.md` and `settings.md` in the documentation folder




## Working with the Program

### Python Dependencies
__pHunt__ runs on python 3.9.12 and depends on packages listed in `pHunt_env.yaml`. All dependencies can be installed using pip:

>pip install -r pHunt_env.yaml

### Dependencies
* `numpy`
* `Biopython`

### Running the program
>import pHunt
>
>pHunt.go()

### Input

pHunt expects a local directory set up the same as this Github directory. This is essential to the program's functioning
since the program expects to find certain files in certain subdirectories.

* **data** folder: Contains the files to be used and accessed by the program. Site where the `output file` will be created

  * Must be in the **data** folder:

    * left motif file
    * right motif file
    * Input sequences file

  * Must be in the **data** folder if `mode_fscr` is set to LLR:

    * promoter_set file
    * background_set file

  * Must be in the **data** folder if `use_GibbsFE` is set to True:

    * GibbsFE_file file

  * For more information about these files, please check out `settings.md`

* **src** folder: Contains all the python code (as well as the json file) necessary for the program's functionality

  * Must be in the **src** folder:

    * pHunt python file
    * GibbsFE python file
    * settings.json

The pHunt program uses a JSON file format (`settings.json`) to collect user input. Below is a sample json file


    "mode_fscr": "norm",
    "LLR_specific_parameters": {
    	"promoter_set": "PromEC_seqs_filtered.fas",
    	"background_set": "negativeset.csv",
    	"step_size_genome": 101,
    	"GibbsFE_n_bins": 30,
    	"GibbsFE_pseudocounts": 1
    },
    "use_GibbsFE": false,
    "GibbsFE_related_parameters": {
    	"GibbsFE_file": "GibbsFE.csv",
    	"GibbsFE_windowsize": 50,
    	"lerg": 12,
    	"rerg": 13
    },
    "input_sequences": "PromEC_extend_hundredbp.fasta",
    "use_GCcont_background": false,
    "left_motif": "Eco-35.fas",
    "left_motif_pseudocounts": 0.25,
    "right_motif": "Eco-10.fas",
    "right_motif_pseudocounts": 0.25,
    "motif_threshold": "original",
    "min_spacer_value": 15,
    "max_spacer_value": 18,
    "output_information": "Pseudotestgenomefitted2_Diagnostic.csv",
    "output_type":"all hits"



To use parameters other than the recommended ones (default at settings.json), edit the settings.json file or create a new json file with the same format. Specifics about what each parameter does and permitted values can be found in `settings.md`

## output
pHunt produces a csv file in the data folder with hits sorted by Final Score for each sequence inputted by the user. The file can either contain *all hits* above the threshold for each sequence or only the *top hits* for each sequence (adjustable in the json file).

## Accessing additional documentation about the program

In the documentation subfolder you will find:

* `user_manual.md`: A file explaining what pHunt does and why it does it
* `benchmarking.md`: A file substantiating the default parameters and exploring the performance differences between certain techniques
* `settings.md`: A file explaining each of the json's parameters
