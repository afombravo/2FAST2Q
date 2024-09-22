# Welcome to 2FAST2Q
A Python3 program that counts sequence occurrences in FASTQ files. 
2FAST2Q is ideal for CRISPRi-Seq, and for extracting and counting any kind of information from reads in the fastq format, such as barcodes in Bar-seq experiments.
2FAST2Q can work with sequence mismatches, Phred-score, and be used to find and extract unknown sequences delimited by known sequences.  
2FAST2Q can extract multiple features per read using either fixed positions or delimiting search sequences.

2FAST2Q can work with any classic CRISPRi experimental setup, or be used for any kind of sequence extraction from FASTQ files.

2FAST2Q is primarily available as a python package (https://pypi.org/project/fast2q/) and can be installed with the following comand:

`pip install fast2q `

In here I share the original Python3 source code, fast2q.py;
The instructions on how to use 2FAST2Q;
And some test data to run the program.

2FAST2Q is published as part of a CRISPRi-seq protocol:
https://www.nature.com/articles/s41596-021-00639-6

A more in depth description of the program is also available:


Bravo AM, Typas A, Veening J. 2022. 2FAST2Q: a general-purpose sequence search and counting program for FASTQ files. PeerJ 10:e14041 : https://peerj.com/articles/14041/



## How to use it

There are two versions of the program, with and without a basic user interface. 


Basic working principle behind 2FAST2Q:
![](https://github.com/afombravo/2FAST2Q/blob/main/graphical_workings.png)



### 1. Instalation
Open a terminal with python3 installed, and download and install the 2FAST2Q Python3 module using pip install: 

`pip install fast2q`

It is then possible to test if 2FAST2Q was correctly installed by running a test with demo data. 
type `python -m fast2q -c -t`

At the end of a successful test, a message displaying "Test successful. 2FAST2Q is working as intended!" should be visible.

### 2. For starting the graphical interface mode:
type `python -m fast2q`

The program will initialize after a few seconds, poping open the folowing window, and starting when 'OK' is selected. See the "inputs" section below for an explanation on these inputs.

![](https://github.com/afombravo/2FAST2Q/blob/main/C_mode.gif)


The default running mode is in "Counter" mode, however the user might want to run 2FAST2Q in 'Extract and Counter' mode where features are not aligned to a reference, but *de novo* extracted from the file based on indicated search sequences. Consider the folowing example:

![](https://github.com/afombravo/2FAST2Q/blob/main/EC_mode.gif)


### 2.1 For starting the non-graphical interface mode:
type `python -m fast2q -c`

When running without specified parameters, 2FAST2Q will assume the current running directory has all the required files:

* one .csv corresponding to features file (not required in 'Extract and Count' mode)
	
* the .FASTQ files

How it looks when running in the default 'Counter' mode:

![](https://github.com/afombravo/2FAST2Q/blob/main/C_mode_cmd.gif)

How it looks when running in the 'Extract and Counter' mode:

![](https://github.com/afombravo/2FAST2Q/blob/main/Ec_mode_cmd.gif)


There are also several optional parameters. For their description and input type. A more in-depth description is provided below:

	 `python -m fast2q -h`

	 `-h, --help  show this help message and exit `

	 `-c [C]      cmd line mode`
	 
	 `-v [V]      prints the current version`

  	 `-t [T]      Runs 2FAST2Q in test mode with example data. `

	 `--s S       The full path to the directory with the sequencing files OR file`

	 `--g G       The full path to the .csv file with the features.`

	 `--o O       The full path to the output directory`
	 
	 `--fn FN     Specify an output compiled file name (default is called compiled)`
	 
	 `--v V       Adds progress bars (default is enabled)`

	 `--m M       number of allowed mismatches (default=1)`

	 `--ph PH     Minimal Phred-score (default=30)`

	 `--st ST     Feature start position in the read (default is 0==1st bp)`

	 `--l L       Feature length`

	 `--us US     Upstream search sequence `

	 `--ds DS     Downstream search sequence `

	 `--msu MSU   mismatches allowed in the upstream sequence`

  	 `--msd MSD   mismatches allowed in the downstream sequence`

  	 `--qsu QSU   Minimal Phred-score (default=30) in the upstream search sequence`

  	 `--qsd QSD   Minimal Phred-score (default=30) in the downstream search sequence`

	 `--mo MO     Running Mode (default=C) [Counter (C) / Extractor + Counter (EC)] `
	 
	 `--cp CP     Number of cpus to be used (default is max(cpu)-2 for >=3 cpus, -1 for >=2 cpus, 1 if 1 cpu `

	 ` --k K       If enabled, keeps all temporary files (default is disabled) `


## Inputs

To run the program, three input paths are required:

### 1  Directory containing the sequencing files (assumed to be the current directory when using the cmd line version and no inputs are given)

A path to the folder with the sequencing files (it doesn´t matter if in .gz or .fastq.gz format as 2fast2q auto determines the correct one). 2FAST2Q will automatically process all the .fastq files that exist in the indicated folder.

### 2  The path to the feature .csv file (optional) (assumed to be the only .csv file in the current directory when using the cmd line version and no inputs are given)
Only needed when searching the fastq file for known sequences, such as with a CRISPRi-Seq experiment.
A path to the .csv file (this format can be obtained using the "save as" option in excel) with the nucleotide sequences of all used features, and their respective names (any name can be given, as long as it doesn’t repeat). See the provided "D39V_guides.csv" sample file. (Optional, only required when running in Counting mode)	

| sgRNA0001 | AATAGCATAGAAATCATACA |
|-----------|----------------------|
| sgRNA0002 | AGTGTTGATTTACCAACGTT |


#### 2.1

2FATS2Q can be used for finding multiple features per read. When such is desirable, the features must be separated by ":", as illustrated here:

| sgRNA0001.1 | AATAGCATAGAAATCATACA:GATTACA |
|-----------|----------------------|
| sgRNA0001 | AATAGCATAGAAATCATACA |


In this case, sgRNA0001.1 corresponds to a double sequence. Only reads containing BOTH sequences will be aligned to this sgRNA. If only the first sequence of the 2 is found, it will align to sgRNA0001, if only the second sequence is found, it will fail to align anywhere. Only the combinations presente in the .csv file will be considered. 
See section 4 for instructions on how to perform multiple sequence searches per read.

For extracting all possible combinations in a file, one can use the "extract and count" mode (extract all found features without alignments) (`--mo EC`). In this case, no .csv is required as input.



### 3 the output directory

A path to the output folder (for safety, a subfolder will always be created on this directory) (2fast2q automatically creates a subdirectory within the current directory when using the cmd line version and no inputs are given)


### 4 Parameters

For extracting all sequences at a certain position in the read select the extractor + Counter (EC) mode. The default is Counter (C) mode only.

Progress Bar. (Default is enabled)

The minimal sequencing phred-score for each nucleotide (default = 30)

The start position of the feature within the read (default = 0, meaning the sequenced feature is located at the first position of the read sequence)

The length of the feature in bp (default = 20)

The number of allowed mismatches per feature (default = 1). When in extract + Count mode, this parameter is ignored as all different sequences are returned.
2FAST2Q mismatch feature calculates HAMMING distance ONLY

Keep temporary files mode (default = no).
When enabled, deletes all temporary files. To keep all files, change to "n" in the graphical mode, or input the parameter `--k` in the cmd lines.

For extracting all sequences at a certain position in the read select the extractor + Counter (EC) mode. The default is Counter (C) mode only.

If the starting position varies within the read, it is possible to search for a delimiting known sequence, and then extract the sequence before/after it.
In this case, it is allowed to input the following: 

 1) A 5' end search sequence, and the amount of bp the program should inventory after.
 2) A 3' end search sequence, and the amount of bp the program should inventory before.
 3) A 5' and 3' end search sequence, the program will return and count everything in between these two.
 4) How many mismatches are allowed in the search sequence

When searching a read for multiple sequences, one can either do so by:
  1) confirguring different fixed positions by separating all start locations with a ",". For example: "0,20,50" - the program will search for 3 sequences per read, starting at position 0,20, and 50, with the predefided sequence length.
  2) configuring different 5' and 3' search sequences, also separated by "," and inputted as pairs: For example: upstream (`--us`) ATCG,GGTGG & downstream (`--ds`) AATC,GCACAC will initiate, per read, searches for any features between the ATCG * AATC and GGTGG * GCACAC sequences. If found, these 2 sequences will be merged separated by ":" and either try to be aligned against any found features in the .csv file (default), or returned as they are if in "extract and count" mode (`--mo EC`)
 

## While Running

=================================

2FAST2Q is coded to maximize any computer's processing power (it runs multiprocessed, so it can process various samples simultaneously). It is therefore advisable to not heavily use the computer while 2FAST2Q is running to avoid constraining the processor.

2FAST2Q will be operational when its logo appears on the window.
Depending on the used computer, 2FAST2Q might take a few minutes to run, especially with large datasets and when using mismatch finding. If no errors are shown, 2FAST2Q is still running. GIVE IT TIME! 


### macOS use WARNING!


When using the graphical user interface option, it's possible that the interface doesn’t close down after pressing OK and "gets stuck". The program is still running, and progress can be monitored by checking the indicated output folder. When the final "compiled.csv" appears on the folder, the program has finished running and can be closed using any means.
A completion message should be given at the end. In any case, the program will be finish when the compiled.csv file is visible in the directory.


## Output

Upon completion, several files should be seen in the indicated output folder (when running in default mode only b, c (only in cmd line mode), and d will be kept):

	a. 	“*_reads.csv” files corresponding to the read counts per feature per inputted sequencing file; 

	b.	A “compiled_stats.csv” containing all the relevant input/output information about the 2FAST2Q analysis; 

	c.	A “compiled.csv” file with the compilation of all the read counts per feature in all the inputted files. Use this latter in the next steps of the data analysis pipeline. 

	d.	A bar plot "reads_plot.png" with the number of total and quality passed reads (absolute), and in percentage ("reads_plot_percentage.png"), per sample; 

	e. 	2 other violin plots with the distribution of the found features per sample are also presented (normalized for reads per milion, and absolute numbers). The interquartile distribution is also ploted for each sample (25%-75%)


### Short Explanation

2FAST2Q will return the read counts for all the features present in the input file. A read will be aligned to its features if the minimum quality score in each nucleotide is >= the indicated phred-score, and if there is less than the indicated allowed mismatches. Like said before, these parameters can be modified by the user.

However, why these parameters?

Base quality filtering using Q>=30 means there is a 0.01% chance of a given nucleotide being miss-sequenced. To assure alignment quality, the program filters out by default any reads that have nucleotides with a Q < 30.


Why the mismatch?

To avoid a too highly stringent cutoff. Allowing a mismatch allows the alignment of reads to their features when just a single nucleotide is wrongly sequenced. Even at a 0.01% chance (Q>=30 default) this event is bound to happen due to the large sample size.

However, there is a safe mechanism in place to prevent 2 or more features with mismatches from being aligned to the same read (the read is discarded in this case, as there is no way of knowing to which feature the read aligns to)

2FAST2Q mismatch feature calculates HAMMING distance ONLY

## Troubleshooting

Running 2FAST2Q with example data :

Download the "D39V_guides.csv" file
Download the "example.fastq.gz"
Run 2FAST2Q


In this example, sgRNA0850 and sgRNA867 share the same sequence; this will appear as a warning message.


=======


The expected example output file is given: "compiled.csv"


# License

This works is distributed under a GNU General Public License v3.0
