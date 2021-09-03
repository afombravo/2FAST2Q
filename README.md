# Welcome to 2FAST2Q
A Python3 program that counts sequence occurrences in raw FASTQ files. 
2FAST2Q is ideal for CRISPRi-Seq, and for extracting and counting any kind of information from Illumina reads, such as barcodes.
2FAST2Q can work with sequence mismatches, and can be used to find sequences delimited by known sequences in unaligned reads.  

2FAST2Q requires absolutely no installation whatsoever, and can work with any classic CRISPRi experimental setup, or be used for any kind of sequence extraction from FASTQ files.

The program is available as a standalone executable on MSwindows and MacOS, and can be downloaded from Zenodo by accessing the following link (keep in mind a newer version might exist, if so, use it): 
https://zenodo.org/record/5410822

2FAST2Q is also available as a python package (https://pypi.org/project/fast2q/)

`pip install fast2q `

In here I share the original Python3 source code, fast2q.py;
The instructions on how to use 2FAST2Q;
And some test data to run the program.


# Before running:

Remove any "Undetermined" or unwanted fastq.gz files from the input folder as the program will attempt to align all the fastq files in the input folder.


# How to use it

There are two versions of the program, with and without a basic user interface. 


There is a graphical interface version for Windows, MacOS, and Linux.
If for some reason the compiled version fails, please use the souce Python code from PyPI (`pip install fast2q`).


# Using the executable files:


# 1.	
Obtain a .csv file (this format can be obtained using the "save as" option in excel) with the nucleotide sequences of all used features, and their respective names (any name can be given, as long as it doesn’t repeat). See the provided "D39V_guides.csv" sample file. (Optional, only required when running in Counting mode)	

| sgRNA0001 | AATAGCATAGAAATCATACA |
|-----------|----------------------|
| sgRNA0002 | AGTGTTGATTTACCAACGTT |


# 2.	
Download the 2FAST2Q software version appropriate to the intended operating system.

# 3. 
Double click the program icon. 

# 4.
The program will initialize and ask, in turn, for directories and file paths. See the "inputs" section below for an explanation of these inputs.



# Using the non executable files (recommendable if the executable file is buggy)

# 1.	
Obtain the .csv file with the features like in step 1 of the previous instructions. (Optional, only required when running in Counting mode)	

# 2.	
Download the 2FAST2Q Python3 module using pip install: 
`pip install fast2q`.

# 3. For starting the graphical interface mode:
type `python -m fast2q`

# 3.1 For starting the non-graphical interface mode:
type `python -m fast2q -c --s "c/path/seqfiledir" --g "c/path/sgrna.csv" --o "c/path/outputfolder" --se ".fastq.gz"`

There are also several optional parameters. For their description and input type. A more in-depth description is provided below:

`python -m fast2q -h`

 `-h, --help  show this help message and exit `

 `-c [C]      cmd line mode`
  
 `--s S       The full path to the directory with the sequencing files`
  
 `--g G       The full path to the .csv file with the features.`
  
 `--o O       The full path to the output directory`
  
 `--se SE     Sequencing file extenction (ie:'.fastq.gz')`
  
 `--m M       number of allowed mismatches (default=1)`
  
 `--ph PH     Minimal Phred-score (default=30)`
  
 `--st ST     feature start position in the read (default is 0==1st bp)`
  
 `--l L       feature length`
  
 `--r R       ram saving mode (only appropriate for mismatch searching) `
 
 `--us US     Upstream search sequence`
 
 `--ds DS     Downstream search sequence`
 
 `--ms MS     mismatches allowed when searching reads with Up/Down stream sequences`
 
 `--mo MO     Running Mode (default=C) [Counter (C) / Extractor + Counter(EC)]`
 
 `--k K       If enabled, keeps all temporary files (default is enabled)`


# Inputs

To run the program, three input paths are required:

# 1  Directory containing the sequencing files

A path to the folder with either:

1. all the compressed sequencing files (at the moment, the program unzips .gz files only)

   or

2. all the uncompressed .fastq files

# 2  The path to the feature .csv file 
(only needed when searching the fastq file for known sequences, such as with a CRISPRi-Seq experiment)

A path to the .csv file with the features. See example "D39V_guides.csv" for layout (remove any headers).

# 3 the output directory

A path to the output folder (for safety, a subfolder will always be created on this directory)

# 4 Parameters

The file extension type (default = fastq.gz) (change to the appropriate extension if uncompressed, for example ".fastq") 

The minimal sequencing phred-score for each nucleotide (default = 30)

The start position of the feature within the read (default = 0, meaning the sequenced feature is located at the first position of the read sequence)

The length of the feature in bp (default = 20)

The number of allowed mismatches per feature (default = 1)

RAM saving mode (default = no) 
Only useful when allowing mismatch search, as search speed is increased by ~40% due to caching. 
When in RAM saving mode, 2FAST2Q should only take a few MB of RAM. 
When NOT in RAM saving mode, several GB might be required.

Keep temporary files mode (default = yes).
When enabled, deletes all temporary files. To keep all files, change to "n" in the graphical mode, or input the parameter `--k` in the cmd lines.

For extracting all sequences at a certain position in the read select the extractor + Counter (EC) mode. The default is Counter (C) mode only.

If the starting position varies within the read, it is possible to search for a delimiting known sequence, and then extract the sequence before/after it.
In this case, it is allowed to input the following: 
 1) A 5' end search sequence, and the amount of bp the program should inventory after.
 2) A 3' end search sequence, and the amount of bp the program should inventory after.
 3) A 5'and 3' end search sequence, the program will return and count everything in between these two.
 4) How many mismatches are allowed in the search sequence
 
# While Running

=================================

2FAST2Q is coded to maximize any computer's processing power (it runs multiprocessed, so it can process various samples simultaneously). It is therefore advisable to not heavily use the computer while 2FAST2Q is running to avoid constraining the processor.

When running 2FAST2Q in the executable form, the initialization sequence might take up to a minute. 2FAST2Q will be operational when "Version X.X.X" appears on the window.
Depending on the used computer, 2FAST2Q might take a few minutes to run, especially with large datasets and when using mismatch finding. If no errors are shown, 2FAST2Q is still running. GIVE IT TIME! 

\\\\

macOS use WARNING!


When using the graphical user interface option, it's possible that the interface doesn’t close down after pressing OK and "gets stuck". The program is still running, and progress can be monitored by checking the indicated output folder. When the final "compiled.csv" appears on the folder, the program has finished running and can be closed using any means.
A completion message should be given at the end. In any case, the program will be finish when the compiled.csv file is visible in the directory.

\\\\

+++++++++

Note on mismatch searching: When performing mismatch searching, especially with feature libraries with thousands of features and/or when large sequencing datasets are used, 2FAST2Q might take 1-2 hours to run. In this case it is advisable to first run 2FAST2Q without mismatch search (see parameters), and check the output. Non mismatch search uses hashing, and thus it is fast. +++++++++

Output

Upon completion, several files should be seen in the indicated output folder (when running in default mode only c, d, and e will be kept):

a. The uncompressed “*.fastq” files;

b. “*_reads.csv” files corresponding to the read counts per feature per inputted sequencing file;

c. A “compiled_stats.csv” containing all the relevant input/output information about the 2FAST2Q analysis;

d. A bar plot "reads_plot.png" representing the total number of reads, and valid reads, per sample;

e. A “compiled.csv” file with the compilation of all the read counts per feature in all the inputted files. Use this latter in the next steps of the data analysis pipeline.


Short Explanation

2FAST2Q will return the read counts for all the features present in the input file. A read will be aligned to its features if the minimum quality score in each nucleotide is >= the indicated phred-score, and if there is less than the indicated allowed mismatches. Like said before, these parameters can be modified by the user.

However, why these parameters?

Base quality filtering using Q>=30 means there is a 0.01% chance of a given nucleotide being miss-sequenced. To assure alignment quality, the program filters out by default any reads that have nucleotides with a Q < 30.


Why the mismatch?

To avoid a too highly stringent cutoff. Allowing a mismatch allows the alignment of reads to their features when just a single nucleotide is wrongly sequenced. Even at a 0.01% chance (Q>=30 default) this event is bound to happen due to the large sample size.

However, there is a safe mechanism in place to prevent 2 or more features with mismatches from being aligned to the same read (the read is discarded in this case, as there is no way of knowing to which feature the read aligns to)


Troubleshooting

Running 2FAST2Q with example data :

Download the "D39V_guides.csv" file
Download the "example.fastq.gz"
Run 2FAST2Q

In this example, sgRNA0850 and sgRNA867 share the same sequence; this will appear as a warning message.


The expected example output file is given: "compiled.csv"
