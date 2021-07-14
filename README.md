# Welcome to Crispery
A Python3.7 based program that counts guideRNA occurences in fastq files

Crispery requires absolutly no instalation whatsoever, and can work with any CRISPRi experimental setup.

The program is available as a standalone executable on MSwindows and MacOS, and can be downloaded from Zenodo by acessing the folowing link: http://doi.org/10.5281/zenodo.5079799

crispery is also available as a python package (https://pypi.org/project/crispery/)

`pip install crispery`

In here I share the original Python3 source code, Crispery.py;
The instructions on how to use crispery;
And some test data to run the program.


# Before running:

Remove any "Undetermined" or unwanted fastq.gz files from the input folder as the program will atempt to align all the fastq files in the input folder.


# How to use it

There are two versions of the program, with and without a basic user interface. 


The windows version only uses the user interface. 
For macOS there are both versions, as the graphical interface can be buggy. If one fails, please try the other one.
For Linux there is only the non graphical interface. 


# Using the executable files (graphical mode only):


# 1.	
Obtain a .csv file (this format can be obtained using the "save as" option in excel) with the nucleotide sequences of all used sgRNAs, and their respective names (any name can be given, as long as it doesnt repeat). See the provided "D39V_guides.csv" sample file.

| sgRNA0001 | AATAGCATAGAAATCATACA |
|-----------|----------------------|
| sgRNA0002 | AGTGTTGATTTACCAACGTT |


# 2.	
Download the Crispery software version appropriate to the intended operating system.

# 3. 
Double click the program icon. 

# 4.
The program will initialize and ask, in turn, for directories and file paths. See the "inputs" section below for an explanation of these inputs.



# Using the non executable files (graphical and non-graphical mode)

# 1.	
Obtain the .csv file with the sgRNAs like in step 1 of the graphical interface instructions.

# 2.	
Download the crispery Python3 module using pip install: 
`pip install crispery`.

# 3. For starting the graphical interface mode:
type `python -m crispery`

# 3.1 For starting the non-graphical interface mode:
type `python -m crispery -c --s "c/path/seqfiledir" --g "c/path/sgrna.csv" --o "c/path/outputfolder" --se ".fastq.gz"`

There are also several optional parameters. For their description and input type. A more indepth description is provided below:

`python -m crispery -h`

 `-h, --help  show this help message and exit `

 `-c [C]      cmd line mode`
  
 `--s S       The full path to the directory with the sequencing files`
  
 `--g G       The full path to the .csv file with the sgRNAs.`
  
 `--o O       The full path to the output directory`
  
 `--se SE     Sequencing file extenction (ie:'.fastq.gz')`
  
 `--m M       number of allowed mismatches (default=1)`
  
 `--ph PH     Minimal Phred-score (default=30)`
  
 `--st ST     guideRNA start position in the read (default is 0==1st bp)`
  
 `--l L       guideRNA length`
  
 `--r R       ram saving mode (only appropriate for mismatch searching) `


# Inputs

To run the program, three input paths are required:

# 1 directory containing the sequencing files

A path to the folder with either:

1. all the compressed sequencing files (at the moment, the program unzips .gz files only)

   or

2. all the uncompressed .fastq files

# 2 the path to the sgRNA .csv file

A path to the .csv file with the sgRNAs. See example "D39V_guides.csv" for layout (remove any headers).

# 3 the output directory

A path to the output folder (for safety, a subfolder will always be created on this directory)

# 4 Parameters

The file extension type (default = fastq.gz) (change to the appropriate extension if uncompressed, for example ".fastq") 

The minimal sequencing phred-score for each nucleotide (default = 30)

The start position of the sgRNA within the read (default = 0, meaning the sequenced sgRNA is located at the first position of the read sequence)

The lenght of the sgRNA in bp (default = 20)

The number of allowed missmatches per sgRNA (default = 1)

RAM saving mode (default = no) 
Only useful when allowing missmatch search, as search speed is increased by ~40% due to caching. 
When in RAM saving mode, Crispery should take about 200MB of RAM. 
When NOT in RAM saving mode, several GB might be required. 

# While Running

=================================

Crispery is coded to maximize any computer's processing power (it runs multiprocessed, so it can process various samples simultaneously). It is therefore advisable to not heavilly use the computer while crispery is running to avoid constraining the processor.

When running Crispery in the compiled form, the initializating sequence might take up to a minute. Crispery will be operational when "Version X.X.X" appears on the window.
Depending on the used computer, crispery might take a few minutes to run. If no errors are shown, crispery is still running. GIVE IT TIME!

If for some reason crispery is closed before completing. please delete the output folder in its entirety as the outputed uncompressed .fastq files might have become corrupted. 

=================================

A completion message will be given at the end

+++++++++
Note on mismatch searching:
When performing mismatch searching, especially when large sgRNA libraries and/or large sequencing datasets are used, crispery might take a few hours to run. 
In this case it is advisable to first run crispery without mismatch search (see parameters), and check the output. Non mismatch search uses hashing, and thus it is fast.
+++++++++

# Output

Upon completion, several files should be seen in the indicated output folder: 

a.	The uncompressed “*.fastq” files; 

b. “*_reads.csv” files corresponding to the read counts per guideRNA per inputted sequencing file; 

c.	A “compiled_stats.csv” containing all the relevant input/output information about the Crispery analysis; 

d. A bar plot "reads_plot.png" representing the total number of reads, and valid reads, per sample; 

e.	A “compiled.csv” file with the compilation of all the read counts per guideRNA in all the inputted files. Use this latter in the next steps of the data analysis pipeline. 


# Short Explanation

Crispery will return the read counts for all the guideRNAs present in the input file. 
A read will be aligned to its guideRNAs if the minimum quality score in each nucleotide is >= the indicated phred-score,
and if there is less than the indicated allowed missmatches. 
Like said before, these parameters can be modified by the user.

However, why these parameters?

Base quality filtering using Q>=30 means there is a 0.01% chance of a given nucleotide being miss-sequenced. 
To assure alignment quality, the program filters out by default any reads that have nucleotides with a Q < 30.

Why the missmatch?

To avoid a too highly stringent cutoff.
Allowing a missmatch allows the alignment of reads to their guideRNAs when just a single nucleotide is wrongly sequenced. 
Even at a 0.01% chance (Q>=30 default) this event is bound to happen due to the large sample size.

However, there is a safe mechanism in place to prevent 2 or more guideRNAs with missmatchs from being aligned to the same read (the read is discarded in this case, as there is no way of knowing to which guideRNA the read aligns to)

# Troubleshooting


Running Crispery with example data:

1. Download the "D39V_guides.csv" file
2. Download the "example.fastq.gz"
3. Run crispery

In this example, sgRNA0850 and sgRNA867 share the same sequence; this will appear as a warning message.

The expected example output file is given: "compiled.csv"



Common Errors:


| Error Message | Probable Cause/Fix |
| ------------- | ------------------ |
|Please enter a valid directory for the folowing request: (see printed requests). Please close the program and start again | You have entered an invalid file or folder location. Check that the correct folder/file was selected. When choosing a file click on the file. When selecting a folder, double click (go into the folder, no files will be visible. Its normal).|
|Please confirm that all the input boxes are filled. Some parameters are missing. Press any key to exit | Click OK on the pop up parameter box. The wrong parameter format was entered. Please restart the program and re-introduce ALL the required paths (you must use the btowse buttons for this), and parameters (or leave at default). Please close the program and start again |
| Only numeric values are accepted in the folowing fields: ... | Introduce parameters that respect the default format (text and integers, when appropriate). |
| No "inputs.txt" file found. Please copy the correct file to the following directory: |The program is not running on the same directory "inputs.txt" is located. Please move the "inputs.txt" file to the directory indicated in the error message. |
| Warning!! X and Y share the same sequence. Only X will be considered valid. |X and Y correspond to guideRNA names. The indicated entries have the same sequence, and only the first will be considered valid. |
| Check the path to the sgRNA file. No file found in the following path: X | Confirm the indicated path (X) to the guideRNA .csv file is correct |
| Check the path to the X files folder. No files of this type found. | Confirm the indicated path to the folder with the sequencing (X) files is correct. |
| Program doesn’t initialize | Confirm the downloaded program is the appropriate one for the current operating system. Contact the Crispery developer if the issue persists. |
| Program crashes, or behaves unexpectedly. | Restart the program and carefully indicate the appropriate required input paths. Check if the program behaves as expected with the provided sample data. Contact the Crispery developer if the issue persists. |
