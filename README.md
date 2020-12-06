# Welcome to Crispery
A Python3.7 based program that counts guideRNA occurences in fastq files

Crispery requires absolutly no instalation whatsoever, and can work with any CRISPRi experimental setup.

 There are 4 program versions:

 A windows executable file (Crispery_windows.exe);
 
 A Linux ELF binary file (Crispery_linux);
 
 A MacOS executable (Crispery_macOS); 
 
 As standalone Python code (Crispery.py)

# How to use it

# 1.	
Obtain a .csv file (this format can be obtained using the "save as" option in excel) with the nucleotide sequences of all used guideRNAs, and their respective names (any name can be given, as long as it doesnt repeat). See the provided "D39V_guides.csv" sample file.

| sgRNA0001 | AATAGCATAGAAATCATACA |
|-----------|----------------------|
| sgRNA0002 | AGTGTTGATTTACCAACGTT |


# 2.	
Download the Crispery software version appropriate to the intended operating system, as well as the provided template “inputs.txt” file, into the same folder.

# 3.
Open the “inputs.txt” file.


# Inputs
All the inputs to the program are given via the "inputs.txt" file.

This file MUST be located in the same folder the program is running on.

To run the program, three input absolute paths are required:

A path to the folder with either:

1. all the compressed sequencing files (at the moment, the program unzips .gz files only)

   or

2. all the umcompressed .fastq files

A path to the .csv file with the sgRNAs. See example file for layout (remove any headers).

A path to the output folder (for safety, a subfolder will then always be created on this directory)

The file extension type (default = fastq.gz) (change to the appropriate extension if uncompressed, for example ".fastq") 

If required, several parameters can also be adjusted :

The number of allowed missmatches per sgRNA (default = 1)

The minimal sequencing phred-score for each nucleotide (default = 30)

The start position of the sgRNA within the read (default = 0, meaning the sequenced sgRNA is located at the first position of the read sequence)

The lenght of the sgRNA in bp (default = 20)


# Before running:

Remove any "Undetermined" or unwanted fastq.gz files from the input folder as the program will atempt to align all the fastq files in the input folder.

Overwrite the template "inputs.txt" with the appropriate information.

# Running 


On MSWindows and MacOS:

Double click the program icon. 



On Linux:

Using the terminal, navigate into the folder where Crispery is located and type the following commands:


 i.	chmod +x ./Crispery
 
 ii.	./Crispery -h



Crispery is coded to maximize any computer's processing power (it runs multiprocessed, so it can process various samples simultaneously). It is therefore advisable to not heavilly use the computer while crispery is running to avoid constraining the processor.


Depending on the used computer, crispery might take a few minutes to run. If no errors are shown, crispery is still running. GIVE IT TIME!

If for some reason crispery is closed before completing. please delete the output folder in its entirity as the outputed uncompressed .fastq files might have become corrupted. 


A completion message will be given at the end


# Output

Upon completion, several files should be seen in the indicated output folder: 

a.	The uncompressed “*.fastq” files; 

b. “*_reads.csv” files corresponding to the read counts per guideRNA per inputted sequencing file; 

c.	A “compiled_stats.txt” containing all the relevant input/output information about the Crispery analysis; 

d.	A “compiled.csv” file with the compilation of all the read counts per guideRNA in all the inputted files. Use this latter in the next steps of the data analysis pipeline. 

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
3. Modify the "inputs.txt" to match your system´s apropriate paths and file names.
4. Run crispery

The expected example output file is given: "compiled.csv"



Common Errors:


| Error Message | Probable Cause/Fix |
| ------------- | ------------------ |
| Check the 'inputs.txt' file. Some parameters are missing. | Check if no lines were added or deleted by mistake in the “inputs.txt” file. Check that all the inputs are delimited by “”. |
| No "inputs.txt" file found. Please copy the correct file to the following directory: |The program is not running on the same directory "inputs.txt" is located. Please move the "inputs.txt" file to the directory indicated in the error message. |
| Warning!! X and Y share the same sequence. Only X will be considered valid. |X and Y correspond to guideRNA names. The indicated entries have the same sequence, and only the first will be considered valid. |
| Check the path to the sgRNA file. No file found in the following path: X | Confirm the indicated path (X) to the guideRNA .csv file is correct |
| Check the path to the X files folder. No files of this type found. | Confirm the indicated path to the folder with the sequencing (X) files is correct. |
| Program doesn’t initialize | Confirm the downloaded program is the appropriate one for the current operating system. Contact the Crispery developer if the issue persists. |
| Program crashes, or behaves unexpectedly. | Check if the program behaves as expected with the provided sample data. Contact the Crispery developer if the issue persists. |
