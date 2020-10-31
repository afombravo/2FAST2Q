# Welcome to Crispy
A Python3.7 based program that counts guideRNA occurences in fastq files

 There are 3 program versions:

 A windows executable file (Crispy.exe);
 
 A Linux ELF binary file (Crispy);
 
 As a standalone Python code (Crispy.py)

# Inputs
All the inputs to the program are given via the "inputs.txt" file.

This file MUST be located in the same folder the program is running on.

To run the program, three input absolute paths are required:
A path to the folder with all the sequencing files (at the moment, the program unzips .gz files only. However, any uncompressed fastq file can be read)

A path to the .csv file with the sgRNAs. See example file for layout.

A path to the output folder (a constant named subfolder will always be created on this directory)

The file extension type (default = fastq.gz) (change to the appropriate extension if uncompressed) 

If required, several parameters can also be adjusted :

The number of allowed missmatches per sgRNA (default = 1)

The minimal sequencing phred-score for each nucleotide (default = 30)

The start position of the sgRNA within the read (default = 0, meaning the sequenced sgRNA is located at the first position of the read sequence)

The lenght of the sgRNA (default = 20)


# Before running:
Remove any "Undetermined" or unwanted fastq.gz files from the input folder as the program will atempt to align all the fastq.gz in the input folder.

Overwrite the template with the appropriate information, and start the program. 
A completion message will be given at the end


# Short Explanation
Crispy will return the read counts for all the guideRNAs present in the input file. 
A read will be aligned to its guideRNAs if the minimum quality score in each nucleotide is >= the indicated phred-score,
and if there is less than the indicated allowed missmatches. 
Like said before, these parameters can be modified by the user.

However, why these parameters?

Base quality filtering using Q>=30 means there is a 0.01% chance of a given nucleotide being "wrongly sequenced". 
To assure alignment quality, the program filters out by default any reads that have nucleotides with a Q < 30.

Why the missmatch?

To avoid a too highly stringent cutoff.
Allowing a missmatch allows the alignment of reads to their guideRNAs when just a single nucleotide is wrongly sequenced. 
Even at a 0.01% chance (Q>=30 default) this event is bound to happen due to the large sample size.

However, there is a safe mechanism in place to prevent 2 or more guideRNAs with missmatchs from being aligned to the same read (the read is discarded in this case, as there is no way of knowing to which guideRNA the read aligns to)

# Output

The program will output a file called compiled.csv with all the compiled information about all the samples.

All conditions also get their own "XX_reads.csv" file.

There is also a compiled_stats.txt file with all the parameter and statistical information.
