import csv
import glob
import os
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor
import multiprocessing 
from platform import system
from regex import match
from time import time
import easygui as ezi
import matplotlib.pyplot as plt
import numpy as np

#####################

class SgRNA:
    
    """ Each sgRNA will have its own class instance, where the read counts will be kept.
    Each sgRNA class is stored in a dictionary with the sequence as its key. 
    See the "guides loader" function """    
    
    def __init__(self, name, sequence, counts):
        self.name = name
        self.sequence = sequence
        self.counts = counts
        
        
def path_finder_seq(folder_path, extension, separator): 
    
    """ Finds the correct file paths from the indicated directories,
    and parses the file names into an ordered list for later use"""
    
    pathing = []
    for filename in glob.glob(os.path.join(folder_path, extension)):
        stop = filename[::-1].find(separator) + 1 
        pathing.append([filename[-stop:]] + [filename] + [os.path.getsize(filename)]) 

    if extension != '*reads.csv':
        
        """sorting by size makes multiprocessing more efficient 
            as the bigger files will be ran first, thus maximizing processor queing """
        
        ordered = [path[:2] for path in sorted(pathing, key=lambda e: e[-1])][::-1]
     
        if ordered == []:
            input(f"Check the path to the {extension[1:]} files folder. No files of this type found.\n Press any key to exit")

    else:
        ordered = [path[:2] for path in sorted(pathing, reverse = False)]

    return ordered

def unzip(file,write_path):
    
    """ loaded from "unpack" function.
    unzips the .gz files into the directory"""
    
    f = gzip.open(file, 'rb')
    
    if not os.path.isfile(write_path): 
        
        with open(write_path, 'wb') as f_out:
            shutil.copyfileobj(f, f_out)
                
def unpack(ordered,directory):
    
    """ gets the names from the fastq.gz files, and parses their new respective
    names and paths to the "unzip" function for .gz unzipping"""
    
    write_path_save,filing = [],[]
    
    for name, filename in ordered: 
        filing.append(filename)
        write_path_save.append(directory + name[:-len(".gz")])
    
    with ThreadPoolExecutor() as executor: # multihreaded for unzipping files
        executor.map(unzip,filing,write_path_save)

    return write_path_save

def guides_loader(guides):
    
    """ parses the sgRNA names and sequences from the indicated sgRNA .csv file.
    Creates a dictionary using the sgRNA sequence as key, with an instace of the 
    SgRNA class as respective value. If duplicated sgRNA sequences exist, 
    this will be caught in here"""
    
    print("\nAligning sgRNAs")
    
    if not os.path.isfile(guides):
        input("\nCheck the path to the sgRNA file.\nNo file found in the following path: {}\nPress any key to exit".format(guides))
    
    sgrna = {}
    
    with open(guides) as current: 
        for line in current:
            line = line[:-1].split(",")
            sequence = line[1].upper()
            sequence = sequence.replace(" ", "")
            
            if sequence not in sgrna:
                sgrna[sequence] = SgRNA(line[0],sequence, 0)
                
            else:
                print("\nWarning!!\n{} and {} share the same sequence. Only {} will be considered valid.\n".format(sgrna[sequence].name, line[0],sgrna[sequence].name))

    return sgrna

def reads_counter(raw, quality_set, start, lenght, sgrna, mismatch):
    
    """ Reads the fastq file on the fly to avoid RAM issues. 
    Each read is assumed to be composed of 4 lines, with the sequense being 
    on line 2, and the basepair quality on line 4. 
    Every read is trimmed based on the indicated sgRNA positioning. 
    The quality of the obtained trimmed read is crossed against the indicated
    Phred score for quality control.
    If the read has a perfect match with a sgRNA, the respective sgRNA gets a 
    read increase of 1 (done by calling .counts from the respective sgRNA class 
    from the sgRNA class dictionary).
    If the read doesnt have a perfect match, it is sent for mismatch comparison
    via the "imperfect_alignment" function.
    """
    
    n = set("N")
    reading = []
    perfect_counter, imperfect_counter, reads = 0,0,0
    guide_len = start + lenght
    
    with open(raw) as current:
        for line in current:
            reading.append(line[:-1])
            
            if len(reading) == 4: #a read always has 4 lines
                
                seq = reading[1][start:guide_len].upper()
                quality = reading[3][start:guide_len]
                
                reading = []
                reads += 1
                
                if (len(quality_set.intersection(quality)) == 0) & \
                    (len(n.intersection(seq)) == 0):
                        
                    if seq in sgrna:
                        sgrna[seq].counts += 1
                        perfect_counter += 1
                    
                    elif mismatch != 0:
                        sgrna, imperfect_counter = imperfect_alignment(seq, sgrna, mismatch, imperfect_counter)
                    
    return reads, perfect_counter, imperfect_counter, sgrna
  

def imperfect_find(seq,guide,diffnumber):

     """ regex function to compare if the inputed sequences are similar to
     the indicated mismatch degree. returns True if they are similar enough"""
    
    if match("(%s" % seq + "){s<=%s" % diffnumber + "}", guide):
        return True
            

def imperfect_alignment(sequence, sgrna, mismatch, counter):
    
    """ for the inputed read sequence, this compares if there is a sgRNA 
    with a sequence that is similar to it, to the indicated mismatch degree
    (see "imperfect_find" function).
    if the read can be atributed to more than 1 sgRNA, the read is discarded.
    If all conditions are meet, the read goes into the respective sgRNA count 
    score"""
    
    found = 0
    for guide in sgrna:
        
        finder = imperfect_find(sequence, guide, mismatch)

        if finder:

            found += 1
            found_guide = guide
            
            if found >=2:
                break           

    if found == 1:
        sgrna[found_guide].counts += 1
        counter += 1

    return sgrna, counter


def aligner(raw, guides, out, quality_set,mismatch,i,o,sgrna,version,separator, start, lenght):

    """ Runs the main read to sgRNA associating function "reads_counter".
    Creates some visual prompts to alert the user that the samples are being
    processed. Some on the fly quality control is possible (such as making sure 
    the total number of samples is correct, getting an estimate of the total
    number of reads per sample, and checking total running time"""
    
    tempo = time()
       
    print(f"Processing file {i+1} out of {o}")

    reads, perfect_counter, imperfect_counter, sgrna = reads_counter(raw, quality_set, start, lenght, sgrna, mismatch)

    master_list = [["#sgRNA"] + ["Reads"]]
    for guide in sgrna:
        master_list.append([sgrna[guide].name] + [sgrna[guide].counts])

    tempo = time() - tempo
    if tempo > 60:
        timing = str(round(tempo / 60, 2)) + " minutes"   
    else:
        timing = str(round(tempo, 2)) + " seconds"

    name = raw[-raw[::-1].find(separator):-len(".fastq")]
    stats_condition = f"#script ran in {timing} for file {name}. {perfect_counter+imperfect_counter} reads out of {reads} were considered valid. {perfect_counter} were perfectly aligned. {imperfect_counter} were aligned with mismatch"
    
    master_list.sort(key = lambda master_list: master_list[0]) #alphabetical sorting
    master_list.insert(0,[stats_condition])
    csvfile = out[:-out[::-1].find(".")-1] + "_reads.csv"
    csv_writer(csvfile, master_list)
    
    print(stats_condition[1:]) # quality control


def csv_writer(path, outfile):
    
    """ writes the indicated outfile into an .csv file in the directory"""
        
    with open(path, "w", newline='') as output: #writes the output
        writer = csv.writer(output)
        writer.writerows(outfile)

def inputs_asker(msgs, input_parameters, directr):
    
    """ Opens the user interface file browsing dialog boxes"""
    
    print(msgs)
    ezi.msgbox(msg=msgs)
    
    if directr:
        input_parameters.append(ezi.diropenbox(msg = msgs))
    else:
        input_parameters.append(ezi.fileopenbox(msg=msgs))
        
    return input_parameters

def input_getter_graphical():

    """ Handles the parsing of the user interface parameters"""
    
    input_parameters = []
    msgs = ["Select the directory containing the sequencing files",
            "Select the sgRNA .csv file",
            "Select the output directory",
            "Running Parameters"]
    
    ### asking for directory with the fastq.gz files
    input_parameters = inputs_asker(msgs[0], input_parameters, True)
    
    ### asking for directory with the sgRNA files
    input_parameters = inputs_asker(msgs[1], input_parameters, False)
    
    ### asking for the output directory
    input_parameters = inputs_asker(msgs[2], input_parameters, True)
    
    ### updating the parameters
    title = "Enter the appropriate running parameters (default)"
    
    fieldNames = ["Sequencing file extension",\
                  "Allowed mismatches",\
                "Minimal sgRNA Phred-score",\
                "sgRNA start position in the read",\
                "sgRNA length"]
        
    default = [".fastq.gz",1,30,0,20]  
    input_parameters.append(ezi.multenterbox(msgs[3],title, fieldNames, default))
    
    ## error messages
    for entry, msg in zip(input_parameters, msgs):
        if (entry == None) & (msg != msgs[3]):
            input(f"Please enter a valid directory for the folowing request: '{msg}'\nPlease close the program and start again")
            break
        
        elif (entry == None) & (msg == msgs[3]):
            input("Please do not cancel the parameters box. Click OK next time.\nPlease close the program and start again")
            break
        
        elif (entry != None) & (msg == msgs[3]) & (len(input_parameters [-1]) != 5):
            input("The wrong parameter format was entered\nPlease restart the program and re-introduce the parameters (or leave at default)\nPlease close the program and start again")
            break
        
    return input_parameters

def input_getter_txt():

    """ Handles the parsing of the .txt only parameters"""
    
    inputs = os.path.join(os.getcwd(), "inputs.txt")
    
    if not os.path.isfile(inputs): 
        input(f'\nNo "inputs.txt" file found. Please copy the correct file to the following directory: {inputs}')
    
    input_parameters = []
    with open(inputs) as current:
        for line in current:
            if ("#" not in line[0]) and ("\n" not in line[0]):
                line = line[:-1].split('"')
                input_parameters.append(line[1])
                
    if len(input_parameters) != 8:
        input("\nCheck the 'inputs.txt' file. Some parameters are missing.\n")
                
    return input_parameters

def initializer():
    
    """ Handles the program initialization process.
    Makes sure the path separators, and the input parser function is correct
    for the used OS.
    Creates the output diretory and handles some parameter parsing"""
    
    version = "1.4.1"
    
    print("\nVersion: {}".format(version))
    
    if system() == 'Windows':
        separator = "\\"
        folder_path, guides, out, (extension, mismatch, phred, start, lenght) = input_getter_graphical()
        
    else:
        separator = "/"
        
        if system() == 'Darwin':
            folder_path, guides, out, (extension, mismatch, phred, start, lenght) = input_getter_graphical()
            
        else:
            folder_path, guides, out, extension, mismatch, phred, start, lenght = input_getter_txt()
    
    extension = f'*{extension}'
    
    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI" #Phred score

    quality_set = set(quality_list[:int(phred)-1])
    
    directory = os.path.join(out, "unpacked")
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    print("\nRunning with parameters:\n{} mismatch allowed\nMinimal Phred Score per bp >= {}\n".format(mismatch, phred))
    print("All data will be saved into {}".format(directory))


    return folder_path, guides, int(mismatch), quality_set, directory, \
        version, int(phred), separator, int(start), int(lenght), extension 


def compiling(directory,phred,mismatch,version,separator):
    
    """ Combines all the individual processed .csv files into one final file.
    Gathers the individual sample statistic and parses it into "run_stats" """
    
    ordered_csv = path_finder_seq(directory, '*reads.csv',separator)
    
    headers = [f"#Crispery version: {version}"] + \
            [f"#Mismatch: {mismatch}"] + \
            [f"#Phred Score: {phred}"]
            
    compiled = {} #dictionary with all the reads per sgRNA
    head = ["#sgRNA"] #name of the samples
    for name, file in ordered_csv:
        head.append(name[name.find(separator)+len(separator):name.find(".csv")])
        with open(file) as current:
            for line in current:
                line = line[:-1].split(",")
                if "#" not in line[0]:
                    if line[0] in compiled:
                        compiled[line[0]] = compiled[line[0]] + [int(line[1])]
                    else:
                        compiled[line[0]] = [int(line[1])]
                        
                elif "#sgRNA" not in line[0]:
                    headers.append(line[0][1:]+"\n")
  
    path = ordered_csv[0][1]
    out_file = path[:-path[::-1].find(separator)-1]
    
    run_stats(headers,out_file,separator)
                        
    final = []
    for sgrna in compiled:
        final.append([sgrna] + compiled[sgrna])
        
    final.insert(0, head)
    
    csvfile = out_file + separator + "compiled.csv"
    csv_writer(csvfile, final)
        
    input("\nAnalysis successfully completed\nAll the reads have been compiled into the compiled.csv file.\nPress any key to exit")

def run_stats(headers, out_file,separator):
    
    """ Manipulates the statistics from all the samples into one file that can
    be used for downstream user quality control aplications. Creates a simple
    bar graph with the number of reads per sample"""
    
    ### parsing the stats from the read files
    global_stat = [["#Sample name", "Running Time", "Running Time unit", \
                   "Total number of reads in sample", "Total number of reads that passed quality control parameters", \
                    "Number of reads that were aligned without mismatches", "Number of reads that were aligned with mismatches"]]
        
    for run in headers:
        if "script ran" in run:
            parsed = run.split()
            global_stat.append([parsed[7][:-1]] + [parsed[3]] + [parsed[4]] + [parsed[12]] + [parsed[8]] + [parsed[16]] + [parsed[20]])
        else:
            global_stat.insert(0,[run])
    
    csvfile = out_file + separator + "compiled_stats.csv"
    csv_writer(csvfile, global_stat)
    
    ### plotting
    header_ofset = 4
    fig, ax = plt.subplots()
    width = 0.4
    for i, (a,b,c,d,e,f,g) in enumerate(global_stat[header_ofset:]):   

        plt.bar(i+width/2, int(d), width,  capsize=5, color = "darkorange", hatch="//")
        plt.bar(i-width/2, int(e), width, capsize=5, color = "slateblue", hatch="\\\\\\")

    ax.set_xticks(np.arange(len([n[0] for n in global_stat[header_ofset:]])))
    ax.set_xticklabels([n[0] for n in global_stat[header_ofset:]])
    plt.xlabel('Sample')
    plt.ylabel('Number of reads')
    plt.xticks(rotation=45)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(["Total number of reads in sample", "Total number of reads that passed quality control parameters"], \
              loc='upper center', bbox_to_anchor=(0.5, -0.4),ncol=1,prop={'size': 8})
        
    plt.gcf().subplots_adjust(bottom=0.4)
    
    plt.savefig(f"{out_file}{separator}reads_plot.png", dpi=300)


def multi(files,guides, write_path_save, quality_set,mismatch,sgrna,version,separator, start, lenght):
    
    """ starts and handles the parallel processing of all the samples by calling 
    multiple instances of the "aligner" function (one per sample) """
    
    cpu = multiprocessing.cpu_count()
    
    if cpu >= 2:
        cpu -= 1
        
    pool = multiprocessing.Pool(processes = cpu)
    
    for i, (name, out) in enumerate(zip(files, write_path_save)):
        pool.apply_async(aligner, args=(name, guides, out, quality_set,mismatch,i,len(files),sgrna,version,separator, start, lenght))
        
    pool.close()
    pool.join()


def input_file_type(ordered, extension, directory):

    """ funnels the sequencing files to either unzipping, or direct processing,
    depending on the file extension they have"""
    
    if '.gz' in extension:
        
        print("\nUnpacking .gz files")
        files = unpack(ordered,directory)
        write_path_save = files
        
    else:
        
        files,write_path_save = [], []
        for name, filename in ordered:
            files.append(filename)
            write_path_save.append(directory + name)
    
    return files, write_path_save

def main():
    
    """ Runs the program by calling all the appropriate functions"""
    
    ### parses all inputted parameters
    folder_path, guides,mismatch, quality_set,directory, \
    version,phred,separator,start, lenght, extension = initializer()
    
    ### parses the names/paths, and orders the sequencing files
    ordered = path_finder_seq(folder_path, extension, separator)
    
    ### parses the sequencing files depending on whether they require unzipping or not
    files, write_path_save = input_file_type(ordered, extension, directory)
    
    ### loads the sgRNAs from the input .csv file. 
    ### Creates a dictionary "sgrna" of class instances for each sgRNA
    sgrna = guides_loader(guides)
    
    ### Processes all the samples by associating sgRNAs to the reads on the fastq files.
    ### Creates one process per sample, allowing multiple samples to be processed in parallel. 
    multi(files,guides, write_path_save, quality_set,mismatch,sgrna,version,separator, start, lenght)
    
    ### Compiles all the processed samples from multi into one file, and creates the run statistics
    compiling(directory,phred,mismatch,version,separator)

##############
    
if __name__ == "__main__":
    multiprocessing.freeze_support() # required to run multiprocess as .exe on windows
    main()