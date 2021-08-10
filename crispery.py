import csv
import glob
import os
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor
import multiprocessing 
from platform import system
from time import time
import matplotlib.pyplot as plt
import numpy as np
from numba import njit
import psutil
import argparse
import datetime
#also needs tkinter (imported inside inputs_initializer())

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
            raise Exception

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
        raise Exception
    
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

def reads_counter(raw, sgrna, param):
    
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
    
    def seq2bin(sequence):
        
        """ Converts a string to binary, and then to 
        a numpy array in int8 format"""
        
        byte_list = bytearray(sequence,'utf8')
        return np.array((byte_list), dtype=np.int8)
    
    def binary_converter(sgrna):

        """ Parses the input sgRNAs into binary dictionaries. Converts all DNA
        sequences to their respective binary array forms. This gives some computing
        speed advantages with mismatches."""
        
        from numba import types
        from numba.typed import Dict
        
        container = Dict.empty(key_type=types.unicode_type,
                               value_type=types.int8[:])

        for sequence in sgrna:
            container[sequence] = seq2bin(sequence)
        return container
    
    def unfixed_starting_place_parser(read,upstream,downstream,mismatch):
        
        """ Determines the starting place of a read trimming based on the
        inputed parameters for upstream/downstream sequence matching"""
        
        read_bin=seq2bin(read)
        start,end = None,None
        
        if (upstream is not None) & (downstream is not None):
            start=border_finder(upstream,read_bin,mismatch)
            end=border_finder(downstream,read_bin,mismatch)
            if (start is not None) & (end is not None):
                start+=len(upstream)

        elif (upstream is not None) & (downstream is None):
            start=border_finder(upstream,read_bin,mismatch)
            if start is not None:
                start+=len(upstream)
                end = start + param['length']
            
        elif (upstream is None) & (downstream is not None):
            end=border_finder(downstream,read_bin,mismatch)
            if end is not None:
                start = end-param['length']

        return start,end
            
    if param['miss'] != 0:
        binary_sgrna = binary_converter(sgrna)
    
    n = set("N")
    quality_set = param['quality_set']
    fixed_start = True
    failed_reads = set()
    reading = []
    perfect_counter, imperfect_counter, reads = 0,0,0
    
    # determining the read trimming starting/ending place
    if (param['upstream'] is None) & (param['downstream'] is None):
        end = param['start'] + param['length']
        start = param['start']
    else:
        fixed_start = False
        if param['upstream'] is not None:
            param['upstream']=seq2bin(param['upstream'].upper())
        if param['downstream'] is not None:
            param['downstream']=seq2bin(param['downstream'].upper())
    
    with open(raw) as current:
        for line in current:
            reading.append(line[:-1])
            
            if len(reading) == 4: #a read always has 4 lines
                
                if not fixed_start:
                    start,end=unfixed_starting_place_parser(reading[1],param['upstream'],param['downstream'],param['miss_search'])
                
                if (fixed_start) or ((start is not None) & (end is not None)):
                    seq = reading[1][start:end].upper()
                    quality = reading[3][start:end]

                    if (len(quality_set.intersection(quality)) == 0) & \
                        (len(n.intersection(seq)) == 0):
                        
                        if param['Running Mode']=='C':
                            if seq in sgrna:
                                sgrna[seq].counts += 1
                                perfect_counter += 1
                            
                            elif param['miss'] != 0:
                                read=seq2bin(seq)
                                
                                if not param['ram']: #keeps track of reads that are already known to not align with confidence
                                    if seq not in failed_reads:
                                        sgrna,imperfect_counter,failed_reads = imperfect_alignment(read,seq,binary_sgrna,param['miss'],imperfect_counter,sgrna,failed_reads,param['ram'])
                                else:
                                    sgrna,imperfect_counter,failed_reads = imperfect_alignment(read,seq,binary_sgrna,param['miss'],imperfect_counter,sgrna,failed_reads,param['ram'])
                        else:
                            if seq not in sgrna:
                                sgrna[seq] = SgRNA(seq,seq, 0)
                            else:
                                sgrna[seq].counts += 1
                                perfect_counter += 1
                reading = []
                reads += 1
                
    # clears the cache RAM from each individual file/process
    failed_reads = set()
    
    return reads, perfect_counter, imperfect_counter, sgrna

@njit
def border_finder(seq,read,mismatch): 
    
    """ Matches 2 sequences (after converting to int8 format)
    based on the allowed mismatches. Used for sequencing searching
    a start/end place in a read"""
    
    s=seq.size
    r=read.size
    for i,bp in enumerate(read):
        comparison = read[i:s+i]
        finder = binary_subtract(seq,comparison,mismatch)
        if finder != 0:
            return i
        if i > r:
            break

@njit
def binary_subtract(array1,array2,mismatch):
    
    """ Used for matching 2 sequences based on the allowed mismatches.
    Requires the sequences to be in numerical form"""
    
    miss=0
    for arr1,arr2 in zip(array1,array2):
        if arr1-arr2 != 0:
            miss += 1
        if miss>mismatch:
            return 0
    return 1

@njit
def sgrna_all_vs_all(binary_sgrna,read,mismatch):
    
    """ Runs the loop of the read vs all sgRNA comparison.
    Sends individually the sgRNAs for comparison.
    Returns the final mismatch score"""
    
    found = 0
    for guide in binary_sgrna:
        if binary_subtract(binary_sgrna[guide],read,mismatch):
            found+=1
            found_guide = guide
            if found>=2:
                return
    if found==1:
        return found_guide
    return

def imperfect_alignment(read,seq,binary_sgrna, mismatch, counter, sgrna,failed_reads,ram):
    
    """ for the inputed read sequence, this compares if there is a sgRNA 
    with a sequence that is similar to it, to the indicated mismatch degree
    if the read can be atributed to more than 1 sgRNA, the read is discarded.
    If all conditions are meet, the read goes into the respective sgRNA count 
    score"""
    
    finder = sgrna_all_vs_all(binary_sgrna, read, mismatch)

    if finder is not None:
        sgrna[finder].counts += 1
        counter += 1
    elif not ram:
        failed_reads.add(seq)
        
    return sgrna, counter, failed_reads

def aligner(raw,out,i,o,sgrna,param):

    """ Runs the main read to sgRNA associating function "reads_counter".
    Creates some visual prompts to alert the user that the samples are being
    processed. Some on the fly quality control is possible (such as making sure 
    the total number of samples is correct, getting an estimate of the total
    number of reads per sample, and checking total running time"""
    
    ram_lock()
    tempo = time()
       
    print(f"Processing file {i+1} out of {o}")

    reads, perfect_counter, imperfect_counter, sgrna = reads_counter(raw,sgrna,param)

    master_list = [["#sgRNA"] + ["Reads"]]
    for guide in sgrna:
        master_list.append([sgrna[guide].name] + [sgrna[guide].counts])

    tempo = time() - tempo
    if tempo > 60:
        timing = str(round(tempo / 60, 2)) + " minutes"   
    else:
        timing = str(round(tempo, 2)) + " seconds"

    name = raw[-raw[::-1].find(param['separator']):-len(".fastq")]
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

def inputs_handler(separator):
    
    """ assertains the correct parsing of the input parameters"""
    
    parameters=inputs_initializer(separator)

    try:
        parameters["start"]=int(parameters["start"])
        parameters["length"]=int(parameters["length"])
        parameters["miss"]=int(parameters["miss"])
        parameters["phred"]=int(parameters["phred"])
        parameters["miss_search"]=int(parameters["miss_search"])
    except Exception:
        input("\nOnly numeric values are accepted in the folowing fields:\nsgRNA read starting place;\nsgRNA length;\nmismatch;\nPhred score;\nmismatches in the search sequence.\n\nPlease try again. Press any key to exit")
        raise Exception    
    
    # parsing the RAM saving choice as a bolean
    if parameters['ram'] == "n":
        parameters['ram'] = False
    else:
        parameters['ram'] = True
        
    if parameters['delete'] == "y":
        parameters['delete'] = True
    else:
        parameters['delete'] = False
        
    if parameters['upstream'] == "None":
        parameters['upstream'] = None
        
    if parameters['downstream'] == "None":
        parameters['downstream'] = None
        
    if "Extractor" in parameters['Running Mode']:
        parameters['Running Mode']="EC"
    else:
        parameters['Running Mode']="C"
        
    if parameters['Running Mode']=='C':
        if len(parameters) != 14:
            input("Please confirm that all the input boxes are filled. Some parameters are missing.\nPress any key to exit")
            raise Exception

    return parameters

def inputs_initializer(separator):
    
    """ Handles the graphical interface, and all the parameter inputs"""
    
    from tkinter import Entry,LabelFrame,Button,Label,Tk,filedialog,StringVar,OptionMenu
    
    def restart():
        root.destroy()
        inputs_initializer()
        
    def submit():
        for arg in temporary:
            parameters[arg] = temporary[arg].get()
        root.destroy()
    
    def directory(column,row,parameter,frame):
        filename = filedialog.askdirectory(initialdir =  separator, title = "Select a folder")
        filing_parser(column,row,filename,parameter,frame)
        
    def file(column,row,parameter,frame):
        filename = filedialog.askopenfilename(initialdir = separator, title = "Select a file", filetypes = \
            (("CSV files","*.csv"),("all files","*.*")) )
        filing_parser(column,row,filename,parameter,frame)
    
    def filing_parser(column,row,filename,parameter,frame):
        place = Entry(frame,borderwidth=5,width=50)
        place.grid(row=row,column=column+1)
        place.insert(0, filename)
        parameters[parameter] = filename

    def browsing(keyword,inputs):
        title1,title2,row,column,function=inputs
        frame=LabelFrame(root,text=title1,padx=5,pady=5)
        frame.grid(row=row,column=column, columnspan=2)
        button = Button(frame, text = title2,command = lambda: function(column,row,keyword,frame))
        button.grid(column = column, row = row)
        button_place = Entry(frame,borderwidth=5,width=50)
        button_place.grid(row=row,column=column+1)
        button_place.insert(0, "Use the Browse button to navigate")
        
    def write_menu(keyword,inputs):
        title,row,column,default=inputs
        start = Entry(root,width=25,borderwidth=5)
        start.grid(row=row,column=column+1,padx=20,pady=5)
        start.insert(0, default)
        placeholder(row,column,title,10,1)
        temporary[keyword]=start
        
    def dropdown(keyword,inputs):
        default,option,row,column = inputs
        file_ext = StringVar()
        file_ext.set(default)
        placeholder(row,column,keyword,20,1)
        drop = OptionMenu(root, file_ext, default, option)
        drop.grid(column = column+1,row = row)
        temporary[keyword]=file_ext
        
    def button_click(row, column, title, function):
        button_ok = Button(root,text=title,padx=12,pady=5, width=15,command=function)
        button_ok.grid(row=row, column=column,columnspan=1)
        
    def placeholder(row, column,title,padx,pady):
        placeholder = Label(root, text=title)
        placeholder.grid(row=row,column=column,padx=padx,pady=pady)
        return placeholder
        
    root = Tk()
    root.title("Crispery Input Parameters Window")
    root.minsize(425, 660)
    parameters,temporary = {},{}  

    browsing_inputs = {"seq_files":["Path to the .fastq(.gz) files folder","Browse",1,0,directory],
                       "sgrna":["Path to the sgRNA .csv file","Browse",2,0,file],
                       "out":["Path to the output folder","Browse",3,0,directory]}

    default_inputs = {"start":["sgRNA start position in the read",6,0,0],
                      "length":["sgRNA length",7,0,20],
                      "miss":["Allowed mismatches",8,0,1],
                      "phred":["Minimal sgRNA Phred-score",9,0,30],
                      "ram":["RAM saving mode [y/n]",10,0,"n"],
                      "delete":["keep intermediary files [y/n]",11,0,"y"],
                      "upstream":["upstream search sequence",12,0,"None"],
                      "downstream":["downstream search sequence",13,0,"None"],
                      "miss_search":["mismatches in the search sequence",14,0,0],}
    
    dropdown_options = {"extension":[".fastq.gz", ".fastq",5,0],
                        "Running Mode":["Counter", "Extractor + Counter",4,0]}

    # Generating the dropdown browsing buttons
    for arg in dropdown_options:
        dropdown(arg,dropdown_options[arg])
    
    # Generating the file/folder browsing buttons
    for arg in browsing_inputs:
        browsing(arg,browsing_inputs[arg])
    
    # Generating the input parameter buttons
    for arg in default_inputs:
        write_menu(arg,default_inputs[arg])
    
    placeholder(0,1,"",0,0)
    placeholder(15,0,"",0,0)
    button_click(16, 0, "OK", submit)
    button_click(16, 1, "Reset", restart)

    root.mainloop()

    return parameters

def initializer(cmd):
    
    """ Handles the program initialization process.
    Makes sure the path separators, and the input parser function is correct
    for the used OS.
    Creates the output diretory and handles some parameter parsing"""
 
    version = "2.0"
    
    print("\nVersion: {}".format(version))
    
    if system() == 'Windows':
        separator = "\\"
    else:
        separator = "/"

    if cmd is None:
        param = inputs_handler(separator)
    else:
        param = cmd
    
    param["extension"] = f'*{param["extension"]}'
    param["version"] = version
    param["separator"] = separator
    
    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI" #Phred score

    param["quality_set"] = set(quality_list[:int(param['phred'])-1])
    
    current_time = datetime.datetime.now().strftime('%S%M%H%d%m%Y')
    param["directory"] = os.path.join(param['out'], f"output_{current_time}")
    if not os.path.exists(param["directory"]):
        os.makedirs(param["directory"])
    
    if psutil.virtual_memory().percent>=60:
        print("\nLow RAM availability detected, file processing may be slow\n")
    
    print(f"\nRunning with parameters:\n{param['miss']} mismatch allowed\nMinimal Phred Score per bp >= {param['phred']}\n")
    print(f"All data will be saved into {param['directory']}")

    return param

def input_parser():
    
    """ Handles the cmd line interface, and all the parameter inputs"""
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c",nargs='?',const=True,help="cmd line mode")
    parser.add_argument("--s",help="The full path to the directory with the sequencing files")
    parser.add_argument("--g",help="The full path to the .csv file with the sgRNAs.")
    parser.add_argument("--o",help="The full path to the output directory")
    parser.add_argument("--se",help="Sequencing file extenction (ie:'.fastq.gz')")
    parser.add_argument("--m",help="number of allowed mismatches (default=1)")
    parser.add_argument("--ph",help="Minimal Phred-score (default=30)")
    parser.add_argument("--st",help="guideRNA start position in the read (default is 0==1st bp)")
    parser.add_argument("--l",help="guideRNA length (default=20bp)")
    parser.add_argument("--r",help="ram saving mode (only appropriate for mismatch searching)")
    parser.add_argument("--us",help="Upstream search sequence")
    parser.add_argument("--ds",help="Downstream search sequence")
    parser.add_argument("--ms",help="mismatches allowed when searching reads with Up/Down stream sequences")
    parser.add_argument("--mo",help="Running Mode (default=C) [Counter (C) / Extractor + Counter (EC)]")
    parser.add_argument("--k",help="If enabled, keeps all temporary files (default is enabled)")
    args = parser.parse_args()

    if args.c is None:
        return None
    
    parameters = {}
    if (args.s is None) or (args.g is None) or (args.o is None) or (args.se is None):
        print(parser.print_usage())
        raise ValueError("\nPlease specify --s,--g,--o, and --se parameters. type -h into cmd for help")
    else:
        parameters['seq_files'],parameters['sgrna'],parameters['out'],\
        parameters['extension'] = args.s, args.g, args.o, args.se

    parameters['ram']=False
    if args.r is not None:
        parameters['ram']=True
        
    parameters['length']=20
    if args.l is not None:
        parameters['length']=args.l
                 
    parameters['start']=0
    if args.st is not None:
        parameters['start']=args.st
        
    parameters['phred']=30
    if args.ph is not None:
        parameters['phred']=args.ph
                 
    parameters['miss']=1
    if args.m is not None:
        parameters['miss']=args.m
        
    parameters['upstream']=None
    if args.us is not None:
        parameters['upstream']=args.us
        
    parameters['downstream']=None
    if args.ds is not None:
        parameters['downstream']=args.ds
                 
    parameters['miss_search']=0
    if args.ms is not None:
        parameters['miss_search']=args.ms
        
    parameters['Running Mode']="C"
    if args.mo is not None:
        if "EC" in args.mo.upper():
            parameters['Running Mode']="EC"
        
    parameters['delete']=True
    if args.k is not None:
        parameters['delete']=False
        
    return parameters


def compiling(directory,phred,mismatch,version,separator,delete,paths):
    
    def deleting(files):
        for file in files:
            os.remove(file)
    
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
    
    if delete:
        for file in paths:
            os.remove(file)
        for file in ordered_csv:
            os.remove(file[1])
        
    input("\nAnalysis successfully completed\nAll the reads have been compiled into the compiled.csv file.\nPress any key to exit")

def run_stats(headers, out_file,separator):
    
    """ Manipulates the statistics from all the samples into one file that can
    be used for downstream user quality control aplications. Creates a simple
    bar graph with the number of reads per sample"""
    
    ### parsing the stats from the read files
    global_stat = [["#Sample name", "Running Time", "Running Time unit", \
                    "Total number of reads in sample", \
                    "Total number of reads that passed quality control parameters", \
                    "Number of reads that were aligned without mismatches", \
                    "Number of reads that were aligned with mismatches"]]
        
    for run in headers:
        if "script ran" in run:
            parsed = run.split()
            global_stat.append([parsed[7][:-1]] + [parsed[3]] + [parsed[4]] + \
                               [parsed[12]] + [parsed[8]] + [parsed[16]] + \
                               [parsed[20]])
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

def ram_lock():

    """ stalls the program until more RAM is available from finishing the 
    processing of other files """
    
    if psutil.virtual_memory().percent >= 98:
        print("\nSub-optimal RAM allocation detected. Please consider running the program in 'RAM saving mode'\n")
    
    while psutil.virtual_memory().percent >= 98:
        pass
    return True

def cpu_counter():
    
    """ counts the available cpu cores, required for spliting the processing
    of the files """
    
    cpu = multiprocessing.cpu_count()
    if cpu >= 2:
        cpu -= 1
    pool = multiprocessing.Pool(processes = cpu)
    
    return pool

def multi(files,write_path_save,sgrna,param):
    
    """ starts and handles the parallel processing of all the samples by calling 
    multiple instances of the "aligner" function (one per sample) """
    
    pool = cpu_counter()
    for i, (name, out) in enumerate(zip(files, write_path_save)):
        pool.apply_async(aligner, args=(name,out,i,len(files),sgrna,param))
        
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
    param = initializer(input_parser())
    
    ### parses the names/paths, and orders the sequencing files
    ordered = path_finder_seq(param["seq_files"], param["extension"], param["separator"])
    
    ### parses the sequencing files depending on whether they require unzipping or not
    files, write_path_save = input_file_type(ordered, param["extension"], param["directory"])
    
    ### loads the sgRNAs from the input .csv file. 
    ### Creates a dictionary "sgrna" of class instances for each sgRNA
    sgrna = {}
    if param['Running Mode']=='C':
        sgrna = guides_loader(param["sgrna"])
    
    ### Processes all the samples by associating sgRNAs to the reads on the fastq files.
    ### Creates one process per sample, allowing multiple samples to be processed in parallel. 
    multi(files,write_path_save,sgrna,param)
    
    ### Compiles all the processed samples from multi into one file, and creates the run statistics
    compiling(param["directory"],param["phred"],param["miss"],param["version"],param["separator"],param["delete"],write_path_save)
    
if __name__ == "__main__":
    multiprocessing.freeze_support() # required to run multiprocess as .exe on windows
    main()
