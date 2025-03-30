import csv
import glob
import os
import gzip
import multiprocessing as mp
import time
import matplotlib.pyplot as plt
import numpy as np
from numba import njit
import psutil
import argparse
import datetime
from tqdm import tqdm
from dataclasses import dataclass
from pathlib import Path
from colorama import Fore
import importlib.resources as resources

#####################

@dataclass
class Features:
    """
    Data class for storing feature information and read counts.

    Attributes
    ----------
    name : str
        The name of the feature (e.g., guide RNA or sequence identifier).
    counts : int
        The count of reads associated with this feature.

    Usage
    -----
    Each feature instance represents a unique sequence, with its associated read 
    counts stored for tracking occurrences. Instances of this class are typically 
    stored in a dictionary where the sequence serves as the key.

    See Also
    --------
    features_loader : Function that creates instances of this class and stores them in a dictionary.
    """
    name: str
    counts: int

def colourful_errors(warning_type,error):
    """
    Print a color-coded error message with a timestamp.

    Parameters
    ----------
    warning_type : str
        Type of warning (e.g., "FATAL", "WARNING").
    error : str
        Error message to be displayed.

    Returns
    -------
    None
    """
    warning_colour = Fore.GREEN
    if warning_type == "FATAL":
        warning_colour = Fore.RED
    elif warning_type == "WARNING":
        warning_colour = Fore.YELLOW

    print(f"{Fore.BLUE} {datetime.datetime.now().strftime('%c')}{Fore.RESET} [{warning_colour}{warning_type}{Fore.RESET}] {error}")

def path_finder(folder_path,extension):
    """
    Find files with a given extension in a directory and retrieve their sizes.

    Parameters
    ----------
    folder_path : str
        Path to the folder where files should be searched.
    extension : list of str
        List of file extensions (e.g., ["*.csv", "*.txt"]) to search for.

    Returns
    -------
    list of list
        A list where each entry is [file_path, file_size].
    """
    pathing = []
    for exten in extension:
        for filename in glob.glob(os.path.join(folder_path, exten)):
            pathing.append([filename] + [os.path.getsize(filename)]) 
    return pathing

def path_parser(folder_path, extension): 
    """
    Retrieve and sort file paths based on their size or filename.

    Parameters
    ----------
    folder_path : str
        Path to the folder where files should be searched.
    extension : list of str
        List of file extensions (e.g., ["*.csv", "*.txt"]) to search for.

    Returns
    -------
    list of str
        List of sorted file paths, either by size or by filename.

    Raises
    ------
    SystemExit
        If no files matching the given extension are found.
    """
    pathing=path_finder(folder_path,extension)
    if extension != ['*reads.csv']:
        ordered = [path[0] for path in sorted(pathing, key=lambda e: e[-1])]

        if ordered == []:
            colourful_errors("FATAL",f"Check the path to the {extension} files folder. No files of this type found.\n")
            exit()

    else:
        ordered = [path[0] for path in sorted(pathing)]

    return ordered

def features_loader(guides):
    
    """
    Load and parse feature names and sequences from a CSV file.

    Parameters
    ----------
    guides : str
        Path to the CSV file containing feature names and sequences.

    Returns
    -------
    dict
        Dictionary where keys are unique feature sequences (str),
        and values are instances of the `Features` class.

    Raises
    ------
    SystemExit
        If the file is missing or improperly formatted.
        If duplicate sequence entries are found, a warning is issued.
    """
    
    colourful_errors("INFO","Loading Features")
    
    if not os.path.isfile(guides):
        colourful_errors("FATAL",f"Check the path to the features file.\nNo .csv file found in the following path: {guides}\n")
        exit()
    
    features = {}
    names = set()
    
    try:
        with open(guides) as current: 
            for line in current:
                line = line.rstrip().split(",")
                sequence = line[1].upper()
                sequence = sequence.replace(" ", "")
                name = line[0]
                
                if name in names:   
                    colourful_errors("WARNING",f"The name {name} seems to appear at least twice. This MIGHT result in unexpected behaviour. Please have only unique name entries in your features.csv file.")

                if sequence not in features:
                    features[sequence] = Features(name, 0)
                    names.add(name)
                    
                else:
                    colourful_errors("WARNING",f"{features[sequence].name} and {name} share the same sequence. Only {features[sequence].name} will be considered valid. {name} will be ignored.")

    except IndexError:
        colourful_errors("FATAL","The given .csv file doesn't seem to be comma separated. Please double check that the file's column separation is ','\n")
        exit()
    
    colourful_errors("INFO",f"{len(features)} different features were provided.")
    return features

def binary_converter(features):

    """
    Convert DNA sequences into a Numba dictionary with int8 binary arrays for faster computation.

    Parameters
    ----------
    features : dict
        Dictionary of DNA sequences to be converted into binary format.

    Returns
    -------
    numba.typed.Dict
        A Numba dictionary where keys are sequence strings and values are 
        numpy arrays of type int8 representing the binary sequence.
    """
    
    from numba import types
    from numba.typed import Dict
    
    container = Dict.empty(key_type=types.unicode_type,
                            value_type=types.int8[:])

    for sequence in features:
        container[sequence] = seq2bin(sequence)
    return container

def sequence_tinder(read_bin,qual,param,i=0):
    
    """
    Extract barcode positions from a binary-encoded sequencing read based on sequence 
    matching and quality thresholds.

    Parameters
    ----------
    read_bin : numpy.ndarray
        Binary-encoded sequence read (converted using `seq2bin`)
    qual : bytes
        ASCII-encoded quality scores for the read.
    param : dict
        Dictionary containing upstream and downstream sequences (binary format),
        mismatch thresholds, quality filters, and read length constraints.
    i : int, optional
        Index for selecting an upstream/downstream sequence when multiple options exist (default is 0).

    Returns
    -------
    tuple[int, int] or tuple[None, None]
        (start, end) indices for trimming the read, or (None, None) if no valid positions are found.
    """

    start,end = None,None
    if (param['upstream'] is not None) & (param['downstream'] is not None):
        start=border_finder(param['upstream_bin'][i],
                            read_bin,
                            param['miss_search_up'])
        
        if start is not None:
            end=border_finder(param['downstream_bin'][i],
                                read_bin,
                                param['miss_search_down'],
                                start_place=start+len(param['upstream_bin'][i]))

            if end is not None:
                qual_up = str(qual[start:start+len(param['upstream_bin'][i])],"utf-8")
                qual_down = str(qual[end:end+len(param['downstream_bin'][i])],"utf-8")
                
                if (len(param['quality_set_up'].intersection(qual_up)) == 0) &\
                    (len(param['quality_set_down'].intersection(qual_down)) == 0):
                    start+=len(param['upstream_bin'][i])
                    return start,end

    elif (param['upstream'] is not None) & (param['downstream'] is None):
        start=border_finder(param['upstream_bin'][i],
                            read_bin,
                            param['miss_search_up'])
        
        if start is not None:
            qual_up = str(qual[start:start+len(param['upstream_bin'][i])],"utf-8")
            
            if len(param['quality_set_up'].intersection(qual_up)) == 0:
                start+=len(param['upstream_bin'][i])
                end = start + param['length']
                return start,end
        
    elif (param['upstream'] is None) & (param['downstream'] is not None):
        end=border_finder(param['downstream_bin'][i],
                            read_bin,
                            param['miss_search_down'])
        
        if end is not None:
            qual_down = str(qual[end:end+len(param['downstream_bin'][i])],"utf-8")
            
            if len(param['quality_set_down'].intersection(qual_down)) == 0:
                start = end-param['length']
                return start,end

    return None,None


def reads_counter(i,o,raw,features,param,cpu,failed_reads,passed_reads,preprocess=False):
    
    """ Reads the fastq file on the fly. 
    Each read is assumed to be composed of 4 lines, with the sequence being 
    on line 2, and the basepair quality on line 4. 
    Every read is trimmed based on the indicated feature positioning. 
    The quality of the obtained trimmed read is crossed against the indicated
    Phred score for quality control.
    If the read has a perfect match with a feature, the respective feature gets a 
    read increase of 1 (done by calling .counts from the respective feature class 
    from the feature class dictionary).
    If the read doesnt have a perfect match, it is sent for mismatch comparison
    via the "mismatch_search_handler" function.
    """
    
    def progress_bar(i,o,raw,cpu):
        def getuncompressedsize(raw):
            try:
                ext = os.path.splitext(raw)[1]
                if ext == ".gz":
                    with gzip.open(raw, 'rt') as f:
                        return sum(1 for _ in f)
                else:
                    with open(raw, 'rt') as f:
                        return sum(1 for _ in f)
            except EOFError:
                return 0

        if o > cpu:
            current = mp.current_process()
            pos = current._identity[0]#-1
        else:
            pos = i+1
        total_file_size = getuncompressedsize(raw)
        tqdm_text = f"Processing file {i+1} out of {o}"
        return tqdm(total=total_file_size,desc=tqdm_text, position=pos,colour="green",leave=False,ascii=True,unit="lines")

    def fastq_parser(current,features,failed_reads,passed_reads,fixed_start):
    
        reading = []
        mismatch = [n+1 for n in range(param['miss'])]
        perfect_counter, imperfect_counter, non_aligned_counter, reads,quality_failed = 0,0,0,0,0
        ram_clearance=ram_lock()

        if param['miss'] != 0:
            binary_features = binary_converter(features)

        try:
            for line in current:
                quality_failed_flag = np.zeros(param['search_iterations'])
                if (not preprocess) & (param['Progress bar']):
                    pbar.update(1)
                reading.append(line.rstrip())
                
                if len(reading) == 4:
                    
                    full_feature = ""
                    for i in range(param['search_iterations']):
                
                        if not fixed_start:

                            start,end=sequence_tinder(seq2bin(str(reading[1],"utf-8")),
                                                        reading[3],
                                                        param,
                                                        i)
                                
                            if (start is not None) & (end is not None):
                                if end < start: #if the end is not found or found before the start
                                    start=None
                                    quality_failed_flag[i] = 1 #this includes reads without search sequences as failed at quality
                            else:
                                quality_failed_flag[i] = 1

                        if fixed_start:
                            start = param['start_positioning'][i]
                            end = param['end_positioning'][i]

                        if (fixed_start) or (start is not None):
                            seq = str(reading[1][start:end].upper(),"utf-8")
                            quality = str(reading[3][start:end],"utf-8") #convert from bin to str

                            if len(param['quality_set'].intersection(quality)) == 0:
                                full_feature += f":{seq}"
                            else:
                                quality_failed_flag[i] = 1          
                    
                    if full_feature != "":
                        seq = full_feature[1:] #remove the first :
                        if param['Running Mode']=='C':
                            if seq in features:
                                features[seq].counts += 1
                                perfect_counter += 1

                            elif mismatch != []:
                                features,imperfect_counter,failed_reads,passed_reads,non_aligned_counter=\
                                mismatch_search_handler(seq,
                                                        mismatch,
                                                        failed_reads,
                                                        binary_features,
                                                        imperfect_counter,
                                                        features,
                                                        passed_reads,
                                                        ram_clearance,
                                                        non_aligned_counter
                                                        )

                            else:
                                non_aligned_counter += 1

                        else:
                            if seq not in features:
                                features[seq] = Features(seq, 1)
                            else:
                                features[seq].counts += 1
                            perfect_counter += 1
                        
                    if quality_failed_flag.all():
                        quality_failed += 1
                    
                    reading = []
                    reads += 1
                    
                    if reads % 1000000 == 0:
                        ram_clearance=ram_lock()
                    
                    if preprocess:
                        if reads == 10000:
                            return reads,perfect_counter,imperfect_counter,features,failed_reads,passed_reads,0,0
            
            if (not preprocess) & (param['Progress bar']):
                pbar.close()
            
        except EOFError:
            colourful_errors("WARNING",f"{raw} is an incomplete or corrupted gzip file. Only partial processing might have occurred.")
            pass

        return reads,perfect_counter,imperfect_counter,features,failed_reads,passed_reads,non_aligned_counter,quality_failed

    fixed_start = True
    # determining the read trimming starting/ending place
    if (param['upstream'] is None) & (param['downstream'] is None):
        param['start_positioning'] = [int(n) for n in param['start'].split(",")]
        param['end_positioning'] = [int(n)+param['length'] for n in param['start_positioning']]
        param['search_iterations'] = len(param['end_positioning'])
        
    else:
        len_up,len_down = 0,0
        fixed_start = False
        if param['upstream'] is not None:
            param['upstream_bin'] = [seq2bin(n.upper()) for n in param['upstream'].split(",")]
            len_up = len(param['upstream_bin'])
        if param['downstream'] is not None:
            param['downstream_bin'] = [seq2bin(n.upper()) for n in param['downstream'].split(",")]
            len_down = len(param['downstream_bin'])
            
        if (param['downstream'] is not None) & (param['upstream'] is not None):
            if len(param['downstream_bin']) != len(param['upstream_bin']):
                colourful_errors("FATAL",f"Up and Downstream sequences must be submitted in concurrent pairs, separated by ,.\n You submitted {len(param['downstream_bin'])} downstream sequences and {len(param['upstream_bin'])} upstream sequences.")
                exit()
                
        param['search_iterations'] = max(len_up,len_down)

    _, ext = os.path.splitext(raw)
    
    if (not preprocess) & (param['Progress bar']):
        pbar = progress_bar(i,o,raw,cpu)
    
    try:
        if ext == ".gz":
            with gzip.open(raw, "rb") as current:
                return fastq_parser(current,features,failed_reads,passed_reads,fixed_start)
        else:
            with open(raw, "rb")  as current:
                return fastq_parser(current,features,failed_reads,passed_reads,fixed_start)
    except EOFError:
        colourful_errors("WARNING",f"{raw} is an incomplete or corrupted gzip file.")
        return None

def seq2bin(sequence):
    """
    Convert a string sequence into a binary numpy array.

    Parameters
    ----------
    sequence : str
        Input string sequence.

    Returns
    -------
    numpy array
        Binary representation of the sequence in int8 format.
    """
    sequence = bytearray(sequence,'utf8')
    return np.array((sequence), dtype=np.int8)

@njit
def binary_subtract(array1,array2,mismatch):
    """
    Compare two binary sequences and check if they match (Hamming distance) within a given mismatch threshold.

    Parameters
    ----------
    array1 : numpy array
        First sequence in numerical (binary) format.
    array2 : numpy array
        Second sequence in numerical (binary) format.
    mismatch : int
        Maximum number of mismatches allowed.

    Returns
    -------
    int
        1 if the sequences match within the allowed mismatches, otherwise 0.
    """
    miss=0
    for arr1,arr2 in zip(array1,array2):
        if arr1-arr2 != 0:
            miss += 1
        if miss>mismatch:
            return 0
    return 1

@njit
def border_finder(seq,read,mismatch,start_place=0): 
    """
    Find the position where a sequence matches a read within a mismatch threshold.

    Parameters
    ----------
    seq : numpy array
        Sequence to search for (converted to int8 format). Make sure the letters' capitalization matches the one from read.
    read : numpy array
        Read in which the sequence should be found (converted to int8 format). Needs to larger than seq.
    mismatch : int
        Maximum number of mismatches allowed.
    start_place : int, optional
        Start index for the search (default is 0).

    Returns
    -------
    int or None
        Index position where the sequence is found, or None if not found.
    """
    s=seq.size
    r=read.size
    fall_over_index = r-s
    for i,bp in enumerate(read[start_place:]): 
        comparison = read[start_place+i:s+start_place+i]
        finder = binary_subtract(seq,comparison,mismatch)
        if i+start_place > fall_over_index:
            return
        if finder != 0:
            return i+start_place

@njit
def features_all_vs_all(binary_features,read,mismatch):
    """
    Compare a read against all sgRNA sequences to find a match within a mismatch threshold.

    Parameters
    ----------
    binary_features : dict
        Dictionary containing sgRNA sequences in binary format.
    read : numpy array
        Read sequence in binary format.
    mismatch : int
        Maximum number of mismatches allowed.

    Returns
    -------
    str or None
        The matching guide if exactly one match is found, otherwise None.
    """    
    found = 0
    r=read.size
    for guide in binary_features:
        g=binary_features[guide].size
        if g==r:
            if binary_subtract(binary_features[guide],read,mismatch):
                found+=1
                found_guide = guide
                if found>=2:
                    return
    if found==1:
        return found_guide

def mismatch_search_handler(seq,mismatch,failed_reads,binary_features,imperfect_counter,features,passed_reads,ram_clearance,non_aligned_counter):
    """
    Perform an imperfect alignment search on a sequencing read and update feature counts.

    Parameters
    ----------
    seq : str
        The input sequencing read.
    mismatch : list of int
        List of mismatch tolerance values to search against.
    failed_reads : set
        Set of previously failed reads to avoid redundant processing.
    binary_features : dict
        Dictionary of feature sequences in binary format.
    imperfect_counter : int
        Counter for reads that align imperfectly.
    features : dict
        Dictionary mapping feature sequences to `Features` class instances.
    passed_reads : dict
        Dictionary storing successfully aligned reads and their matched features.
    ram_clearance : bool
        If True, failed reads will be stored to conserve memory.
    non_aligned_counter : int
        Counter for non-aligned reads.

    Returns
    -------
    tuple
        (features, imperfect_counter, failed_reads, passed_reads, non_aligned_counter)
        Updated data structures and counters after processing the read.
    """
    if seq in failed_reads:
        non_aligned_counter += 1
        return features,imperfect_counter,failed_reads,passed_reads,non_aligned_counter
    
    if seq in passed_reads:
        features[passed_reads[seq]].counts += 1
        imperfect_counter += 1
        return features,imperfect_counter,failed_reads,passed_reads,non_aligned_counter
    
    read=seq2bin(seq)                         
    for miss in mismatch:

        feature = features_all_vs_all(binary_features, read, miss)
        
        if feature is not None:
            features[feature].counts += 1
            imperfect_counter += 1
            passed_reads[seq] = feature
            return features,imperfect_counter,failed_reads,passed_reads,non_aligned_counter
    
        else: 
            if miss == mismatch[-1]:
                #if function reaches here its because nothing was aligned anywhere
                if ram_clearance:
                    failed_reads.add(seq)
                non_aligned_counter += 1
                return features,imperfect_counter,failed_reads,passed_reads,non_aligned_counter

def aligner(raw,i,o,features,param,cpu,failed_reads,passed_reads):

    """ Runs the main read to sgRNA associating function "reads_counter".
    Creates some visual prompts to alert the user that the samples are being
    processed. Some on the fly quality control is possible (such as making sure 
    the total number of samples is correct, getting an estimate of the total
    number of reads per sample, and checking total running time"""
    
    tempo = time.perf_counter()

    reads, perfect_counter, imperfect_counter, features,failed_reads,passed_reads,non_aligned_counter,quality_failed = \
        reads_counter(i,o,raw,features,param,cpu,failed_reads,passed_reads)

    master_list = []
    [master_list.append([features[guide].name] + [features[guide].counts]) for guide in features]

    tempo = time.perf_counter() - tempo
    if tempo > 3600:
        timing = str(round(tempo / 3600, 2)) + " hours"   
    elif tempo > 60:
        timing = str(round(tempo / 60, 2)) + " minutes"   
    else:
        timing = str(round(tempo, 2)) + " seconds"
    
    name = Path(raw).stem
    path,_ = os.path.splitext(raw)
    if ".fastq" in name:
         name = Path(name).stem
         path,_ = os.path.splitext(path)

    stats_condition = f"#script ran in {timing} for file {name}. {perfect_counter+imperfect_counter} reads out of {reads} were aligned. {perfect_counter} were perfectly aligned. {imperfect_counter} were aligned with mismatch. {non_aligned_counter} passed quality filtering but were not aligned. {quality_failed} did not pass quality filtering."
    
    if not param['Progress bar']:
        colourful_errors("INFO",f"Sample {name} was processed in {timing}")
        
    try:
        master_list.sort(key = lambda master_list: int(master_list[0])) #numerical sorting
    except ValueError:
        master_list.sort(key = lambda master_list: master_list[0]) #alphabetical sorting
    
    master_list.insert(0,["#Feature"] + ["Reads"])
    master_list.insert(0,[stats_condition])
    
    csvfile = os.path.join(param["directory"], name+"_reads.csv")
    csv_writer(csvfile, master_list)

    return failed_reads,passed_reads

def csv_writer(path, outfile):
    
    """ writes the indicated outfile into an .csv file in the directory"""
        
    with open(path, "w", newline='') as output: 
        writer = csv.writer(output)
        writer.writerows(outfile)

def inputs_handler():
    
    """ assertains the correct parsing of the input parameters"""
    
    parameters=inputs_initializer()

    try:

        parameters["length"]=int(parameters["length"])
        parameters["miss"]=int(parameters["miss"])
        parameters["phred"]=int(parameters["phred"])
        
        if parameters['Search Features'] == "Custom":
            parameters["miss_search_up"]=int(parameters["miss_search_up"])
            parameters["miss_search_down"]=int(parameters["miss_search_down"])
            parameters["qual_up"]=int(parameters["qual_up"])
            parameters["qual_down"]=int(parameters["qual_down"])
        else:
            parameters["start"]="0"
            parameters["miss_search_up"]=0
            parameters["miss_search_down"]=0
            parameters["qual_up"]=30
            parameters["qual_down"]=30
            parameters['upstream'] = None
            parameters['downstream'] = None
            
    except Exception:
        colourful_errors("FATAL",f"Please confirm you have provided the correct parameters.\nOnly numeric values are accepted in the folowing fields:\n-Feature read starting place;\n-Feature length;\n-mismatch;\n-Phred score.\n")
        exit()
    
    # avoids getting -1 and actually filtering by highest phred score by mistake
    if int(parameters["phred"]) == 0:
        parameters["phred"] = 1

    if int(parameters["qual_up"]) == 0:
        parameters["qual_up"] = 1
        
    if int(parameters["qual_down"]) == 0:
        parameters["qual_down"] = 1
    
    if parameters['delete'] == "y":
        parameters['delete'] = True
    else:
        parameters['delete'] = False
        
    if parameters['Progress bar'] == "Yes":
        parameters['Progress bar'] = True
    else:
        parameters['Progress bar'] = False
        
    if parameters['upstream'] == "None":
        parameters['upstream'] = None
        
    if parameters['downstream'] == "None":
        parameters['downstream'] = None
        
    if "Extractor" in parameters['Running Mode']:
        parameters['Running Mode']="EC"
    else:
        parameters['Running Mode']="C"
        
    parameters["cmd"] = False
    parameters['cpu'] = False

    return parameters

def inputs_initializer():
    """ Handles the graphical interface and all the parameter inputs """
    
    from tkinter import Entry, LabelFrame, Button, Label, Tk, filedialog, StringVar, OptionMenu, Toplevel, DISABLED, NORMAL
    
    def restart():
        root.quit()
        root.destroy()
        inputs_initializer()
        
    def submit():
        for arg in temporary:
            if isinstance(temporary[arg], Entry):
                value = temporary[arg].get()
                if ("Use the Browse button to navigate, or paste a link" not in value) and (value != ""):
                    parameters[arg] = value
            elif isinstance(temporary[arg], StringVar):
                parameters[arg] = temporary[arg].get()
        root.quit()
        root.destroy()
    
    def directory(column, row, parameter, frame):
        filename = filedialog.askdirectory(title="Select a folder")
        filing_parser(column, row, filename, parameter, frame)
        
    def file(column, row, parameter, frame):
        filename = filedialog.askopenfilename(title="Select a file", filetypes=(("CSV files", "*.csv"), ("all files", "*.*")))
        filing_parser(column, row, filename, parameter, frame)
    
    def filing_parser(column, row, filename, parameter, frame):
        place = Entry(frame, borderwidth=5, width=50)
        place.grid(row=row, column=column+1)
        place.insert(0, filename)
        parameters[parameter] = filename

    def browsing(keyword, inputs, placeholder="Use the Browse button to navigate, or paste a link"):
        title1, title2, row, column, function = inputs
        frame = LabelFrame(root, text=title1, padx=5, pady=5)
        frame.grid(row=row, column=column, columnspan=2)
        button = Button(frame, text=title2, command=lambda: function(column, row, keyword, frame))
        button.grid(column=column, row=row)
        button_place = Entry(frame, borderwidth=5, width=50)
        button_place.grid(row=row, column=column+1, columnspan=2)
        button_place.insert(0, placeholder)
        temporary[keyword] = button_place
        
    def write_menu(keyword, inputs):
        title, row, column, default = inputs
        start = Entry(root, width=25, borderwidth=5)
        start.grid(row=row, column=column+1, padx=5, pady=2)
        start.insert(0, default)
        placeholder(row, column, title, 10, 1)
        temporary[keyword] = start
        return start  # Return the Entry widget for reference
        
    def dropdown(keyword, inputs, command=None):
        default, option, row, column = inputs
        variable = StringVar()
        variable.set(default)
        placeholder(row, column, keyword, 10, 1)
        options = [default] + (option if isinstance(option, list) else [option])
        drop = OptionMenu(root, variable, *options, command=command)
        drop.grid(column=column+1, row=row)
        temporary[keyword] = variable
        return variable  # Return the variable for tracing if needed
        
    def button_click(row, column, title, function):
        button_ok = Button(root, text=title, padx=5, pady=2, width=15, command=function)
        button_ok.grid(row=row, column=column, columnspan=1)
        
    def placeholder(row, column, title, padx, pady):
        label = Label(root, text=title)
        label.grid(row=row, column=column, padx=padx, pady=pady)
        return label
    
    # Callback function to enable/disable "Feature length" Entry widget
    def variable_length_callback(*args):
        value = args[0]
        if value == 'Yes':
            length_entry.config(state=DISABLED)
            variable_len_popup()
        else:
            length_entry.config(state=NORMAL)
    
    # Variable length feature popup
    def variable_len_popup():
        variable_window = Toplevel(root)
        variable_window.title("Variable Length Options")
        variable_window.minsize(300, 200)

        input_labels = {
            "upstream": ["Upstream search sequence", "None"],
            "downstream": ["Downstream search sequence", "None"],
            "miss_search_up": ["Mismatches in the upstream sequence", "0"],
            "miss_search_down": ["Mismatches in the downstream sequence", "0"],
            "qual_up": ["Minimal upstream sequence Phred-score", "30"],
            "qual_down": ["Minimal downstream sequence Phred-score", "30"],
        }

        for i, (key, (label_text, default_value)) in enumerate(input_labels.items()):
            Label(variable_window, text=label_text).grid(row=i, column=0, padx=10, pady=5)
            var = StringVar()
            entry = Entry(variable_window, width=25, textvariable=var)
            entry.grid(row=i, column=1, padx=10, pady=5)
            var.set(default_value)
            temporary[key] = var

        Button(variable_window, text="Submit", command=variable_window.destroy).grid(row=len(input_labels), column=0, columnspan=2, pady=10)
    
    def create_second_menu():
        var = StringVar()
        var.set("No")
        placeholder(7, 0, "Variable length feature?", 10, 1)
        OptionMenu(root, var, "No", "Yes", command=variable_length_callback).grid(column=1, row=7)
        temporary['Variable length feature?'] = var
    
    # Callback function to trigger "Search Features" popup
    def search_features_callback(*args):
        value = temporary['Search Features'].get()
        if value == 'Custom':
            search_features_popup()

    # Popup for Search Features
    def search_features_popup():
        search_window = Toplevel(root)
        search_window.title("Search Features")
        search_window.minsize(400, 300)

        search_features = {
            "start": ["Feature start position in the read", 0],
            "upstream": ["Upstream search sequence", "None"],
            "downstream": ["Downstream search sequence", "None"],
            "miss_search_up": ["Mismatches in the upstream sequence", 0],
            "miss_search_down": ["Mismatches in the downstream sequence", 0],
            "qual_up": ["Minimal upstream sequence Phred-score", 30],
            "qual_down": ["Minimal downstream sequence Phred-score", 30],
        }

        for i, (key, (label_text, default_value)) in enumerate(search_features.items()):
            Label(search_window, text=label_text).grid(row=i, column=0, padx=10, pady=5)
            var = StringVar()
            entry = Entry(search_window, width=25, textvariable=var)
            entry.grid(row=i, column=1, padx=10, pady=5)
            var.set(default_value)
            temporary[key] = var

        Button(search_window, text="Submit", command=search_window.destroy).grid(row=len(search_features), column=0, columnspan=2, pady=10)

    # Create Search Features OptionMenu
    def create_search_features_menu():
        var = StringVar()
        var.set("Default")
        placeholder(10, 0, "Search Features Options", 10, 1)
        OptionMenu(root, var, "Default", "Custom", command=search_features_callback).grid(column=1, row=10)
        temporary['Search Features'] = var
    
    ####
    
    root = Tk()
    root.title("2FAST2Q Input Parameters Window")
    root.geometry("450x500")  # Set the window to a tighter size
    parameters, temporary = {}, {}  
    
    browsing_inputs = {
        "seq_files": ["Path to the .fastq(.gz) files folder", "Browse", 1, 0, directory],
        "feature": ["Path to the features .csv file", "Browse", 2, 0, file],
        "out": ["Path to the output folder", "Browse", 3, 0, directory]
    }
    
    default_inputs = {
        "out_file_name": ["Output File Name", 4, 0, "compiled"],
        "length": ["Feature length", 8, 0, 20],
        "miss": ["Allowed mismatches in features", 5, 0, 1],
        "phred": ["Minimal feature Phred-score in features", 6, 0, 30]}
    
    dropdown_options = {
        "Running Mode": ["Counter", "Extractor + Counter", 9, 0],
        "Progress bar": ["Yes", "No", 11, 0],
        "Delete intermediary files": ["Yes", "No", 12, 0]
    }
    
    # Dropdown for variable length feature
    create_second_menu()
    create_search_features_menu()
    
    # Generating the dropdown browsing buttons
    for arg in dropdown_options:
        dropdown(arg, dropdown_options[arg])
    
    # Generating the file/folder browsing buttons
    [browsing(arg, browsing_inputs[arg]) for arg in browsing_inputs]
    
    # Generating the input parameter buttons
    length_entry = None  # Initialize variable to hold reference to "length" Entry widget
    for arg in default_inputs:
        entry = write_menu(arg, default_inputs[arg])
        if arg == "length":
            length_entry = entry  # Save reference to "length" Entry widget
    
    placeholder(0, 1, "", 0, 0)
    button_click(15, 0, "OK", submit)
    button_click(15, 1, "Reset", restart)
    
    root.mainloop()
    
    parameters["delete"] = parameters.pop("Delete intermediary files")

    return parameters

def initializer(cmd):
    
    """ Handles the program initialization process.
    Makes sure the path separators, and the input parser function is correct
    for the used OS.
    Creates the output diretory and handles some parameter parsing"""
    
    print(f"\n{Fore.RED} Welcome to:{Fore.RESET}\n")
    print(f" {Fore.RED}██████╗{Fore.RESET} ███████╗ █████╗ ███████╗████████╗{Fore.RED}██████╗{Fore.RESET}  ██████╗ ")
    print(f" {Fore.RED}╚════██╗{Fore.RESET}██╔════╝██╔══██╗██╔════╝╚══██╔══╝{Fore.RED}╚════██╗{Fore.RESET}██╔═══██╗")
    print(f" {Fore.RED} █████╔╝{Fore.RESET}█████╗  ███████║███████╗   ██║    {Fore.RED}█████╔╝{Fore.RESET}██║   ██║")
    print(f" {Fore.RED}██╔═══╝ {Fore.RESET}██╔══╝  ██╔══██║╚════██║   ██║   {Fore.RED}██╔═══╝ {Fore.RESET}██║▄▄ ██║")
    print(f" {Fore.RED}███████╗{Fore.RESET}██║     ██║  ██║███████║   ██║   {Fore.RED}███████╗{Fore.RESET}╚██████╔╝")
    print(f" {Fore.RED}╚══════╝{Fore.RESET}╚═╝     ╚═╝  ╚═╝╚══════╝   ╚═╝   {Fore.RED}╚══════╝{Fore.RESET} ╚══▀▀═╝ ")
                                                          
    print(f"\n{Fore.GREEN} Version: {Fore.RESET}{version}")
    
    param = inputs_handler() if cmd is None else cmd
    
    if "test_mode" in param:
        colourful_errors("WARNING",f"Running test mode!\n")
    
    if (param["upstream"] == None) or (param["downstream"] == None):
        if "Variable length feature?" in param:
            if param["Variable length feature?"] == "Yes":
                colourful_errors("WARNING",f"You have selected variable length features, but such mode is only implemented when 2 search sequences are provided. Only features {param['length']}bp long will be considered.\n")
                exit()
        
    param["version"] = version
    
    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI" #Phred score in order of probabilities

    param["quality_set"] = set(quality_list[:int(param['phred'])-1])
    param["quality_set_up"] = set(quality_list[:int(param['qual_up'])-1])
    param["quality_set_down"] = set(quality_list[:int(param['qual_down'])-1])
    
    current_time = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    param["directory"] = os.path.join(param['out'], f"2FAST2Q_output_{current_time}")

    if psutil.virtual_memory().percent>=75:
        colourful_errors("WARNING","Low RAM availability detected, file processing may be slow\n")
    
    print(f"\n{Fore.YELLOW} -- Parameters -- {Fore.RESET}")
    
    if param['Running Mode']=='C':
        print(f"\n {Fore.GREEN}Mode: {Fore.RESET}Align and count")
        print(f" {Fore.GREEN}Allowed mismatches per alignement: {Fore.RESET}{param['miss']}")
        print(f" {Fore.GREEN}Minimal Phred Score per bp >= {Fore.RESET}{param['phred']}")
        if (param['downstream'] is None) or (param['upstream'] is None):
            print(f" {Fore.GREEN}Feature length: {Fore.RESET}{param['length']}")
            print(f" {Fore.GREEN}Read alignment start position: {Fore.RESET}{param['start']}\n")
    else:
        print(f"\n {Fore.GREEN}Mode: {Fore.RESET}Extract and count")
        print(f" {Fore.GREEN}Minimal Phred Score per bp >= {Fore.RESET}{param['phred']}")

        if param['upstream'] is not None:
            print(f" {Fore.GREEN}Upstream search sequence: {Fore.RESET}{param['upstream']}")
            print(f" {Fore.GREEN}Mismatches allowed in the upstream search sequence: {Fore.RESET}{param['miss_search_up']}")
            print(f" {Fore.GREEN}Minimal Phred-score in the upstream search sequence: {Fore.RESET}{param['qual_up']}")
            
        if param['downstream'] is not None:
            print(f" {Fore.GREEN}Downstream search sequence: {Fore.RESET}{param['downstream']}")
            print(f" {Fore.GREEN}Mismatches allowed in the downstream search sequence: {Fore.RESET}{param['miss_search_down']}")
            print(f" {Fore.GREEN}Minimal Phred-score in the downstream search sequence: {Fore.RESET}{param['qual_down']}")

        if (param['upstream'] is None) or (param['downstream'] is None):
            print(f" {Fore.GREEN}Finding features with the folowing length: {Fore.RESET}{param['length']}bp\n")

    print(f" {Fore.GREEN}All data will be saved into {Fore.RESET}{param['directory']}")
    print(f"\n{Fore.YELLOW} ---- {Fore.RESET}")
    
    return param

def input_parser():
    
    """ Handles the cmd line interface, and all the parameter inputs"""
    
    global version
    version = "2.7.7"
    
    def current_dir_path_handling(param):
        if param[0] is None:
            parameters[param[1]]=os.getcwd()
            if param[1] == 'feature':
                file = path_finder(os.getcwd(), ["*.csv"])
                if parameters['Running Mode']!="EC":
                    if len(file) > 1:
                        colourful_errors("FATAL","There is more than one .csv in the current directory. If not directly indicating a path for the features .csv, please have only 1 .csv file in the directory.\n")
                        exit()
                    if len(file) == 1:
                        parameters[param[1]]=file[0][0]
        else:
            parameters[param[1]]=param[0]
        return parameters
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c",nargs='?',const=True,help="cmd line mode")
    parser.add_argument("-t",nargs='?',const=True,help="Runs 2FAST2Q in test mode with example data.")
    parser.add_argument("-v",nargs='?',const=True,help="prints the installed version")
    parser.add_argument("--s",help="The full path to the directory with the sequencing files OR file")
    parser.add_argument("--g",help="The full path to the .csv file with the sgRNAs.")
    parser.add_argument("--o",help="The full path to the output directory")
    parser.add_argument("--fn",nargs='?',const="compiled",help="Specify an output compiled file name (default is called compiled)")
    parser.add_argument("--pb",nargs='?',const=False,help="Adds progress bars (default is enabled)")
    parser.add_argument("--m",help="number of allowed mismatches (default=1)")
    parser.add_argument("--ph",help="Minimal Phred-score (default=30)")
    parser.add_argument("--st",help="Feauture start position in the read (default is 0==1st bp)")
    parser.add_argument("--l",help="Feature length (default=20bp)")
    parser.add_argument("--us",help="Upstream search sequence")
    parser.add_argument("--ds",help="Downstream search sequence")
    parser.add_argument("--msu",help="mismatches allowed in the upstream sequence")
    parser.add_argument("--msd",help="mismatches allowed in the downstream sequence")
    parser.add_argument("--qsu",help="Minimal Phred-score (default=30) in the upstream search sequence")
    parser.add_argument("--qsd",help="Minimal Phred-score (default=30) in the downstream search sequence")
    parser.add_argument("--mo",help="Running Mode (default=C) [Counter (C) / Extractor + Counter (EC)]")
    parser.add_argument("--cp",help="Number of cpus to be used (default is max(cpu)-2 for >=3 cpus, -1 for >=2 cpus, 1 if 1 cpu")
    parser.add_argument("--k",nargs='?',const=False,help="If enabled, keeps all temporary files (default is disabled)")
    args = parser.parse_args()
    
    if args.v is not None:
        print(f"\nVersion: {version}\n")
        exit()
    
    #if its not running on command window mode
    if args.c is None:
        return None
    
    parameters = {}
    parameters["cmd"] = True
    if args.t is None:
        paths_param = [[args.s,'seq_files'],
                       [args.g,'feature'],
                       [args.o,'out']]
    else:
        parameters["test_mode"] = True
        paths_param = [[resources.files("fast2q").joinpath('data/example.fastq.gz'), 'seq_files'],
                        [resources.files("fast2q").joinpath('data/D39V_guides.csv'), 'feature'],
                        [os.getcwd(), 'out']]
       
    parameters['out_file_name'] = "compiled"
    if args.fn is not None:
        parameters['out_file_name'] = args.fn

    parameters['length']=20
    if args.l is not None:
        parameters['length']=int(args.l)
        
    parameters['Progress bar']=True
    if args.pb is not None:
        parameters['Progress bar']=False
                 
    parameters['start']="0"
    if args.st is not None:
        parameters['start']=args.st
        
    parameters['phred']=30
    if args.ph is not None:
        parameters['phred']=int(args.ph)
        if int(parameters["phred"]) == 0:
            parameters["phred"] = 1
                 
    parameters['miss']=1
    if args.m is not None:
        parameters['miss']=int(args.m)
        
    parameters['upstream']=None
    if args.us is not None:
        parameters['upstream']=args.us
        
    parameters['downstream']=None
    if args.ds is not None:
        parameters['downstream']=args.ds
                 
    parameters['miss_search_up']=0
    if args.msu is not None:
        parameters['miss_search_up']=int(args.msu)
        
    parameters['miss_search_down']=0
    if args.msd is not None:
        parameters['miss_search_down']=int(args.msd)
        
    parameters['qual_up']=30
    if args.qsu is not None:
        parameters['qual_up']=int(args.qsu)
        
    parameters['qual_down']=30
    if args.qsd is not None:
        parameters['qual_down']=int(args.qsd)
        
    parameters['Running Mode']="C"
    if args.mo is not None:
        if "EC" in args.mo.upper():
            parameters['Running Mode']="EC"
        
    parameters['delete']=True
    if args.k is not None:
        parameters['delete']=False
        
    parameters['cpu']=False
    if args.cp is not None:
        parameters['cpu']=int(args.cp)
        
    for param in paths_param:
        parameters = current_dir_path_handling(param)
        
    return parameters

def compiling(param):

    """ Combines all the individual processed .csv files into one final file.
    Gathers the individual sample statistic and parses it into "run_stats" """
    
    ordered_csv = path_parser(param["directory"], ['*reads.csv'])

    headers = [f"#2FAST2Q version: {param['version']}"] + \
            [f"#Mismatch: {param['miss']}"] + \
            [f"#Phred Score: {param['phred']}"] + \
            [f"#Feature Length: {param['length']}"] + \
            [f"#Feature start position in the read: {param['start']}"] + \
            [f"#Running mode: {param['Running Mode']}"] + \
            [f"#Upstream search sequence: {param['upstream']}"] + \
            [f"#Downstream search sequence: {param['downstream']}"] + \
            [f"#Mismatches in the upstream search sequence: {param['miss_search_up']}"] + \
            [f"#Mismatches in the downstream search sequence: {param['miss_search_down']}"] + \
            [f"#Minimal Phred-score in the upstream search sequence: {param['qual_up']}"] + \
            [f"#Minimal Phred-score in the downstream search sequence: {param['qual_down']}"]
            
    headers = headers[::-1]

    compiled = {} #dictionary with all the reads per feature
    head = ["#Feature"] #name of the samples
    for i, file in enumerate(ordered_csv):
        path,_ = os.path.splitext(file)
        path = Path(path).stem
        path = path[:-len("_reads")]
        head.append(path)
        with open(file) as current:
            for line in current:
                line = line.rstrip().split(",")
                if "#" not in line[0]:
                    
                    if line[0] in compiled: 
                        compiled[line[0]] = compiled[line[0]] + [int(line[1])]
                    else:
                        compiled[line[0]] = [0]*i + [int(line[1])]
                        
                elif "#Feature" not in line[0]:
                    headers.append(line[0][1:]+"\n")
    
        #important in extract and count mode for creating entries with 0 reads
        for entry in compiled:
            if len(compiled[entry])<i+1:
                compiled[entry] = compiled[entry] + [0]*(i+1-len(compiled[entry]))

    run_stats(headers,param,compiled,head)

    final = []
    [final.append([feature] + compiled[feature]) for feature in compiled] 
    final.insert(0, head)
    
    csvfile = os.path.join(param["directory"],f"{param['out_file_name']}.csv")
    csv_writer(csvfile, final)

    if param["delete"]:
        for file in ordered_csv:
            os.remove(file)
            
    colourful_errors("INFO","Analysis successfully completed")

    print(f"\n{Fore.GREEN} If you find 2FAST2Q useful, please consider citing:{Fore.MAGENTA}\n Bravo AM, Typas A, Veening J. 2022. \n 2FAST2Q: a general-purpose sequence search and counting program for FASTQ files. PeerJ 10:e14041\n DOI: 10.7717/peerj.14041\n{Fore.RESET}")

    if "test_mode" in param:
        colourful_errors("WARNING","Test successful. 2FAST2Q is working as intended!\n")

def run_stats(headers, param, compiled, head):
    
    """ Manipulates the statistics from all the samples into one file that can
    be used for downstream user quality control aplications. Creates a simple
    bar graph with the number of reads per sample"""
    
    ### parsing the stats from the read files
    global_stat = [["#Sample name", "Running Time", "Running Time unit", \
                    "Total number of reads in sample", \
                    "Total number of reads that were aligned", \
                    "Number of reads that were aligned without mismatches", \
                    "Number of reads that were aligned with mismatches",\
                    "Number of reads that passed quality filtering but were not aligned",\
                    'Number of reads that did not pass quality filtering.']]
    
    header_ofset = 1
    for run in headers:
        if "script ran" in run:
            parsed = run.split()
            global_stat.append([parsed[7][:-1]] + [parsed[3]] + [parsed[4]] + \
                               [parsed[12]] + [parsed[8]] + [parsed[15]] + \
                               [parsed[19]] + [parsed[24]] + [parsed[32]])
        else:
            global_stat.insert(0,[run])
            header_ofset+=1
            
    csvfile = os.path.join(param["directory"],f"{param['out_file_name']}_stats.csv")
    csv_writer(csvfile, global_stat)
    
    ######## for bar plots with absolute number of reads
    
    fig, ax = plt.subplots(figsize=(12, int(len(global_stat)/4)))
    width = .75
    for i, (_,_,_,total_reads,aligned,_,_,not_aligned,_) in enumerate(global_stat[header_ofset:]):   
        
        plt.barh(i, int(total_reads), width,  capsize=5, color = "#FFD25A", hatch="//",edgecolor = "black", linewidth = .7)
        plt.barh(i, int(aligned), width,  capsize=5, color = "#FFAA5A", hatch="\\",edgecolor = "black", linewidth = .7)
        plt.barh(i, int(not_aligned), width, capsize=5, color = "#F56416", hatch="x",edgecolor = "black", linewidth = .7)
    
    ax.set_yticks(np.arange(len([n[0] for n in global_stat[header_ofset:]])))
    ax.set_yticklabels([n[0] for n in global_stat[header_ofset:]])
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=16)
    plt.xlabel('Number of reads',size=20)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.set_xscale('log')
    ax.set_xlim(xmin=1)
    ax.legend(["Total reads in sample", "Aligned reads","Reads that passed quality filtering but failed to align"], \
              loc='right',bbox_to_anchor=(1.1, 1),ncol=3,prop={'size': 12})
    plt.tight_layout()
    file = os.path.join(param["directory"],f"{param['out_file_name']}_reads_plot.png")
    plt.savefig(file, dpi=300, bbox_inches='tight')
    
    ######## for bar plots with relative (percentage) number of reads
    
    fig, ax = plt.subplots(figsize=(12, int(len(global_stat)/4)))
    width = .75

    for i, (_,_,_,total_reads,aligned,_,_,not_aligned,q_failed) in enumerate(global_stat[header_ofset:]):   

        aligned = int(aligned)/int(total_reads)*100
        not_aligned = int(not_aligned)/int(total_reads)*100
        q_failed = int(q_failed)/int(total_reads)*100
        
        plt.barh(i, aligned, width,  capsize=5, color = "#6290C3", hatch="\\",edgecolor = "black", linewidth = .7)
        plt.barh(i, not_aligned, width, capsize=5, left=aligned, color = "#F1FFE7", hatch="//",edgecolor = "black", linewidth = .7)
        plt.barh(i, q_failed, width, capsize=5, left=not_aligned+aligned, color = "#FB5012", hatch="||",edgecolor = "black", linewidth = .7)
    
    ax.set_yticks(np.arange(len([n[0] for n in global_stat[header_ofset:]])))
    ax.set_yticklabels([n[0] for n in global_stat[header_ofset:]])
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=16)
    plt.xlabel('% of reads per sample',size=20)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.set_xscale('log')
    ax.set_xlim(xmin=1)
    ax.legend(["Aligned reads","Reads that passed quality filtering but failed to align","Reads that did not pass quality filtering"], \
              loc='right',bbox_to_anchor=(1.1, 1),ncol=3,prop={'size': 12})
    plt.tight_layout()
    file = os.path.join(param["directory"],f"{param['out_file_name']}_reads_plot_percentage.png")
    plt.savefig(file, dpi=300, bbox_inches='tight')
    
    ########
    """ For plotting a violin plot distribution """
    
    distributions = {}
    for feature in compiled:
        for i,read in enumerate(compiled[feature]):
            if head[i+1] in distributions:
                distributions[head[i+1]].append(read)
            else:
                distributions[head[i+1]] = [read]
    
    def violin(data,head,normalized=False):
        fig, ax = plt.subplots(figsize=(12, int(len(global_stat))/2))
        if not normalized:
            ax.set_title('Reads per feature distribution',size=20)
        else:
            ax.set_title('Reads per feature (RPM normalized) distribution',size=20)
        plt.xlabel('Reads per feature',size=20)
        parts=ax.violinplot(data, points=200, widths=1, showmeans=False, showmedians=False,showextrema=False,vert=False)
        
        for pc in parts['bodies']:
            pc.set_facecolor('#D43F3A')
            pc.set_edgecolor('black')
            pc.set_alpha(1)
        
        quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
        inds = np.arange(1, len(medians) + 1)
        ax.scatter(medians,inds, marker='o', color='white', s=40, zorder=3)
        ax.hlines(inds,quartile1, quartile3, color='k', linestyle='-', lw=8)
        #ax.hlines(inds,whiskers_min, whiskers_max, color='k', linestyle='-', lw=2.5)
        ax.set_yticks(np.arange(len(head[1:]))+1)
        ax.set_yticklabels(head[1:])
        plt.tick_params(axis='y', which='major', labelsize=20)
        plt.tick_params(axis='x', which='major', labelsize=20)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #ax.set_xscale('log')
        ax.set_xlim(xmin=1)
        if not normalized:
            file = os.path.join(param["directory"],f"{param['out_file_name']}_distribution_plot.png")
        else:
            file = os.path.join(param["directory"],f"{param['out_file_name']}_distribution_normalized_RPM_plot.png")
        plt.savefig(file, dpi=300, bbox_inches='tight')
    
    data = []
    for entry in distributions:
        data.append(distributions[entry])
    violin(data,head)
    
    ## normalized RPM
    try:
        data = np.array(data)
        data1 = []
        for i,entry in enumerate(data):
            if sum(entry)>0:
                data1.append(entry/sum(entry)*1000000) #RPM
        violin(data1,head,normalized=True)
    except ValueError:
        pass

def ram_lock():
    """
    Check system memory usage and prevent excessive RAM consumption.

    Returns
    -------
    bool
        False if system memory usage is at or above 95%, otherwise True.
    """
    if psutil.virtual_memory().percent >= 95:
        return False
    return True

def cpu_counter(cpu):
    """
    Determine the number of available CPU cores and create a multiprocessing pool.

    Parameters
    ----------
    cpu : int or None
        Desired number of CPU cores for processing. If not provided or invalid, 
        the function automatically selects an optimal number.

    Returns
    -------
    tuple
        (mp.Pool, int) - A multiprocessing pool configured with the selected number of CPUs, 
        and the final CPU count used.
    """
    available_cpu = mp.cpu_count()
    if type(cpu) is not int:
        cpu = available_cpu
        if cpu >= 3:
            cpu -= 2
        if cpu == 2:
            cpu -= 1
    else:
        if cpu > available_cpu:
            cpu = available_cpu

    pool = mp.Pool(processes = cpu, 
                    initargs=(mp.RLock(),), 
                    initializer=tqdm.set_lock)
    
    return pool,cpu

def multiprocess_merger(start,failed_reads,passed_reads,files,features,param,pool,cpu,preprocess=False):
    
    """ runs all samples in blocks equal to the amount of cpus in the PC. 
    Before starting a new block, the hash tables for the failed and passed reads 
    are updated. This confers speed advantages for the next block. """ 
    
    result = []
    for i, name in enumerate(files[start:start+cpu]):
        
        result.append(pool.apply_async(aligner, args=(name,\
                                                      i+start,\
                                                      len(files),\
                                                      features,\
                                                      param,\
                                                      cpu,\
                                                      failed_reads,\
                                                      passed_reads)))
    
    if param["miss"] != 0:
        compiled = [x.get() for x in result]
        failed_reads_compiled,passed_reads_compiled = zip(*compiled)
        return hash_reads_parsing(failed_reads_compiled,passed_reads_compiled,failed_reads,passed_reads)
    
    else:
        return failed_reads,passed_reads

def hash_reads_parsing(failed_reads_compiled,passed_reads_compiled,failed_reads,passed_reads):
    
    """ parsed the results from all the processes, merging the individual failed and
    passed reads into one master file, that will be subsquently used for the new
    samples's processing """ 
    
    for failed,passed in zip(failed_reads_compiled,passed_reads_compiled):
        failed_reads = set.union(failed_reads,failed)
        passed_reads = {**passed_reads,**passed}

    return failed_reads,passed_reads

def hash_preprocesser(files,features,param,pool,cpu):
    
    """ For the smallest files, we processe them for the first x amount of reads to
    initialize the failed reads and passed reads hash tables. This will confer some 
    speed advantages, as subsquent files normally share the same reads that dont
    align to anything, or reads with mismatches that indeed align"""
    
    colourful_errors("INFO","Please standby for the initialization procedure.")
    result=[]

    for name in files[:cpu]:
        result.append(pool.apply_async(reads_counter, args=(0,0,name,features,param,cpu,set(),{},True)))
    
    compiled = [x.get() for x in result]

    _,_,_,_,failed_reads_compiled,passed_reads_compiled,_,_ = zip(*compiled)

    return hash_reads_parsing(failed_reads_compiled,passed_reads_compiled,set(),{})

def aligner_mp_dispenser(files,features,param):
    
    """ starts and handles the parallel processing of all the samples by calling 
    multiple instances of the "aligner" function (one per sample) """
    
    if not os.path.exists(param["directory"]):
        os.makedirs(param["directory"])
    
    start,failed_reads,passed_reads = 0,set(),{}
    pool,cpu = cpu_counter(param["cpu"])
    
    if param["miss"] != 0:
        failed_reads,passed_reads=hash_preprocesser(files,features,param,pool,cpu)
    
    colourful_errors("INFO",f"Processing {len(files)} files. Please hold.")

    for start in range(0,len(files),cpu):
        failed_reads,passed_reads = \
        multiprocess_merger(start,failed_reads,passed_reads,files,\
                            features,param,pool,cpu)
        
    pool.close()
    pool.join()

def main():
    
    """ Runs the program by calling all the appropriate functions"""
    
    ### parses all inputted parameters
    param = initializer(input_parser())
    
    ### parses the names/paths, and orders the sequencing files
    files = [param["seq_files"]]
    if os.path.splitext(param["seq_files"])[1] == '':
        files = path_parser(param["seq_files"], ["*.gz","*.fastq"])

    ### loads the features from the input .csv file. 
    ### Creates a dictionary "feature" of class instances for each sgRNA
    features = {}
    if param['Running Mode']=='C':
        features = features_loader(param["feature"])
    
    ### Processes all the samples by associating sgRNAs to the reads on the fastq files.
    ### Creates one process per sample, allowing multiple samples to be processed in parallel. 
    aligner_mp_dispenser(files,features,param)
    
    ### Compiles all the processed samples from multi into one file, and creates the run statistics
    compiling(param)
    
if __name__ == "__main__":
    mp.freeze_support() # required to run multiprocess as .exe on windows
    main()
