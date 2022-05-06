try:
    import csv
    import glob
    import os
    import gzip
    import shutil
    from concurrent.futures import ThreadPoolExecutor
    import multiprocessing as mp
    from platform import system
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
    import struct
    import tkinter #(imported inside inputs_initializer()). imported here just to assertain
except ImportError as e:
    input(f"{e}\nConsider installing it using 'pip install {e.name}', and try again.\nPress enter to exit")
    raise Exception

#####################

@dataclass
class Features:
    
    """ Each feature will have its own class instance, where the read counts will be kept.
    Each feature class is stored in a dictionary with the name as its key. 
    See the "guides loader" function """   
    
    name: str
    counts: int

def path_finder(folder_path,extension):
    
    """ Finds the correct file paths from the indicated directories,
    and parses the file names into a list for later use"""
    
    pathing = []
    for exten in extension:
        for filename in glob.glob(os.path.join(folder_path, exten)):
            pathing.append([filename] + [os.path.getsize(filename)]) 
    return pathing

def path_parser(folder_path, extension): 
    
    """ parses the file names and paths into an ordered list for later use"""

    pathing=path_finder(folder_path,extension)
    if extension != '*reads.csv':
        
        """sorting by size makes multiprocessing more efficient 
            as the bigger files will be ran first, thus maximizing processor queing """
        
        ordered = [path[0] for path in sorted(pathing, key=lambda e: e[-1])]#[::-1]

        if ordered == []:
            input(f"Check the path to the {extension[1:]} files folder. No files of this type found.\nPress enter to exit")
            raise Exception

    else:
        ordered = [path[0] for path in sorted(pathing, reverse = False)]

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

def features_loader(guides):
    
    """ parses the features names and sequences from the indicated features .csv file.
    Creates a dictionary using the feature sequence as key, with an instace of the 
    Features class as respective value. If duplicated features sequences exist, 
    this will be caught in here"""
    
    print("\nLoading Features")
    
    if not os.path.isfile(guides):
        input("\nCheck the path to the sgRNA file.\nNo file found in the following path: {}\nPress enter to exit".format(guides))
        raise Exception
    
    features = {}
    
    try:
        with open(guides) as current: 
            for line in current:
                if "\n" in line:
                    line = line[:-1]
                line = line.split(",")
                sequence = line[1].upper()
                sequence = sequence.replace(" ", "")
                
                if sequence not in features:
                    features[sequence] = Features(line[0], 0)
                    
                else:
                    print("\nWarning!!\n{} and {} share the same sequence. Only {} will be considered valid.".format(features[sequence].name, line[0],features[sequence].name))
    except IndexError:
        input("\nThe given .csv file doesn't seem to be comma separated. Please double check that the file's column separation is ','\nPress enter to exit")
        raise Exception
        
    print(f"\n{len(features)} different features were provided.")
    return features

def reads_counter(i,o,raw,features,param,cpu,failed_reads,passed_reads,preprocess=False):
    
    """ Reads the fastq file on the fly to avoid RAM issues. 
    Each read is assumed to be composed of 4 lines, with the sequense being 
    on line 2, and the basepair quality on line 4. 
    Every read is trimmed based on the indicated feature positioning. 
    The quality of the obtained trimmed read is crossed against the indicated
    Phred score for quality control.
    If the read has a perfect match with a feature, the respective feature gets a 
    read increase of 1 (done by calling .counts from the respective feature class 
    from the feature class dictionary).
    If the read doesnt have a perfect match, it is sent for mismatch comparison
    via the "imperfect_alignment" function.
    """

    def binary_converter(features):

        """ Parses the input features into numba dictionaries with int8 layout. 
        Converts all DNA sequences to their respective binary array forms. 
        This gives some computing speed advantages with mismatches."""
        
        from numba import types
        from numba.typed import Dict
        
        container = Dict.empty(key_type=types.unicode_type,
                               value_type=types.int8[:])

        for sequence in features:
            container[sequence] = seq2bin(sequence)
        return container
    
    def unfixed_starting_place_parser(read,upstream,downstream,mismatch,lenght):
        
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
                end = start + lenght
            
        elif (upstream is None) & (downstream is not None):
            end=border_finder(downstream,read_bin,mismatch)
            if end is not None:
                start = end-lenght

        return start,end
    
    def progress_bar(i,o,raw,cpu):
        
        def getuncompressedsize(raw):
            # only works in files < 4gb compressed
            with open(raw, 'rb') as f:
                f.seek(-4, 2)
                return struct.unpack('I', f.read(4))[0]
        
        if o > cpu:
            current = mp.current_process()
            pos = current._identity[0]#-1
        else:
            pos = i+1
        total_file_size = getuncompressedsize(raw)
        tqdm_text = f"Processing file {i+1} out of {o}"
        return tqdm(total=total_file_size,desc=tqdm_text, position=pos,colour="green",leave=False,ascii=True,unit="characters")
    
    def fastq_parser(current,end,start,features,failed_reads,passed_reads):
        quality_set = param['quality_set']
        reading = []
    
        mismatch = [n+1 for n in range(param['miss'])]
        perfect_counter, imperfect_counter, reads = 0,0,0
        ram_clearance=ram_lock()

        for line in current:
            if not preprocess:
                pbar.update(len(line))
            reading.append(line[:-1])
            
            if len(reading) == 4: #a read always has 4 lines
                
                if not fixed_start:
                    start,end=unfixed_starting_place_parser(reading[1],upstream,downstream,mismatch_search,lenght)
                
                if (fixed_start) or ((start is not None) & (end is not None)):
                    seq = reading[1][start:end].upper()
                    quality = str(reading[3][start:end],"utf-8") #convert from bin to str
    
                    if (len(quality_set.intersection(quality)) == 0):
                        
                        if param['Running Mode']=='C':
                            if seq in features:
                                features[seq].counts += 1
                                perfect_counter += 1
                            
                            elif mismatch != []:
                                features,imperfect_counter,failed_reads,passed_reads=\
                                mismatch_search_handler(seq,mismatch,\
                                                        failed_reads,binary_features,\
                                                        imperfect_counter,features,\
                                                        passed_reads,ram_clearance)
                        else:
                            if seq not in features:
                                features[seq] = Features(seq, 1)
                            else:
                                features[seq].counts += 1
                                perfect_counter += 1
                reading = []
                reads += 1
                
                # keeps RAM under control by avoiding overflow
                if reads % 1000000 == 0:
                    ram_clearance=ram_lock()
                
                if preprocess:
                    if reads == 200000:
                        return reads,perfect_counter,imperfect_counter,features,failed_reads,passed_reads
        
        if not preprocess:
            pbar.close()

        return reads,perfect_counter,imperfect_counter,features,failed_reads,passed_reads
    
    if param['miss'] != 0:
        binary_features = binary_converter(features)
        
    fixed_start = True
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
        upstream,downstream,mismatch_search,lenght = \
        param['upstream'],param['downstream'],param['miss_search'],param['length']
    
    _, ext = os.path.splitext(raw)

    if not preprocess:
        pbar = progress_bar(i,o,raw,cpu)
    
    if ext == ".gz":
        with gzip.open(raw, "rb") as current:
            return fastq_parser(current,end,start,features,failed_reads,passed_reads)
    else:
        with open(raw, "rb")  as current:
            return fastq_parser(current,end,start,features,failed_reads,passed_reads)

def seq2bin(sequence):
    
    """ Converts a string to binary, and then to 
    a numpy array in int8 format"""
    
    if type(sequence) != bytes:
        sequence = bytearray(sequence,'utf8')
    else:
        sequence = bytearray(sequence)
    return np.array((sequence), dtype=np.int8)

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
def border_finder(seq,read,mismatch): 
    
    """ Matches 2 sequences (after converting to int8 format)
    based on the allowed mismatches. Used for sequencing searching
    a start/end place in a read"""
    
    s=seq.size
    r=read.size
    fall_over_index = r-s-1
    for i,bp in enumerate(read): #range doesnt exist in njit
        comparison = read[i:s+i]
        finder = binary_subtract(seq,comparison,mismatch)
        if i > fall_over_index:
            return
        if finder != 0:
            return i

@njit
def features_all_vs_all(binary_features,read,mismatch):
    
    """ Runs the loop of the read vs all sgRNA comparison.
    Sends individually the sgRNAs for comparison.
    Returns the final mismatch score"""
    
    found = 0
    for guide in binary_features:
        if binary_subtract(binary_features[guide],read,mismatch):
            found+=1
            found_guide = guide
            if found>=2:
                return
    if found==1:
        return found_guide
    return #not needed, but here for peace of mind

def mismatch_search_handler(seq,mismatch,failed_reads,binary_features,imperfect_counter,features,passed_reads,ram_clearance):
    
    """Converts a read into numpy int 8 form. Runs the imperfect alignment 
    function for all number of inputed mismatches."""
    
    read=seq2bin(seq)                         
    for miss in mismatch:

        # we already know this read is going to pass
        if seq in passed_reads:
            features[passed_reads[seq]].counts += 1
            imperfect_counter += 1
            break
        
        elif seq not in failed_reads:
            features,imperfect_counter,feature=\
                imperfect_alignment(read,binary_features,\
                                    miss, imperfect_counter,\
                                    features)
            if ram_clearance:
                if feature is None:
                    failed_reads.add(seq)
                    break
                else:
                    passed_reads[seq] = feature

    return features,imperfect_counter,failed_reads,passed_reads

def imperfect_alignment(read,binary_features, mismatch, counter, features):
    
    """ for the inputed read sequence, this compares if there is a feature 
    with a sequence that is similar to it, to the indicated mismatch degree"""

    feature = features_all_vs_all(binary_features, read, mismatch)
    
    if feature is not None:
        features[feature].counts += 1
        counter += 1

    return features,counter,feature

def aligner(raw,i,o,features,param,cpu,failed_reads,passed_reads):

    """ Runs the main read to sgRNA associating function "reads_counter".
    Creates some visual prompts to alert the user that the samples are being
    processed. Some on the fly quality control is possible (such as making sure 
    the total number of samples is correct, getting an estimate of the total
    number of reads per sample, and checking total running time"""
    
    tempo = time.perf_counter()

    reads, perfect_counter, imperfect_counter, features,failed_reads,passed_reads = \
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

    stats_condition = f"#script ran in {timing} for file {name}. {perfect_counter+imperfect_counter} reads out of {reads} were considered valid. {perfect_counter} were perfectly aligned. {imperfect_counter} were aligned with mismatch"
    
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
        
    with open(path, "w", newline='') as output: #writes the output
        writer = csv.writer(output)
        writer.writerows(outfile)

def inputs_handler():
    
    """ assertains the correct parsing of the input parameters"""
    
    parameters=inputs_initializer()

    try:
        parameters["start"]=int(parameters["start"])
        parameters["length"]=int(parameters["length"])
        parameters["miss"]=int(parameters["miss"])
        parameters["phred"]=int(parameters["phred"])
        parameters["miss_search"]=int(parameters["miss_search"])
    except Exception:
        input("\nPlease confirm you have provided the correct parameters.\nOnly numeric values are accepted in the folowing fields:\n-Feature read starting place;\n-Feature length;\n-mismatch;\n-Phred score;\n-mismatches in the search sequence.\n\nPlease try again. Press enter to exit")
        raise Exception    
    
    # avoids getting -1 and actually filtering by highest phred score by mistake
    if int(parameters["phred"]) == 0:
        parameters["phred"] = 1
    
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
        if len(parameters) != 13:
            input("Please confirm that all the input boxes are filled. Some parameters are missing.\nPress enter to exit")
            raise Exception
            
    parameters["cmd"] = False
    parameters['cpu'] = False

    return parameters

def inputs_initializer():
    
    """ Handles the graphical interface, and all the parameter inputs"""
    
    from tkinter import Entry,LabelFrame,Button,Label,Tk,filedialog,StringVar,OptionMenu
    separator = "\\" if system() == 'Windows' else "/"
    
    def restart():
        root.quit()
        root.destroy()
        inputs_initializer()
        
    def submit():
        for arg in temporary:
            parameters[arg] = temporary[arg].get()
        root.quit()
        root.destroy()
    
    def directory(column,row,parameter,frame):
        filename = filedialog.askdirectory(initialdir = separator, title = "Select a folder")
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

    def browsing(keyword,inputs,placeholder="Use the Browse button to navigate"):
        title1,title2,row,column,function=inputs
        frame=LabelFrame(root,text=title1,padx=5,pady=5)
        frame.grid(row=row,column=column, columnspan=2)
        button = Button(frame, text = title2,command = lambda: function(column,row,keyword,frame))
        button.grid(column = column, row = row)
        button_place = Entry(frame,borderwidth=5,width=50)
        button_place.grid(row=row,column=column+1)
        button_place.insert(0, placeholder)
        
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
    root.title("2FAST2Q Input Parameters Window")
    root.minsize(425, 600)
    parameters,temporary = {},{}  

    browsing_inputs = {"seq_files":["Path to the .fastq(.gz) files folder","Browse",1,0,directory],
                       "feature":["Path to the features .csv file","Browse",2,0,file],
                       "out":["Path to the output folder","Browse",3,0,directory]}

    default_inputs = {"start":["Feature start position in the read",6,0,0],
                      "length":["Feature length",7,0,20],
                      "miss":["Allowed mismatches",8,0,1],
                      "phred":["Minimal feature Phred-score",9,0,30],
                      "delete":["Delete intermediary files [y/n]",11,0,"y"],
                      "upstream":["upstream search sequence",12,0,"None"],
                      "downstream":["downstream search sequence",13,0,"None"],
                      "miss_search":["mismatches in the search sequence",14,0,0],}
    
    dropdown_options = {"Running Mode":["Counter", "Extractor + Counter",4,0]}

    # Generating the dropdown browsing buttons
    [dropdown(arg,dropdown_options[arg]) for arg in dropdown_options]
    
    # Generating the file/folder browsing buttons
    [browsing(arg,browsing_inputs[arg]) for arg in browsing_inputs]
    
    # Generating the input parameter buttons
    [write_menu(arg,default_inputs[arg]) for arg in default_inputs]
    
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
 
    version = "2.4"
    
    print("\nVersion: {}".format(version))

    param = inputs_handler() if cmd is None else cmd

    param["version"] = version
    
    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI" #Phred score

    param["quality_set"] = set(quality_list[:int(param['phred'])-1])
    
    current_time = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    param["directory"] = os.path.join(param['out'], f"2FAST2Q_output_{current_time}")
    
    if not os.path.exists(param["directory"]):
        os.makedirs(param["directory"])
    
    if psutil.virtual_memory().percent>=75:
        print("\nLow RAM availability detected, file processing may be slow\n")
    
    if param['Running Mode']=='C':
        print(f"\nRunning in align and count mode with the following parameters:\n{param['miss']} mismatch allowed\nMinimal Phred Score per bp >= {param['phred']}\nFeature length: {param['length']}\nRead alignment start position: {param['start']}\n")
    else:
        print(f"\nRunning in extract and count mode with the following parameters:\nMinimal Phred Score per bp >= {param['phred']}\n")
        
        if param['upstream'] is not None:
            print(f"Upstream search sequence: {param['upstream']}\n")
        if param['downstream'] is not None:
            print(f"Downstream search sequence: {param['downstream']}\n")
        
        if (param['upstream'] is not None) and (param['downstream'] is not None):
            print(f"Mismatches allowed in the search sequence: {param['miss_search']}\n")
            
        elif (param['upstream'] is not None) or (param['downstream'] is not None):
            print(f"Mismatches allowed in the search sequence: {param['miss_search']}\nReturned feature lenght: {param['length']}\n")

    print(f"All data will be saved into {param['directory']}")

    return param

def input_parser():
    
    """ Handles the cmd line interface, and all the parameter inputs"""
    
    def current_dir_path_handling(param):
        if param[0] is None:
            parameters[param[1]]=os.getcwd()
            if param[1] == 'feature':
                file = path_finder(os.getcwd(), ["*.csv"])
                if len(file) > 1:
                    input("There is more than one .csv in the current directory. If not directly indicating a path for sgRNA.csv, please only have 1 .csv file.") 
                    raise Exception
                if len(file) == 1:
                    parameters[param[1]]=file[0][0]
        else:
            parameters[param[1]]=param[0]
        return parameters
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c",nargs='?',const=True,help="cmd line mode")
    parser.add_argument("--s",help="The full path to the directory with the sequencing files OR file")
    parser.add_argument("--g",help="The full path to the .csv file with the sgRNAs.")
    parser.add_argument("--o",help="The full path to the output directory")
    parser.add_argument("--m",help="number of allowed mismatches (default=1)")
    parser.add_argument("--ph",help="Minimal Phred-score (default=30)")
    parser.add_argument("--st",help="Feauture start position in the read (default is 0==1st bp)")
    parser.add_argument("--l",help="Feature length (default=20bp)")
    parser.add_argument("--us",help="Upstream search sequence")
    parser.add_argument("--ds",help="Downstream search sequence")
    parser.add_argument("--ms",help="mismatches allowed when searching reads with Up/Down stream sequences")
    parser.add_argument("--mo",help="Running Mode (default=C) [Counter (C) / Extractor + Counter (EC)]")
    parser.add_argument("--cp",help="Number of cpus to be used (default is max(cpu)-2 for >=3 cpus, -1 for >=2 cpus, 1 if 1 cpu")
    parser.add_argument("--k",nargs='?',const=False,help="If enabled, keeps all temporary files (default is disabled)")
    args = parser.parse_args()
    
    #if its not running on command window mode
    if args.c is None:
        return None
    
    parameters = {}
    parameters["cmd"] = True
    paths_param = [[args.s,'seq_files'],
                   [args.g,'feature'],
                   [args.o,'out']]
    
    for param in paths_param:
        parameters = current_dir_path_handling(param)

    parameters['length']=20
    if args.l is not None:
        parameters['length']=int(args.l)
                 
    parameters['start']=0
    if args.st is not None:
        parameters['start']=int(args.st)
        
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
                 
    parameters['miss_search']=0
    if args.ms is not None:
        parameters['miss_search']=int(args.ms)
        
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
            [f"#Mismatches in search sequence: {param['miss_search']}"]
    headers = headers[::-1]

    compiled = {} #dictionary with all the reads per feature
    head = ["#Feature"] #name of the samples
    for i, file in enumerate(ordered_csv):
        path,_ = os.path.splitext(file)
        path = path[:-len("_reads")]
        head.append(path)
        with open(file) as current:
            for line in current:
                line = line[:-1].split(",")
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

    run_stats(headers,param)
                        
    final = []
    [final.append([feature] + compiled[feature]) for feature in compiled] 
    final.insert(0, head)
    
    csvfile = os.path.join(param["directory"],"compiled.csv")
    csv_writer(csvfile, final)

    if param["delete"]:
        for file in ordered_csv:
            os.remove(file)

    if not param["cmd"]:
        input("\nAnalysis successfully completed\nAll the reads have been compiled into the compiled.csv file.\nPress enter to exit")

def run_stats(headers, param):
    
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
            
    csvfile = os.path.join(param["directory"],"compiled_stats.csv")
    csv_writer(csvfile, global_stat)
    
    ### plotting
    header_ofset = 10
    fig, ax = plt.subplots()
    width = 0.4
    for i, (a,b,c,d,e,f,g) in enumerate(global_stat[header_ofset:]):   

        plt.bar(i+width/2, int(d), width,  capsize=5, color = "#FFD25A", hatch="//",edgecolor = "black", linewidth = .7)
        plt.bar(i-width/2, int(e), width, capsize=5, color = "#FFAA5A", hatch="\\\\\\",edgecolor = "black", linewidth = .7)

    ax.set_xticks(np.arange(len([n[0] for n in global_stat[header_ofset:]])))
    ax.set_xticklabels([n[0] for n in global_stat[header_ofset:]])
    plt.ylabel('Number of reads')
    plt.xticks(rotation=90)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(["Reads in sample", "Reads that passed quality control parameters"], \
              loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=1,prop={'size': 8})
    plt.tight_layout()
    file = os.path.join(param["directory"],"reads_plot.png")
    plt.savefig(file, dpi=300, bbox_inches='tight')

def ram_lock():

    """ stops increasing the hash tables (failed and passed reads) 
    size to avoid using RAM that isnt there """

    if psutil.virtual_memory().percent >= 95:
        return False
    return True

def cpu_counter(cpu):
    
    """ counts the available cpu cores, required for spliting the processing
    of the files """
    
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
        return hash_reads_parsing(result,failed_reads_compiled,passed_reads_compiled,failed_reads,passed_reads)
    
    else:
        return failed_reads,passed_reads

def hash_reads_parsing(result,failed_reads_compiled,passed_reads_compiled,failed_reads,passed_reads):
    
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
    
    print("\nInitializing hash tables.")
    result=[]
    for name in files[:cpu]:
        result.append(pool.apply_async(reads_counter, args=(0,0,name,features,param,cpu,set(),{},True)))
    
    compiled = [x.get() for x in result]
    throw1,throw2,throw3,throw4,failed_reads_compiled,passed_reads_compiled = zip(*compiled)
    
    return hash_reads_parsing(result,failed_reads_compiled,passed_reads_compiled,set(),{})

def aligner_mp_dispenser(files,features,param):
    
    """ starts and handles the parallel processing of all the samples by calling 
    multiple instances of the "aligner" function (one per sample) """
    
    start,failed_reads,passed_reads = 0,set(),{}
    pool,cpu = cpu_counter(param["cpu"])
    
    if (len(files)>cpu) & (param["miss"] != 0):
        failed_reads,passed_reads=hash_preprocesser(files,features,param,pool,cpu)
    
    print(f"\nProcessing {len(files)} files. Please hold.")
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
