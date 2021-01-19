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

#####################

class SgRNA:
    def __init__(self, name, sequence, counts):
        
        self.name = name
        self.sequence = sequence
        self.counts = counts
        
def path_finder_seq(folder_path, extension, separator): 
    
    """ Finds the correct file paths """
    
    pathing = []
    for filename in glob.glob(os.path.join(folder_path, extension)):
        stop = filename[::-1].find(separator) + 1 # / for linux, use \\ for windows
        pathing.append([filename[-stop:]] + [filename] + [os.path.getsize(filename)]) 

    if extension != '*reads.csv':
        #sorting by size makes multiprocessing more efficient 
        #as the bigger files will be ran first
        
        ordered = [path[:2] for path in sorted(pathing, key=lambda e: e[-1])][::-1]
     
        if ordered == []:
            input(f"Check the path to the {extension[1:]} files folder. No files of this type found.\n Press any key to exit")

    else:
        ordered = [path[:2] for path in sorted(pathing, reverse = False)]

    return ordered

def unzip(file,write_path):
    
    f = gzip.open(file, 'rb')
    
    if not os.path.isfile(write_path): 
        
        with open(write_path, 'wb') as f_out:
            shutil.copyfileobj(f, f_out)
                

def unpack(ordered,directory):
    
    write_path_save,filing = [],[]
    
    for name, filename in ordered: 
        filing.append(filename)
        write_path_save.append(directory + name[:-len(".gz")])
    
    with ThreadPoolExecutor() as executor: # multhreaded
        executor.map(unzip,filing,write_path_save)

    return write_path_save


def guides_loader(guides):
    
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
    
    reads = {}
    n = set("N")
    reading = []
    perfect_counter, imperfect_counter = 0, 0
    guide_len = start + lenght
    
    with open(raw) as current:
        for line in current:
            reading.append(line[:-1])
            
            if len(reading) == 4: #a read always has 4 lines
                
                seq = reading[1][start:guide_len].upper()
                quality = reading[3][start:guide_len]
                
                reading = []
                perfect_counter += 1
                
                if (len(quality_set.intersection(quality)) == 0) & \
                    (len(n.intersection(seq)) == 0):
                        
                    if seq in sgrna:
                        sgrna[seq].counts += 1
                    
                    elif mismatch != 0:
                    
                        sgrna, imperfect_counter = imperfect_alignment(seq, sgrna, mismatch, imperfect_counter)
                    
    return reads, perfect_counter, imperfect_counter
  

def imperfect_find(seq,guide,diffnumber):

    if match("(%s" % seq + "){s<=%s" % diffnumber + "}", guide):
        return True
            

def imperfect_alignment(sequence, sgrna, mismatch, counter):
    
    found = 0
    for guide in sgrna:
        
        finder = imperfect_find(sequence, guide, mismatch)

        if finder:

            found += 1
            found_guide = guide
            
            if found >=2:
                break           

    if found == 1: #only appends sequences that align to the gRNA set 
        sgrna[found_guide].counts += 1
        counter += 1

    return sgrna, counter


def aligner(raw, guides, out, quality_set,mismatch,i,o,sgrna,version,separator, start, lenght):

    tempo = time()
       
    print("Processing file {} out of {}".format(i+1,o))

    reads, perfect_counter, imperfect_counter = reads_counter(raw, quality_set, start, lenght, sgrna, mismatch)

    master_list = [["#sgRNA"] + ["Reads"]]
    total_reads = 0
    
    for guide in sgrna:
        master_list.append([sgrna[guide].name] + [sgrna[guide].counts])
        total_reads += sgrna[guide].counts

    master_list.sort(key = lambda master_list: master_list[0]) #alphabetical sorting
    
    tempo = time() - tempo

    if tempo > 60:
        timing = str(round(tempo / 60, 2)) + " minutes"
        
    else:
        timing = str(round(tempo, 2)) + " seconds"

    name = raw[-raw[::-1].find(separator):-len(".fastq")]
    
    stats_condition = "#script ran in {} for file {}. {} reads out of {} were considered valid. {} were perfectly aligned. {} were aligned with mismatch".format(timing,name,total_reads,perfect_counter,total_reads-imperfect_counter,imperfect_counter)

    master_list.insert(0,[stats_condition])
    
    csvfile = out[:-out[::-1].find(".")-1] + "_reads.csv"
                          
    csv_writer(csvfile, master_list)
        
    print(stats_condition[1:]) # quality control



def csv_writer(path, outfile):
    
    with open(path, "w", newline='') as output: #writes the output
        writer = csv.writer(output)
        writer.writerows(outfile)

def inputs_asker(msgs, input_parameters, directr):
    
    print(msgs)
    ezi.msgbox(msg=msgs)
    
    if directr:
        input_parameters.append(ezi.diropenbox(msg = msgs))
    else:
        input_parameters.append(ezi.fileopenbox(msg=msgs))
        
    return input_parameters

def input_getter_graphical():

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
    
    version = "1.4"
    
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
    
    ordered_csv = path_finder_seq(directory, '*reads.csv',separator)
    
    headers = [f"Parameters:\n\nVersion: {version}\n"] + \
            [f"mismatch: {mismatch}\n"] + \
            [f"Phred Score: {phred}\n\n"]
            
    compiled = {}
    head = ["#sgRNA"]
    for name, file in ordered_csv:
        head.append(name[name.find(separator)+len(separator):name.find(".csv")])
        with open(file) as current:
            for line in current:
                line = line[:-1].split(",")
                if "#" not in line[0]:
                    if line[0] in compiled.keys():
                        compiled[line[0]] = compiled[line[0]] + [int(line[1])]
                    else:
                        compiled[line[0]] = [int(line[1])]
                        
                elif "#sgRNA" not in line[0]:
                    headers.append(line[0][1:]+"\n")
    
    path = ordered_csv[0][1]
    out_file = path[:-path[::-1].find(separator)-1]
    
    text_file = out_file + separator + "compiled_stats.txt"
    
    with open(text_file, "w") as text_file:
        for item in headers:
            text_file.write(item)
                        
    final = []
    for sgrna in compiled:
        final.append([sgrna] + compiled[sgrna])
        
    final.insert(0, head)
    
    csvfile = out_file + separator + "compiled.csv"
    
    csv_writer(csvfile, final)
        
    input("\nAnalysis successfully completed\nAll the reads have been compiled into the compiled.csv file.\nPress any key to exit")



def multi(files,guides, write_path_save, quality_set,mismatch,sgrna,version,separator, start, lenght):

    cpu = multiprocessing.cpu_count()
    
    if cpu >= 2:
        cpu -= 1
        
    pool = multiprocessing.Pool(processes = cpu)
    
    for i, (name, out) in enumerate(zip(files, write_path_save)):
        pool.apply_async(aligner, args=(name, guides, out, quality_set,mismatch,i,len(files),sgrna,version,separator, start, lenght))
        
    pool.close()
    pool.join()


def input_file_type(ordered, extension, directory):

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
    
    
    folder_path, guides,mismatch, quality_set,directory, \
    version,phred,separator,start, lenght, extension = initializer()

    ordered = path_finder_seq(folder_path, extension, separator)

    files, write_path_save = input_file_type(ordered, extension, directory)

    sgrna = guides_loader(guides)

    multi(files,guides, write_path_save, quality_set,mismatch,sgrna,version,separator, start, lenght)
    
    compiling(directory,phred,mismatch,version,separator)


if __name__ == "__main__":
    multiprocessing.freeze_support() # required to run multiprocess as exe on windows
    main()