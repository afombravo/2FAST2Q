import csv
import time
import glob, os
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor
import multiprocessing 
import platform
import regex

#### multiprocessing doesnt run on spyder. use terminal for debugging

#####################
    
class SgRNA:
    def __init__(self, name, sequence, counts):
        
        self.name = name
        self.sequence = sequence
        self.counts = counts
        
        
class Read:
    def __init__(self, name, sequence, quality):
        
        self.name = name
        self.sequence = sequence
        self.quality = quality

def imperfect_find(seq,guide,diffnumber):

    if regex.match("(%s" % seq + "){s<=%s" % diffnumber + "}", guide):
        return True
            
def path_finder_seq(folder_path, extenction, separator): 

    pathing = []
    for filename in glob.glob(os.path.join(folder_path, extenction)):
        stop = filename[::-1].find(separator) + 1 # / for linux, use \\ for windows
        pathing.append([filename[-stop:]] + [filename] + [os.path.getsize(filename)]) 

    if extenction == '*.fastq.gz':
        #sorting by size makes multiprocessing more efficient 
        #as the small files will be ran first
        
        ordered = [path[:2] for path in sorted(pathing, key=lambda e: e[-1], reverse = False)]
     
        if ordered == []:
            input("Check the path to the fastq.gz files folder. No files found.\n Press any key to exit")

    else:
        ordered = [path[:2] for path in sorted(pathing, reverse = False)]

    return ordered

def unzip(file,write_path):
    
    f = gzip.open(file, 'rb')
    
    if not os.path.isfile(write_path): 
        
        with open(write_path, 'wb') as f_out:
            shutil.copyfileobj(f, f_out)
                
    return 

def unpack(ordered,directory):
    
    write_path_save,filing = [],[]
    
    for name, filename in ordered: 
        filing.append(filename)
        write_path_save.append(directory + name[:-len(".gz")])
        
    with ThreadPoolExecutor() as executor: # multhreaded
        executor.map(unzip,filing,write_path_save)

    return write_path_save


def guides_loader(guides):
    
    if not os.path.isfile(guides):
        input("\nCheck the path to the sgRNA file.\nNo file found in the folowing path: {}\nPress any key to exit".format(guides))
    
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

def reads_loader(raw,quality_set, start, lenght):
    
    reads = {}
    n = set("N")
    reading = []
    loading = 0
    guide_len = start + lenght
    
    with open(raw) as current:
        for line in current:
            reading.append(line[:-1])
            
            if len(reading) == 4: #a read always has 4 lines
                
                name = reading[0]
                seq = reading[1][start:guide_len].upper()
                quality = reading[3][start:guide_len]
                
                reading = []
                loading += 1
                
                if (len(quality_set.intersection(quality)) == 0) & \
                    (len(n.intersection(seq)) == 0):
                        
                    reads[name] = Read(name, seq, quality)
                    
    return reads, loading
  
    
def perfect_alignment(reads,sgrna):

    for read in list(reads):
                
        if reads[read].sequence in sgrna.keys():
            sgrna[reads[read].sequence].counts += 1
            del reads[read]

    return sgrna,reads


def imperfect_alignment(reads,sgrna,missmatch,counter):

    for read in reads:
                
        found = 0
        
        for guide in sgrna.keys():
            
            finder = imperfect_find(reads[read].sequence, guide, missmatch)

            if finder:

                found += 1
                found_guide = guide
                
                if found >=2:
                    break
                
        if found == 1: #only appends sequences that align to the gRNA set 
            sgrna[found_guide].counts += 1
            counter += 1

    return sgrna, counter


def aligner(raw,guides,quality_set,missmatch,i,o,sgrna,version,separator, start, lenght):

    tempo = time.time()
       
    print("Processing file {} out of {}".format(i+1,o))

    reads, unloading = reads_loader(raw,quality_set, start, lenght)

    sgrna, reads = perfect_alignment(reads,sgrna)
    
    counter = 0
    if missmatch != 0:

        sgrna, counter = imperfect_alignment(reads,sgrna,missmatch,counter)

    master_list = [["#sgRNA"] + ["Reads"]]
    total_reads = 0
    
    for guide in sgrna:
        master_list.append([sgrna[guide].name] + [sgrna[guide].counts])
        total_reads += sgrna[guide].counts

    master_list.sort(key = lambda master_list: master_list[0]) #alphabetical sorting
    
    tempo = time.time() - tempo

    if tempo > 60:
        timing = str(round(tempo / 60, 2)) + " minutes"
        
    else:
        timing = str(round(tempo, 2)) + " seconds"

    name = raw[-raw[::-1].find(separator):-len(".fastq")]

    master_list.insert(0,["#script ran in {} for file {}. {} reads out of {} were considered valid. {} were perfectly aligned. {} were aligned with missmatch".format(timing,name,total_reads,unloading,total_reads-counter,counter)])
    
    csvfile = raw[:-raw[::-1].find(".")-1] + "_reads.csv"
                          
    csv_writer(csvfile, master_list)
        
    print("\nscript ran in {} for file {}. {} reads out of {} were considered valid. {} were perfectly aligned. {} were aligned with missmatch\n".format(timing,name,total_reads,unloading,total_reads-counter,counter)) # quality control

    return


def csv_writer(path, outfile):
    
    with open(path, "w", newline='') as output: #writes the output
        writer = csv.writer(output)
        writer.writerows(outfile)
        
    return

def input_getter():
    
    inputs = os.path.join(os.getcwd(), "inputs.txt")
    
    input_parameters = []
    with open(inputs) as current:
        for line in current:
            if ("#" not in line[0]) and ("\n" not in line[0]):
                line = line[:-1].split('"')
                input_parameters.append(line[1])
                
    return input_parameters

def initializer():
    
    version = "1.3.0"
    
    print("\nVersion: {}".format(version))
    
    folder_path, guides, out, missmatch, phred, start, lenght = input_getter()
    
    quality_list = '!"#$%&' + "'()*+,-/0123456789:;<=>?@ABCDEFGHI" #Phred score

    quality_set = set(quality_list[:int(phred)-1])
    
    directory = os.path.join(out, "unpacked")
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    print("\nRunning with parameters:\n{} Missmatch allowed\nMinimal Phred Score per bp >= {}\n".format(missmatch, phred))
    print("All data will be saved into {}".format(directory))
    
    if platform.system() == 'Windows':
        separator = "\\"
    else:
        separator = "/"
    
    return folder_path, guides, int(missmatch), quality_set, directory, \
        version, int(phred), separator, int(start), int(lenght)


def compiling(directory,phred,missmatch,version,separator):
    
    ordered_csv = path_finder_seq(directory, '*reads.csv',separator)
    
    headers = ["Parameters:\n\nVersion: {}\n".format(version)] + \
            ["Missmatch: {}\n".format(missmatch)] + \
            ["Phred Score: {}\n\n".format(phred)]
            
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

    return


def multi(files,guides,quality_set,missmatch,sgrna,version,separator, start, lenght):

    cpu = multiprocessing.cpu_count()
    
    if cpu >= 2:
        cpu -= 1
        
    pool = multiprocessing.Pool(processes = cpu)
    
    for i,name in enumerate(files):
        pool.apply_async(aligner, args=(name,guides,quality_set,missmatch,i,len(files),sgrna,version,separator, start, lenght))
        
    pool.close()
    pool.join()
        
    return


def main():
    
    folder_path, guides,missmatch, quality_set,directory, \
    version,phred,separator,start, lenght = initializer()

    print("\nUnpacking .gz files")
    
    ordered = path_finder_seq(folder_path, '*.fastq.gz',separator)
    files = unpack(ordered,directory)
    
    print("\nAligning sgRNAs")
    
    sgrna = guides_loader(guides)
    
    multi(files,guides,quality_set,missmatch,sgrna,version,separator, start, lenght)
    
    compiling(directory,phred,missmatch,version,separator)
    
    return


if __name__ == "__main__":
    multiprocessing.freeze_support() # required to run multiprocess as exe on windows
    main()