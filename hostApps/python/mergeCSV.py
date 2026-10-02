import sys
import os
import datetime
import pandas as pd
import numpy as np

# -----------------------------------------------
# --- methods
# -----------------------------------------------

def make_file_outname():
    """ Define shared output name """

    return "mergedData/"+datetime.datetime.now(tz=datetime.UTC).strftime("%Y%m%d-%Hh%Mm%S")

def make_log():
    """Create log file with details about generation and which input files were used."""

    #get total measurement duration
    durations = [f.split('_')  for f in ins]
    durations_s = [i[:-1] for i in [i for r in durations for i in r] if i[-1]=="s"] #flatten durations and pull out entries that end in s after stripping the s
    #try to convert remaining value to floats to add. if it can't be converted, then it's not part of the timing 
    total_duration = 0
    for i in durations_s:
        try:
            total_duration+=float(i)
        except ValueError:
            continue

    with open(fout+".log", 'w') as f:
        f.write("INPUT FILES:\n")
        for inputs in ins:
            f.write(inputs+"\n")
        f.write("\n")
        f.write("TOTAL DURATION:\n")
        f.write(str(total_duration)+" s")


def check_pix_map():
    """Confirm that all pixel masks are the same. If so, save a new version wih the new dataset name"""

    enabled = [ [] for _ in range(len(ins)) ]

    for i,fin in enumerate(ins):
        try:
            with open(fin+'.txt', 'r') as f:
                for line in f:
                    if "Enabling" in line:
                        enabled[i].append(int(line.split(" ")[1]))
        except FileNotFoundError:
            pass

    #check that same map was used for all inputs
    ident_maps = [enabled.count(enabled[0])==len(enabled)]

    if not ident_maps:
        print("Must merge files that have identical pixel masks - these files do not.")
        sys.exit()

    #create new output with expected file name
    with open(ins[0]+'.txt','r') as firstfile, open(fout+'.txt','a') as secondfile:
        for line in firstfile:
                secondfile.write(line)

def check_holdoff():
    """Confirm that the same holdoff time was used for every input. If so, create new file with that time and new dataset name."""

    holdoff = []

    for i,fin in enumerate(ins):
        try:
            with open(fin+'_config.txt', 'r') as f:
                for line in f:
                    if "HOLD_TIME_NS" in line:
                        holdoff.append(float(line.split('=')[-1]))
        except FileNotFoundError:
            pass

    #check that same map was used for all inputs
    ident_holdoff = [holdoff.count(holdoff[0])==len(holdoff)]

    if not ident_holdoff:
        print("Must merge files that have identical holdoff times - these files do not.")
        sys.exit()

    #create new output with expected file name
    with open(ins[0]+'_config.txt','r') as firstfile, open(fout+'_config.txt','a') as secondfile:
        for line in firstfile:
                secondfile.write(line)

def merge_csv():  
    """Merge CSVs into one in chunks. Save to file of shared dataset name."""

    CHUNK_SIZE = 50000
    frameIdx_increment = 0
    last_largest_frameIdx = 0
    for i,csv_file_name in enumerate(ins):
        print(f"Merging {csv_file_name}")
        chunk_container = pd.read_csv(csv_file_name+'.csv', chunksize=CHUNK_SIZE, delimiter=';')
        firstChunk = True
        for chunk in chunk_container:
            #increment frameIdx so it looks like the merged file is one long, single acquisition
            if chunk['frameIdx'].min()+frameIdx_increment >= last_largest_frameIdx:
                last_largest_frameIdx = chunk['frameIdx'].max()+frameIdx_increment
            else:
                frameIdx_increment = last_largest_frameIdx+1
                last_largest_frameIdx += chunk['frameIdx'].max()
            chunk['frameIdx']+=frameIdx_increment
            #save output to new file. Only write header at the beginning
            header_bool = True if (i==0 and firstChunk) else False
            chunk.to_csv(fout+'.csv', mode="a", index=False, header=header_bool, sep=";")
            firstChunk=False

# -----------------------------------------------
# --- main
# -----------------------------------------------
def main():

    #define output name
    global fout 
    fout = make_file_outname()

    #confirm that same pixel map was used for all and create new output if so
    check_pix_map()

    #confirm that same holdoff time was used for all and create new output if so. ONLY HOLDOFF VALUE IS VALID FOR MERGED FILE
    check_holdoff()

    #merge CSVs into one
    merge_csv()

    #create note/log 
    make_log()

# -----------------------------------------------
# --- call to main
# -----------------------------------------------
if __name__ == "__main__":

    #check that inputs exist
    inputs = sys.argv[1:]
    inputs_that_exist = []
    for inp in inputs:
        if not os.path.exists(inp+".csv"):
            print(f"Could not find {inp} - skip this input")
        else:
            inputs_that_exist.append(inp)

    #check that output dir exists
    if not os.path.isdir('mergedData'):
        os.mkdir('mergedData/')

    global ins 
    ins=inputs_that_exist

    main()