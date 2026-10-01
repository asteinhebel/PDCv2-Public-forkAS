import sys
import os
import datetime
import pandas as pd
import numpy as np

# -----------------------------------------------
# --- methods
# -----------------------------------------------

def make_file_outname(ins):
    return datetime.datetime.now(tz=datetime.UTC).strftime("%Y%m%d-%Hh%Mm%S")

def checkPixMap(ins, fout):

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

def checkHoldoff(ins, fout):
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

def mergeCSV(ins, fout):    
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
def main(input_list):

    #define output name
    fout = make_file_outname(input_list)

    #confirm that same pixel map was used for all and create new output if so
    checkPixMap(input_list, fout)

    #confirm that same holdoff time was used for all and create new output if so. ONLY HOLDOFF VALUE IS VALID FOR MERGED FILE
    checkHoldoff(input_list, fout)

    #merge CSVs into one
    mergeCSV(input_list, fout)

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

    main(inputs_that_exist)