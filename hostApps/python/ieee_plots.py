# run like ieee_plots.py <file path>

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------
# --- methods
# -----------------------------------------------

def check_path(path_in):
    """ Confirm path exists, return valid path """
    if not os.path.exists(path_in):
        print(f"Could not find {path_in} - skip this input")
        return
    else:
        return path_in


# -----------------------------------------------
# --- main
# -----------------------------------------------
def main(file_in):
    df = pd.read_csv(file_in, usecols=["energy","psd", "nPdcsInFrame", "pdcIdx"])

    df3 = df[df['pdcIdx']==3]
    df3_3or4 = df3[df3['nPdcsInFrame']>2]

    #PDC3 PSD
    fig, ax = plt.subplots()
    h = plt.hist2d(df3['energy'], df3['psd'], bins=[300, 100], range=[[0,3150],[0,1]]) #3150 = max of lead file
    fig.colorbar(h[3], ax=ax)
    print(len(df3))
    plt.show()
    plt.clf()

    #PDC3 PSD, only frames with 3 or 4 PDCs
    fig, ax = plt.subplots()
    h=plt.hist2d(df3_3or4['energy'], df3_3or4['psd'], bins=[300,100], range=[[0,3150],[0,1]])
    fig.colorbar(h[3], ax=ax)
    print(len(df3_3or4))
    plt.show()

    # prop of counts in top blob
    df3_3or4_top = df3_3or4[df3_3or4['psd']>0.65]
    print(len(df3_3or4_top)/len(df3_3or4))


# -----------------------------------------------
# --- call to main
# -----------------------------------------------
if __name__ == "__main__":

    #check that inputs exist
    input = check_path(sys.argv[1]) 

    main(input)