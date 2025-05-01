import pandas as pd
import numpy as np
import sys, glob, os
import matplotlib.pyplot as plt


def plotScans(dict, vals, save, name='scan'):
    for pdc in range(4):
        fig=plt.figure()
        plt.plot(np.array([list(dict.values())[i][f'SPAD_TCR{pdc}'] for i in range(len(dict.values()))]).T, label=vals)
        plt.legend(loc='best')
        plt.yscale('log')
        plt.xlabel('SPAD index')
        plt.ylabel('TCR (cps)')
        plt.title(f'{name} PDC {pdc}')
        if save:
            plt.savefig(f'plots/timing_{name}_TCR_pdc{pdc}.png')
        else:
            plt.show()  

    


printvals = False
savePlots = True

try:
    fIn = sorted(glob.glob(sys.argv[1]+"/*.csv"))
    dfs_In = [pd.read_csv(i, sep=';') for i in fIn]
except IndexError:
    print("Input a directory with output *csv files ")
    sys.exit()

#duplicate and rename nominal to explicitly include parameter values
fIn.append(fIn[-1])
fIn.append(fIn[-1])
dfs_In.append(dfs_In[-1].copy())
dfs_In.append(dfs_In[-1].copy())
fIn[-1] = fIn[-1][:-4]+"_flag10.csv" #rename
fIn[-2] = fIn[-2][:-4]+"_rech4.csv" #rename
fIn[-3] = fIn[-3][:-4]+"_holdoff150.csv" #rename

#Group list into categories of what variable is changing
var = [n.split('_')[-1][:-4] for n in fIn]
dict_dfs = dict(zip(var,dfs_In))

#Isolate each variable
var_flag = [i for i in var if 'flag' in i]
dict_flag = {k: dict_dfs[k] for k in var_flag}
var_rech = [i for i in var if 'rech' in i]
dict_rech = {k: dict_dfs[k] for k in var_rech}
var_holdoff = [i for i in var if 'holdoff' in i]
dict_holdoff = {k: dict_dfs[k] for k in var_holdoff}


#plot
plotScans(dict_flag, var_flag, savePlots, name="flag")
plotScans(dict_rech, var_rech, savePlots, name="rech")
plotScans(dict_holdoff, var_holdoff, savePlots, name="holdoff")
