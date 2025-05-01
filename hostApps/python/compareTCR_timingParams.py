import pandas as pd
import numpy as np
import sys, glob, os
import matplotlib.pyplot as plt
import re


def plotScans(dict, vals, save, name='scan',show=True):
    medians=[]
    for pdc in range(4):
        #transpose data
        dataToPlot = np.array([list(dict.values())[i][f'SPAD_TCR{pdc}'] for i in range(len(dict.values()))]).T
        fig=plt.figure()
        plt.plot(dataToPlot, label=vals)
        tmp_medians = np.median(dataToPlot, axis=0)
        medians.append(tmp_medians)
        for i,m in enumerate(tmp_medians):
            plt.axhline(y=m, linestyle='--',color=f"C{i}", label=f"{m:.2f}")
        plt.yscale('log')
        plt.xlabel('SPAD index')
        plt.ylabel('TCR (cps)')
        plt.title(f'{name} PDC {pdc}')
        plt.legend(loc='best',ncol=2)
        if save:
            plt.savefig(f'plots/timing_{name}_TCR_pdc{pdc}.png')
        elif show:
            plt.show()
        else:
            pass 
        plt.close()
    return medians
    
def plotMedians(meds, vals_full,save):
    #get values out of vars and turn them into ints
    vals_ind=[re.search(r"\d", v).start() for v in vals_full] #overkill - could just compute for the first entry and then apply to all because they all have the same leading string
    vals = [int(v[i:]) for v,i in zip(vals_full,vals_ind)]
    name = vals_full[0][:vals_ind[0]]
    
    if len(meds[0])!=len(vals):
        print("Error in plotting medians")
        sys.exit()

    fig=plt.figure()
    for i,m in enumerate(meds):
        #plot sorted arrays
        plt.plot(sorted(vals),[x for _, x in sorted(zip(vals, m))], 'o-',label=f"PDC{i}")
    plt.legend(loc='best')
    plt.title(name)
    plt.xlabel('Scan value')
    plt.ylabel('Median TCR')
    if save:
        plt.savefig(f'plots/timing_{name}_TCR__medians.png')
    else:
        plt.show()
    plt.close()








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
medians_flag = plotScans(dict_flag, var_flag, savePlots, name="flag", show=False)
medians_rech = plotScans(dict_rech, var_rech, savePlots, name="rech", show=False)
medians_holdoff = plotScans(dict_holdoff, var_holdoff, savePlots, name="holdoff", show=False)


plotMedians(medians_flag, var_flag, savePlots)
plotMedians(medians_rech, var_rech, savePlots)
plotMedians(medians_holdoff, var_holdoff, savePlots)