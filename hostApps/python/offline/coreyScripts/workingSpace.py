#This file is to mess around with things that I don't want to include in the tcr plotter file
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.stats import poisson
from scipy.optimize import curve_fit
from matplotlib import colors
from tcrPlotterV2 import tcrPlotterMultipleRuns
from scipy.special import gamma

dir_path = "../../../../../data"
myFiles = os.listdir(dir_path)

dflist = []

for file in myFiles:
    filePath = os.path.join(dir_path, file)
    if os.path.isfile(filePath): #check that its a file
        dflist.append(pd.read_csv(filePath, sep=';'))


'''pdcAverageTcr = [[] for i in range(len(dflist))]
for i in range(len(dflist)): #loop over runs
    for j in range(4): #loop over pdcs
        count = 0
        pdcAverageTcr[i].append(dflist[i][f'SPAD_TCR{j}'].sum() / len(dflist[i][f'SPAD_TCR{j}']))

        for k in dflist[i][f'SPAD_TCR{j}']:
            if k >= pdcAverageTcr[i][j]:
                dflist[i][f'SPAD_TCR{j}'] = dflist[i][f'SPAD_TCR{j}'].replace(k, np.nan)
                count += 1
        print(f'file {i}, PDC{j}: SPADs removed = {count}. percentage = {count * 100 / 64}')'''

ipdc = 0
ifile = 3
differenceList = []
meanList = []
medianList = []
maxIndexList = []
tcrList = []

differenceList.append(dflist[ifile][f'SPAD_TCR{ipdc}'].mean() - dflist[ifile][f'SPAD_TCR{ipdc}'].median())
meanList.append(dflist[ifile][f'SPAD_TCR{ipdc}'].mean())
medianList.append(dflist[ifile][f'SPAD_TCR{ipdc}'].median())

count = 1
for k in range(len(dflist[ifile][f'SPAD_TCR{ipdc}'])):
    #dflist[ifile][f'SPAD_TCR{ipdc}'] = dflist[ifile][f'SPAD_TCR{ipdc}'].replace(dflist[ifile][f'SPAD_TCR{ipdc}'].max(), np.nan)
    maxIndexList.append(dflist[ifile][f'SPAD_TCR{ipdc}'].idxmax())
    tcrList.append(dflist[ifile][f'SPAD_TCR{ipdc}'].max())
    dflist[ifile].loc[dflist[ifile][f'SPAD_TCR{ipdc}'].idxmax(), f'SPAD_TCR{ipdc}'] = np.nan
    count += 1
    if count > 63:
        break
    differenceList.append(dflist[ifile][f'SPAD_TCR{ipdc}'].mean() - dflist[ifile][f'SPAD_TCR{ipdc}'].median())
    meanList.append(dflist[ifile][f'SPAD_TCR{ipdc}'].mean())
    medianList.append(dflist[ifile][f'SPAD_TCR{ipdc}'].median())
    if differenceList[count-1] < 0:
        break


'''plt.plot([i for i in range(len(differenceList))], differenceList)
plt.show()'''
p = poisson.pmf([i for i in range(600)], 300)
plt.figure(0)
plt.plot([i for i in range(600)], p)
def model_func(x,a):
                #y = [i+1 for i in x]
                return (a ** x) * math.exp(-a) / gamma(x+1)
lst=[]
for i in [i for i in range(600)]:
      lst.append(model_func(i,300))
plt.figure(1)
plt.plot([i for i in range(600)], lst)

plt.show()
print(sum(p))
