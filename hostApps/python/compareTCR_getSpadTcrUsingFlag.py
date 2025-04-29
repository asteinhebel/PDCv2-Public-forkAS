import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

nmbPDCs = 4
savePlots = True
try:
    if nmbPDCs<4:
        df = [pd.read_csv(i) for i in sys.argv[1].split(',')]
    else:
        df = [pd.read_csv(i, sep=';') for i in sys.argv[1].split(',')]
except IsADirectoryError:
    print("Input an array of file names as an argument, like F1,F2,F3")
    sys.exit()

#Make label array of PDC number
pdcName = [i.split('/')[5] for i in sys.argv[1].split(',')]

#plot individual curves
for i,d in enumerate(df):
    plt.figure()
    for j in range(nmbPDCs):
        plt.plot(d[f'SPAD_percent{j}'], d[f'SPAD_distribution{j}'], label=j)
    plt.legend(loc='best')
    plt.xlabel('Percent [%]')
    plt.yscale('log')
    plt.title(pdcName[i])
    if savePlots:
        plt.savefig(f"plots/compareTCR_turnonPerc_{pdcName[i]}.png")
    else:
        plt.show()
    plt.close()

#calculate 2nd derivative
for i,d in enumerate(df):
    for j in range(nmbPDCs):
        d[f'SPAD_dist_diff{j}'] = d[f'SPAD_distribution{j}'].diff()
        d[f'SPAD_dist_diffdiff{j}'] = d[f'SPAD_dist_diff{j}'].diff()
        d[f'SPAD_turnonPerc{j}'] = d[f'SPAD_percent{j}'].iloc[np.argwhere(d[f'SPAD_dist_diffdiff{j}']>25)[0]]

#plot all curves together
plt.figure()
for i,d in enumerate(df):
    for j in range(nmbPDCs):
        k=(i*3)+j
        plt.plot(d[f'SPAD_percent{j}'], d[f'SPAD_distribution{j}'], label=pdcName[i]+" "+str(j),color=f"C{k}")
        plt.scatter(d[f'SPAD_turnonPerc{j}'], d[f'SPAD_distribution{j}'],s=50, color=f"C{k}")
plt.legend(loc='best')
plt.yscale('log')
plt.xlabel('Percent [%]')
if savePlots:
    plt.savefig("plots/compareTCR_turnonPerc_all.png")
else:
    plt.show()
plt.close()



