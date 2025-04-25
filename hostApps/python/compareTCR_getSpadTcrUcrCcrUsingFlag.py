import pandas as pd
import numpy as np
import sys, glob, os
import matplotlib.pyplot as plt

nmbPDCs = 4
printvals = False

try:
    fIn = sorted(glob.glob(sys.argv[1]+"/*.csv"))
    #print(fIn)
    dfs_In = [pd.read_csv(i, sep=';') for i in fIn]
    df = pd.concat(dfs_In, axis=1)
except IndexError:
    print("Input a directory with output *csv files ")
    sys.exit()

#Add cols for ratios
for i in range(nmbPDCs):
    df[f'PDC{i}_UCR/TCR (%)'] = df[f'PDC{i}_UCR (cps)']/df[f'PDC{i}_TCR (cps)']*100.
    df[f'PDC{i}_dead'] = df[f'PDC{i}_TOT'].where(df[f'PDC{i}_TOT']==0)
    df[f'PDC{i}_noUCRCCR'] = df[f'PDC{i}_UCR (cps)'].where(df[f'PDC{i}_UCR (cps)']==-1)

if printvals:
    #Print number dead and problematic ZPP fit pixels
    print(df.count())

    # Print average and median values
    print("DF Mean values")
    for i in range(4):
        print(f"PDC{i}")
        print(df[f'PDC{i}_CCR (%)'][df[f'PDC{i}_CCR (%)']>0].mean())
        print(df[f'PDC{i}_UCR/TCR (%)'][df[f'PDC{i}_UCR/TCR (%)']>0].mean())

    print("DF median values")
    for i in range(4):
        print(f"PDC{i}")
        print(df[f'PDC{i}_CCR (%)'][df[f'PDC{i}_CCR (%)']>0].median())
        print(df[f'PDC{i}_UCR/TCR (%)'][df[f'PDC{i}_UCR/TCR (%)']>0].median())


#pull in TCR only measurement to compare the two scripts
pathIn = '/'.join(sys.argv[1].split('/')[:-1])+"/" #remove last folder in path

if os.path.exists(pathIn+"getSpadTcrUsingFlag"):
    pathIn+="getSpadTcrUsingFlag/"
elif os.path.exists(pathIn+"getSpadTcrUsingFlag_AS"):
    pathIn+="getSpadTcrUsingFlag_AS/"
else:
    print("Cannot find folder with TCR-only run")
    sys.exit()

try:
    fin = glob.glob(pathIn+'*.csv')
    df_tcr = pd.read_csv(fin[0], sep=';') 
except FileNotFoundError:
    print("TCR file not found")
    sys.exit()

#create combo df
measTime_tcr = 0.2
combo_df = pd.DataFrame()
combo_df["tuc_0"]=df['PDC0_TCR (cps)']
combo_df["tuc_1"]=df['PDC1_TCR (cps)']
combo_df["tuc_2"]=df['PDC2_TCR (cps)']
combo_df["tuc_3"]=df['PDC3_TCR (cps)']
combo_df['t_0']=df_tcr['SPAD_TCR0']/measTime_tcr
combo_df['t_1']=df_tcr['SPAD_TCR1']/measTime_tcr
combo_df['t_2']=df_tcr['SPAD_TCR2']/measTime_tcr
combo_df['t_3']=df_tcr['SPAD_TCR3']/measTime_tcr
combo_df['t/tuc_0']=combo_df['t_0']/combo_df['tuc_0']
combo_df['t/tuc_1']=combo_df['t_0']/combo_df['tuc_1']
combo_df['t/tuc_2']=combo_df['t_0']/combo_df['tuc_2']
combo_df['t/tuc_3']=combo_df['t_0']/combo_df['tuc_3']

fig = plt.figure()
plt.plot([combo_df[f't/tuc_{i}'] for i in range(4)], 'o-',label=[str(i) for i in range(64)])
plt.legend(loc='best')
plt.yscale('log')
plt.show()