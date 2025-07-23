#This is meant to be a compilation of all the working plots that I have to analyze the PDCs
#Corey Fox
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
from matplotlib import colors

nPdcs = 4
nSpads = 64
colorlist = ['r', 'b', 'g', 'orange', 'y', 'k']

#Create dflist by taking files from a folder specified by dir_path. Not all functions utilize this dflist yet.
dir_path = "data"
myFiles = os.listdir(dir_path)

dflist = []

for file in myFiles:
    filePath = os.path.join(dir_path, file)
    if os.path.isfile(filePath): #check that its a file
        dflist.append(pd.read_csv(filePath, sep=';'))
        
#plots a single run across all PDCs on the same graph. Gives scatter plot of TCR and percent distribution
def tcrPlotter():
    df = pd.read_csv("20250602_13h37m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv" , sep=';')
    colors = ['r', 'b', 'g', 'y']

    for j in range(2):
        plt.figure(j + 1)
        for i in range(4):
            if j == 0:
                y_err = df[f'SPAD_TCR{i}'] ** 0.5
                plt.scatter(df[f'SPAD_idx{i}'], df[f'SPAD_TCR{i}'], facecolors = 'none', edgecolors = colors[i])
                plt.errorbar(df[f'SPAD_idx{i}'], df[f'SPAD_TCR{i}'], yerr = y_err, fmt = 'none', ecolor = colors[i])
                plt.plot(df[f'SPAD_idx{i}'], [df[f'SPAD_TCR{i}'].mean()] * 64, '--', color = colors[i])
                plt.title('TCR vs SPAD Index')
            else:
                plt.plot(df[f'SPAD_percent{i}'], df[f'SPAD_distribution{i}'], color = colors[i])
                plt.plot(list(range(101)), [df[f'SPAD_TCR{i}'].mean()] * 101, '--', color = colors[i])
            plt.yscale('log')
            plt.legend(['PDC1', 'Mean1', 'PDC2', 'Mean2', 'PDC3', 'Mean3', 'PDC4', 'Mean4'])
            
    plt.show()

#Plots a different graph for each small PDC
#Use this to analyze a PDC across runs, or to look at multiple big PDCs
def tcrPlotterMultipleRuns(dflist):
    plt.rcParams.update({'font.size': 12})
    '''dflist = [
        pd.read_csv("20250602_16h37m37_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep = ';'),
        pd.read_csv("20250602_16h32m37_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep = ';'),
        pd.read_csv("20250602_16h34m18_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep = ';'),
        pd.read_csv("20250602_16h36m05_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep = ';')
    ]'''
    fs = 16 #Fontsize for the plot
    colorlist = ['r', 'b', 'g', 'orange', 'y', 'k']
    nmbRuns = len(dflist)
    nmbPdcs = 4
    legend = []
    colors = []

    for i in range(nmbRuns):
        legend.append(f'Run{i + 1}')
        legend.append(f'Mean{i + 1}')
        colors.append(colorlist[i % len(colorlist)])

    for j in range(nmbPdcs): #j = pdc index
        for i in range(nmbRuns): #i = run index
            plt.figure(j)
            y_err = dflist[i][f'SPAD_TCR{j}'] ** 0.5
            plt.scatter(dflist[i][f'SPAD_idx{j}'], dflist[i][f'SPAD_TCR{j}'], facecolors = 'none', edgecolors=colors[i])
            #plt.errorbar(dflist[i][f'SPAD_idx{j}'], dflist[i][f'SPAD_TCR{j}'], yerr = y_err, fmt = 'none', ecolor = colors[i])
            #plt.plot(dflist[i][f'SPAD_idx{j}'], [dflist[i][f'SPAD_TCR{j}'].mean()]*64, '--', color=colors[i])
            plt.title(f'Tile{j} ({nmbRuns} Runs)')
            plt.yscale('log')

        plt.legend(legend, bbox_to_anchor=(1.02, 1.03))
        if j == 1:
            plt.legend(legend, bbox_to_anchor=(1.02, 0.9))
        plt.ylabel('TCR',fontsize=fs)
        plt.xlabel('SPAD Index',fontsize=fs)
        plt.tick_params(length = 10, width = 2)
        plt.tick_params(length = 10, width = 1, which = 'minor')

    plt.show()

#Gives a histogram of TCR using threshold method. One run at a time, each small PDC gets a figure.
def noScreamersThresholdHist():
    i=0
    for ipdc in range(nPdcs):
        dflist = [
        pd.read_csv( "20250602_13h37m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250602_16h32m37_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250603_16h47m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
        ]
        for x in dflist[i].index:
            if dflist[i].loc[x, f'SPAD_TCR{ipdc}'] > 500:
                dflist[i].drop(x,inplace = True)
        plt.figure(ipdc)
        plt.hist(dflist[i][f'SPAD_TCR{ipdc}'], bins = 'fd')

    plt.show()

#Plots histogram of before and after using threshold method. One run at a time, each small PDC gets a figure.
#Add quality factor when you get the chance
def noScreamersThresholdHistBnA():
    i=0 #Specify which file using i value
    threshold = 400
    print(f'Any SPAD with a TCR greater than {threshold} will be removed')
    for ipdc in range(nPdcs):
        print('-------------------------------------------------------')
        nScreamers = 0
        #Need to initialize dflist again here so that it takes the original dataframe at the start of each loop
        dflist = [
        pd.read_csv( "20250602_13h37m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250602_16h32m37_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250603_16h47m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
        ]

        plt.figure(ipdc)
        plt.subplot(1,2,1)
        plt.hist(dflist[i][f'SPAD_TCR{ipdc}'], bins = 'fd')
        plt.xscale('log')
        plt.title(f'Before, PDC{ipdc} (Log x-scale)')

        for x in dflist[i].index:
            if dflist[i].loc[x, f'SPAD_TCR{ipdc}'] > threshold:
                dflist[i].drop(x,inplace = True)
                nScreamers += 1

        plt.figure(ipdc)
        plt.subplot(1,2,2)
        plt.hist(dflist[i][f'SPAD_TCR{ipdc}'], bins = 'fd')
        plt.title(f'After, PDC{ipdc}')
        print(f'PDC{ipdc} has {nScreamers} Screamers\n{(1-nScreamers/nSpads)*100}% of SPADs remain')

    plt.show()

#Gives scatter plot of TCR without screamers using percent method. Plots multiple runs, figures separated by small PDC.
def noScreamersPerPlotMultipleRuns(dflist, index):
    nScreamers = []
    df = dflist
    percent = 100
    lst=[[] for ipdc in range(nPdcs)]
    dfSort = [[] for ipdc in range(nPdcs)]
    for ipdc in range(nPdcs):
        for ispad in range(nSpads):
            dfSort[ipdc].append(df[f'SPAD_TCR{ipdc}'][ispad])
            lst[ipdc].append(df[f'SPAD_TCR{ipdc}'][ispad])
        dfSort[ipdc].sort()

    for ipdc in range(nPdcs):
        n = 0 #number of screamers in a PDC, resets at each new PDC.
        i = 0
        while i < nSpads:
            if df[f'SPAD_percent{ipdc}'][i] >= percent:
                lst[ipdc].remove(dfSort[ipdc][i])
                n+=1
            i += 1
        nScreamers.append(n)

    for ipdc in range(nPdcs):
        plt.figure(ipdc)
        plt.scatter(range(0,nSpads - nScreamers[ipdc]), lst[ipdc], facecolors = 'none', edgecolors = colorlist[index])
        plt.plot(range(0,nSpads - nScreamers[ipdc]), [np.mean(lst[ipdc])] * (len(lst[ipdc])), '--', color = colorlist[index])
        plt.legend([1,1,2,2,3,3], bbox_to_anchor = (0.99,1))   
        plt.title(f'PDC{ipdc}')

#Used in conjunction with noScreamersPerPlot to get each plot
def getPlot(func):
    dflist = [
        pd.read_csv( "20250602_13h37m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250602_16h32m37_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250603_16h47m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
        ]
    
    for i in range(len(dflist)):
        func(dflist[i], i)
    plt.show()

#Gives scatter plot and percent distribution of a single run with every small PDC plotted at once using percent method.
def noScreamersPerPlot():
    df = pd.read_csv( "20250602_13h37m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
    nScreamers = []
    percent = 90
    lst=[[] for ipdc in range(nPdcs)] #Values used to plot
    dfSort = [[] for ipdc in range(nPdcs)] #list of sorted values to remove from lst
    for ipdc in range(nPdcs): #Get values for lst and dfsort
        for ispad in range(nSpads):
            dfSort[ipdc].append(df[f'SPAD_TCR{ipdc}'][ispad])
            lst[ipdc].append(df[f'SPAD_TCR{ipdc}'][ispad])
        dfSort[ipdc].sort()

    for ipdc in range(nPdcs):
        n = 0 #number of screamers in a PDC, resets at each new PDC.
        i = 0

        while i < nSpads: #remove SPAD with corresponding tcr if it is over the percent threshold
            if df[f'SPAD_percent{ipdc}'][i] >= percent:
                lst[ipdc].remove(dfSort[ipdc][i])
                n+=1
            i += 1
        nScreamers.append(n)

    for ipdc in range(nPdcs): #plots
        figure = 1
        while figure == 1:
            plt.figure(figure)
            plt.scatter(range(0,nSpads - nScreamers[ipdc]), lst[ipdc], facecolors = 'none', edgecolors = colorlist[ipdc])
            plt.plot(range(0,nSpads - nScreamers[ipdc]), [np.mean(lst[ipdc])] * (len(lst[ipdc])), '--', color = colorlist[ipdc])
            plt.legend(['PDC1', 'Mean1', 'PDC2', 'Mean2', 'PDC3', 'Mean3', 'PDC4', 'Mean4'])
            plt.title('TCR vs SPAD Index (Screamers Removed)')
            figure += 1
        plt.figure(figure)
        plt.plot(np.linspace(0,90,nSpads - nScreamers[ipdc]), sorted(lst[ipdc]), color = colorlist[ipdc])

    plt.legend(['PDC1', 'PDC2', 'PDC3', 'PDC4'])    
    plt.show()

#Uses percent method. Gives a scatter plot and percentage distribution for every small PDC
def noScreamersPerPlotBnA():
    df = pd.read_csv( "20250602_13h37m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
    nScreamers = []
    percent = 90
    lst=[[] for ipdc in range(nPdcs)]
    dfSort = [[] for ipdc in range(nPdcs)]
    for ipdc in range(nPdcs):
        for ispad in range(nSpads):
            dfSort[ipdc].append(df[f'SPAD_TCR{ipdc}'][ispad])
            lst[ipdc].append(df[f'SPAD_TCR{ipdc}'][ispad])
        dfSort[ipdc].sort()

    for ipdc in range(nPdcs):
        n = 0 #number of screamers in a PDC, resets at each new PDC.
        i = 0

        figure = 3
        while figure == 3: #Plot everything before taking out the screamers 
            plt.figure(figure)
            plt.scatter(range(0,nSpads), lst[ipdc], facecolors = 'none', edgecolors = colorlist[ipdc])
            plt.plot(range(0,nSpads), [np.mean(lst[ipdc])] * (len(lst[ipdc])), '--', color = colorlist[ipdc])
            plt.yscale('log')
            plt.legend(['PDC1', 'Mean1', 'PDC2', 'Mean2', 'PDC3', 'Mean3', 'PDC4', 'Mean4'])
            plt.title('TCR vs SPAD Index (Before)')
            figure += 1
        plt.figure(figure)
        plt.plot(np.linspace(0,90,nSpads), sorted(lst[ipdc]), color = colorlist[ipdc])
        plt.yscale('log')
        plt.title('Before')
        plt.legend(['PDC1', 'PDC2', 'PDC3', 'PDC4'])    

        while i < nSpads: #Take out the screamers
            if df[f'SPAD_percent{ipdc}'][i] >= percent:
                lst[ipdc].remove(dfSort[ipdc][i])
                n+=1
            i += 1
        nScreamers.append(n)

    for ipdc in range(nPdcs): #Plot without the screamers (After)
        figure = 1
        while figure == 1:
            plt.figure(figure)
            plt.scatter(range(0,nSpads - nScreamers[ipdc]), lst[ipdc], facecolors = 'none', edgecolors = colorlist[ipdc])
            plt.plot(range(0,nSpads - nScreamers[ipdc]), [np.mean(lst[ipdc])] * (len(lst[ipdc])), '--', color = colorlist[ipdc])
            plt.legend(['PDC1', 'Mean1', 'PDC2', 'Mean2', 'PDC3', 'Mean3', 'PDC4', 'Mean4'])
            plt.title('TCR vs SPAD Index (Screamers Removed)')
            figure += 1
        plt.figure(figure)
        plt.plot(np.linspace(0,90,nSpads - nScreamers[ipdc]), sorted(lst[ipdc]), color = colorlist[ipdc])
        plt.title('After')

    plt.legend(['PDC1', 'PDC2', 'PDC3', 'PDC4'])    
    plt.show()

#Iterative filter for screamers based on a multiple of the median. Gives a scatter plot and histogram for one small PDC for one run
def noScreamersMedianIterative():
    df = pd.read_csv( "20250602_13h37m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
    n = 0 #initialize screamer count
    i = 0 #gives the index of PDC to look at
    newdf = df
    cont = 'y'

    while cont == 'y':
        print(f'Median before iteration is {df[f'SPAD_TCR{i}'].median()}')

        for x in newdf.index:
            if newdf.loc[x, f'SPAD_TCR{i}'] > df[f'SPAD_TCR{i}'].median()*1.3:
                newdf.drop(x,inplace = True)
                n+=1

        print(f'Number of screamers after iteration is {n}')
        print(f'Median after iteration is {df[f'SPAD_TCR{i}'].median()}')

        plt.figure(1)
        plt.scatter(newdf.index, newdf[f'SPAD_TCR{i}'], facecolors = 'none', edgecolors = 'r')
        plt.plot([0,64], [newdf[f'SPAD_TCR{i}'].median()] * 2, '--', color = 'r')
        plt.figure(2)
        plt.hist(newdf[f'SPAD_TCR{i}'], bins = 'fd')
        plt.show()

        cont = input("continue? \n")

#Plots 1 run across all 4 small PDCs. Gives both a histogram and a gaussian curve fit
#first parameter is the height of the peak
#second parameter is the position of the center of the peak
#third parameter is the standard deviation
def thresholdGaussianFit():
    i=0 #indicate which file to use in dflist
    for ipdc in range(nPdcs):
        dflist = [
        pd.read_csv( "20250602_13h37m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250602_16h32m37_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250603_16h47m57_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
        ]
        for x in dflist[i].index:
            if dflist[i].loc[x, f'SPAD_TCR{ipdc}'] > 400:
                dflist[i].drop(x,inplace = True)

        plt.figure(ipdc)
        #plt.subplot(1,2,1)
        hist = plt.hist(dflist[i][f'SPAD_TCR{ipdc}'], bins = 'fd')
        plt.title(f'After, PDC{ipdc}')

        midpoint = []
        for k in range(len(hist[1])):
            if k != len(hist[1]) - 1:
                midpoint.append((hist[1][k] + hist[1][k+1]) / 2)
        #plt.subplot(1,2,2)
        #plt.scatter(midpoint, hist[0])

        def model_func(x,a,b,c):
            return a * np.exp(-(x - b)**2 / (2 * c**2))
        
        popt, pcov = curve_fit(model_func, midpoint, hist[0], p0 = [15,250,6])
        a_opt, b_opt, c_opt = popt
        plt.plot(midpoint, model_func(midpoint, *popt), 'r-')
        perr = np.sqrt(np.diag(pcov))
        #plt.title(f'{perr}')
        perr = np.sqrt(np.diag(pcov))
        print(f'PDC{ipdc}: {perr}')

    plt.show()

#Treats all PDCs from one file like a 64x4 array of spads
#Removes screamers by saying how many you want to take out
#Gives histogram with gauss fit and scatter plot 
def sumPdcsGaussFit():
    df = pd.read_csv("20250602_16h12m04_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
    lst = []
    for ipdc in range(4):
        for i in range(64):
            lst.append(df[f'SPAD_TCR{ipdc}'][i])

    removedlst=[]
    for i in range(50):
        removedlst.append(max(lst))
        lst.remove(max(lst))
    lowest = min(removedlst)
    print(lowest)

    plt.figure(1)
    plt.scatter(np.linspace(1,len(lst),len(lst)), lst, facecolors = 'none', edgecolors='r')
    plt.xlabel('SPAD Index')
    plt.ylabel('TCR')
    plt.yscale('log')
    plt.title('Sum of PDCs')

    plt.figure(2)
    hist = plt.hist(lst, bins = 'auto')
    plt.xlabel('TCR')
    plt.ylabel('Occurrences')
    plt.title('Sum of PDCs (PDC020)')

    midpoint = []
    for k in range(len(hist[1])):
        if k != len(hist[1]) - 1:
            midpoint.append((hist[1][k] + hist[1][k+1]) / 2)

    def model_func(x,a,b,c):
                return a * np.exp(-(x - b)**2 / (2 * c**2))
    
    popt, pcov = curve_fit(model_func, midpoint, hist[0], p0 = [30,250,10])
    a_opt, b_opt, c_opt = popt
    sigmaThresh = (lowest - b_opt) / c_opt
    print(f'n sigma = {sigmaThresh}')
    chiSquare = 0
    for i in range(len(midpoint)):
        chiSquare += ((model_func(midpoint, *popt)[i] - hist[0][i]) ** 2) / model_func(midpoint, *popt)[i]
    print(f'Chi-squared = {chiSquare}')
    print(f'sigma = {c_opt}')
    print(f'number of bins = {len(midpoint)}')

    plt.plot(midpoint, model_func(midpoint, *popt), 'r-')
    plt.show()

#Plots ccr of multiple runs
#Data must be in 'PDC0_CCR (%)' format
def ccrPlotter(dflist):
    fs = 16
    plt.rcParams.update({'font.size': 12})

    for i in range(len(dflist)):
        plt.scatter(dflist[i]['PDC0_SPAD_idx'], dflist[i]['PDC0_CCR (%)'], facecolors = 'none', edgecolors=colorlist[i], label = f'Run {i+1}')
        plt.xlabel('SPAD index', fontsize=fs)
        plt.ylabel('CCR [%]', fontsize = fs)
        plt.title('PDC0 CCR')
        plt.tick_params(length = 10, width = 2)
        plt.legend()
    plt.show()

#Gives the number of spads that get a negative ccr value based on the holdoff time
#have to order dflist by holdoff time to get correct table
#returns the dataframe, if you want to view you have to print it
#Data must be in CCR format
def getNegativeCcrValues(dflist):
    ipdc = 0 #which pdc
    data = {
        'holdoff-value': [14.4, 140, 1400, 14000, 18100],
        'number-of-negative-CCR-SPADs': []
    }

    for i in range(len(dflist)):
        n = 0
        for j in range(len(dflist[i][f'PDC{ipdc}_SPAD_idx'])):
            if dflist[i][f'PDC{ipdc}_CCR (%)'][j] < 0:
                n+=1
        data['number-of-negative-CCR-SPADs'].append(n)
        
    return(pd.DataFrame(data))

#Displays 64x64 SPAD array tcr heatmap
#Assumes that all 64x64 SPADs are turned on 
def fullArrayHeatMap():
    matrix = [[0] * 64 for i in range(64)]
    df = pd.read_csv("PDC11_25V_TCRmap.txt", sep = '\t')

    for i in range(64): #row index
        for j in range(64): #column index
            matrix[i][j] = df.loc[(df['Xcoord'] == j) & (df['Ycoord'] == i)]['SPAD_TCR0']

    plot = plt.imshow(matrix, cmap='viridis', norm = colors.LogNorm(), origin = 'lower')
    plt.ylabel('Y-coordinate')
    plt.xlabel('X-coordinate')
    cbar = plt.colorbar(plot)
    cbar.set_label('TCR (Log scale)')
    plt.show()

#plots ucr of multiple runs
def ucrPlotter():
    for i in range(len(dflist)):
        plt.scatter(dflist[i]['PDC0_SPAD_idx'], dflist[i]['PDC0_UCR (cps)'], facecolors = 'none', edgecolors=colorlist[i])
        plt.title('PDC0 UCR')
    plt.show()

def getNegativeUcrValues(dflist):
    ipdc = 0 #which pdc
    data = {
        'holdoff-value': [14.4, 140, 1400, 14000, 18100],
        'number-of-negative-UCR-SPADs': []
    }

    for i in range(len(dflist)):
        n = 0
        for j in range(len(dflist[i][f'PDC{ipdc}_SPAD_idx'])):
            if dflist[i][f'PDC{ipdc}_UCR (cps)'][j] <= 0:
                n+=1
        data['number-of-negative-UCR-SPADs'].append(n)

    return(pd.DataFrame(data))

def tcrUcrCcrPlotter(dflist):
    plt.figure(1)
    for i in range(len(dflist)):
        plt.scatter(dflist[i]['PDC0_SPAD_idx'], dflist[i]['PDC0_CCR (%)'], facecolors = 'none', edgecolors=colorlist[i])
        plt.title('PDC0 CCR')

    plt.figure(2)
    for i in range(len(dflist)):
        plt.scatter(dflist[i]['PDC0_SPAD_idx'], dflist[i]['PDC0_UCR (cps)'], facecolors = 'none', edgecolors=colorlist[i])
        plt.title('PDC0 UCR')
    
    plt.figure(3)
    for i in range(len(dflist)):
        plt.scatter(dflist[i]['PDC0_SPAD_idx'], dflist[i]['PDC0_TCR (cps)'], facecolors = 'none', edgecolors=colorlist[i])
        plt.title('PDC0 TCR (log scale)')
        plt.yscale('log')
    
    plt.show()

#Gives scatter plot and histogram for one PDC summed across multiple runs
def sumRunsGaussFit(dflist):
    nRemove = 14 #how many to remove. each increment of 1 removes 1 from each run, so it multiplies by the number of files.
    threshold = 380 #Any SPADs above this value will be removed
    ipdc = 3 #Which PDC to look at
    whichMethod = 2 # number method = 0; threshold method = 1; percent method = 2
    percent = 85
    fs=16 #Font size for the plots

    plt.rcParams.update({'font.size': 12})
    plt.tick_params(length = 10, width = 2)

    #Put all runs into a list to make a histogram from
    lst = []
    plt.figure(0)
    for i in range(len(dflist)): #scatter plot before screamers removed
        plt.subplot(1,2,1)
        plt.scatter(dflist[i][f'SPAD_idx{ipdc}'], dflist[i][f'SPAD_TCR{ipdc}'], facecolors = 'none', edgecolors = colorlist[i])
        plt.ylabel('TCR',fontsize=fs)
        plt.xlabel('SPAD index',fontsize=fs)
        plt.yscale('log')

    nScreamers = 0
    match whichMethod: #remove screamers from dataframe based on chosen method
        case 0: 
            print(f'For each run, {nRemove} of the highest TCR SPADs will be removed from PDC{ipdc}')
            for k in range(nRemove): #number of screamers method
                for i in range(len(dflist)):
                    dflist[i].drop(dflist[i][f'SPAD_TCR{ipdc}'].idxmax(), inplace = True)
            print(f'{nRemove * len(dflist)} removed in total')
        case 1:
            print(f'Any SPADs with a TCR above {threshold} will be removed from PDC{ipdc}')
            for i in range(len(dflist)): #threshold method
                dflist[i] = dflist[i][dflist[i][f'SPAD_TCR{ipdc}'] < threshold]
                nScreamers += 64 - len(dflist[i])
                print(f'{64 - len(dflist[i])} removed for run {i}')
            print(f'{nScreamers} removed in total')
        case 2:
            print(f'For each run, the SPADs with the highest {100-percent}% of TCR will be removed from PDC{ipdc}')
            for i in range(len(dflist)): #percent method
                newdf = dflist[i]
                newdf = newdf[newdf[f'SPAD_percent{ipdc}'] > percent]
                nScreamers += len(newdf)
                for j in newdf[f'SPAD_idx{ipdc}']:
                    dflist[i] = dflist[i][dflist[i][f'SPAD_TCR{ipdc}'] != newdf[f'SPAD_distribution{ipdc}'][j]] #keeps any spads <= the percentage
                
                #calculate lowest value that is thrown out
                if i == 0:
                    lowestThresh = newdf[f'SPAD_distribution{ipdc}'].min()
                    lowest = dflist[i][f'SPAD_TCR{ipdc}'].min()
                if newdf[f'SPAD_distribution{ipdc}'].min() < lowestThresh:
                    lowestThresh = newdf[f'SPAD_distribution{ipdc}'].min()
                if dflist[i][f'SPAD_TCR{ipdc}'].min() < lowest:
                    lowest = dflist[i][f'SPAD_TCR{ipdc}'].min()

            print(f'{nScreamers} removed in total')   

    for i in range(len(dflist)): #loop over every run
        for j in dflist[i][f'SPAD_idx{ipdc}']: #loop over every SPAD
            lst.append(dflist[i][f'SPAD_TCR{ipdc}'][j]) #give lst all the values from the modified dataframes

    for i in range(len(dflist)): #Make a scatter plot after screamers are removed
        plt.figure(0)                          
        plt.subplot(1,2,2)
        plt.scatter(dflist[i][f'SPAD_idx{ipdc}'], dflist[i][f'SPAD_TCR{ipdc}'], facecolors = 'none', edgecolors = colorlist[i])
        plt.ylabel('TCR',fontsize=fs)
        plt.xlabel('SPAD index',fontsize=fs)
        plt.title('TCR, Screamers Removed')
        plt.tick_params(length = 10, width = 2)

    #Create histogram and calculate midpoint of each bin    
    plt.figure(1)
    hist = plt.hist(lst, bins = 'auto', ec = 'b')
    plt.xlabel('')
    midpoint = []
    for k in range(len(hist[1])):
        if k != len(hist[1]) - 1:
            midpoint.append((hist[1][k] + hist[1][k+1]) / 2)

    def model_func(x,a,b,c): #Gaussian distribution
        return a * np.exp(-(x - b)**2 / (2 * c**2))
    
    popt, pcov = curve_fit(model_func, midpoint, hist[0], p0 = [30,250,10])
    a_opt, b_opt, c_opt = popt

    if whichMethod == 2:
        sigmaThresh = (lowestThresh - b_opt) / c_opt
        sigmaLowest = (b_opt - lowest) / c_opt
    else:
        sigmaThresh=0

    chiSquare = 0
    for i in range(len(midpoint)): #calculate chi squared
        chiSquare += ((model_func(midpoint, *popt)[i] - hist[0][i]) ** 2) / model_func(midpoint, *popt)[i]

    print('----------------------------')
    print(f'Chi-Squared = {chiSquare}')
    print(f'Reduced chi-squared = {chiSquare/(len(midpoint)-1)}')
    print(f'Sigma = {c_opt}')
    print(f'Number of bins = {len(midpoint)}')
    print(f'sigma threshold = {sigmaThresh}')
    print(f'sigma lowest = {sigmaLowest}')
    

    plt.plot(midpoint, model_func(midpoint, *popt), 'r-')
    plt.xlabel('TCR',fontsize=fs)
    plt.title(f'Histogram of PDC{ipdc} TCR')
    plt.show()

#Plots chi-squared and standard devation based on a gaussian model fit vs the number of screamers removed across all 4 small PDCs of one big PDC
#Only looks at one run at a time
#Need to clean up the code at some point. It works, just very messy
def plotChiSigmaSumPdcs():

    def model_func(x,a,b,c):
        return a * np.exp(-(x - b)**2 / (2 * c**2))
    nscreamersTot = 35 #Anything above this value it is not gauranteed that chi-squared and sigma can be calculated
    chiList = []
    reducedChiList = []
    #Try reduced chi-squared. chi-square/bins-1
    #Calculate (threshold-mean)/sigma. threshold is lowest value you threw out. mean is also from the fit
    sigmaList = []
    screamersDict = {
        'PDC0': [],
        'PDC1': [],
        'PDC2': [],
        'PDC3': []
    }
    for nscreamers in range(nscreamersTot): #Each instance of this for loop creates a value for chi-square and sigma to be plotted
        if nscreamers > 8:
            #initialize list with screamers
            df = pd.read_csv("20250602_16h12m04_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
            df2 = pd.read_csv("20250602_16h12m04_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
            lst = []
            lst2 = []
            for ipdc in range(4):
                for i in range(64):
                    lst.append(df[f'SPAD_TCR{ipdc}'][i])
            for ipdc in range(4):
                for i in range(64):
                    lst2.append(df2[f'SPAD_TCR{ipdc}'][i])
            lstdf = pd.DataFrame(lst)
            

            #remove screamers
            for i in range(nscreamers):
                lstdf.drop(lstdf.idxmax())
                #if nscreamers == nscreamersTot - 1:
                    #print(lst2.index(max(lst)))
                    #print(lst.index(max(lst)))
                lst.remove(max(lst))
        
            plt.figure(0)
            hist = plt.hist(lst, bins = 'auto')
            midpoint = []
            for k in range(len(hist[1])):
                if k != len(hist[1]) - 1:
                    midpoint.append((hist[1][k] + hist[1][k+1]) / 2)
        
            popt, pcov = curve_fit(model_func, midpoint, hist[0], p0 = [30,250,10])
            a_opt, b_opt, c_opt = popt

            chiSquare = 0
            reducedChi = 0
            for i in range(len(midpoint)):
                chiSquare += ((model_func(midpoint, *popt)[i] - hist[0][i]) ** 2) / model_func(midpoint, *popt)[i]
                reducedChi += ((model_func(midpoint, *popt)[i] - hist[0][i]) ** 2) / model_func(midpoint, *popt)[i] / (len(midpoint) - 1)
            chiList.append(chiSquare)
            reducedChiList.append(reducedChi)
            sigmaList.append(c_opt)

    #This block calculates which SPADs are removed for each PDC. The total number removed across PDCs will be = to nscreamerstot
    df = pd.read_csv("20250602_16h12m04_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
    lst = []
    for ipdc in range(4):
        for i in range(64):
            lst.append(df[f'SPAD_TCR{ipdc}'][i])
    lstdf = pd.DataFrame(lst)
    df = pd.DataFrame(lst)
    for i in range(nscreamersTot): 
        lstdf.drop(lstdf[0].idxmax(), inplace = True)
    for i in lstdf.index:
        df.drop(i, inplace = True)
    pdclst = [pd.DataFrame(columns = [f'PDC{i}']) for i in range(4)]
    for i in range(4):
        for j in range(64):
            if j + 64 * i in df.index:
                pdclst[i].loc[j] = df[0][j + 64 * i]
        print(pdclst[i])
        print(f'number of screamers from PDC{i} = {len(pdclst[i])}')
        print('--------------------------------')

    plt.close(0)
    plt.figure(1)
    plt.plot(np.linspace(9,nscreamersTot, nscreamersTot-9), chiList, color='b', label = 'chi-square')
    plt.plot(np.linspace(9,nscreamersTot, nscreamersTot-9), reducedChiList, color='r', label = 'reduced chi-square')
    plt.plot(np.linspace(9,nscreamersTot, nscreamersTot-9), sigmaList, color='orange', label = 'standard deviation')
    plt.grid(True)
    plt.ylim(0,100)
    plt.legend()
    plt.xlabel('Number of Screamers Removed')
    plt.title('Sum of PDC020')
    plt.show()

#Provides A single plot to show the evolution of the average tcr over multiple runs
#Each line is a different tile of the pdc
def avgPlotter(dflist, nSpads = 64, nPdcs = 4, colorlist = ['r', 'b', 'g', 'orange', 'y', 'k']):
        avgList = [[] for ipdc in range(nPdcs)]
        plt.rcParams.update({'font.size': 12})

        for ipdc in range(nPdcs): 
                
                for i in range(len(dflist)):
                        avgList[ipdc].append(sum(dflist[i][f'SPAD_TCR{ipdc}'])/len(dflist[i][f'SPAD_TCR{ipdc}']))

                yerr=[i**0.5 for i in avgList[ipdc]]
                plt.plot([i+1 for i in range(len(dflist))], avgList[ipdc], 'o', markersize=4, color = colorlist[ipdc], label = f'Tile{ipdc}')
                plt.plot([i+1 for i in range(len(dflist))], avgList[ipdc], '--', color = colorlist[ipdc])
                plt.errorbar([i+1 for i in range(len(dflist))], avgList[ipdc], yerr=yerr, fmt='none', ecolor=colorlist[ipdc])
        
        fs = 16
        plt.xticks([i+1 for i in range(len(dflist))])
        plt.ylabel('Mean TCR',fontsize=fs)
        plt.xlabel('Run #',fontsize=fs)
        plt.legend(fontsize=fs)
        plt.title('PDC020')
        plt.tick_params(length = 10, width = 2)
        plt.show()

#Provides a single plot of the evolution of the average tcr as one of the timing parameters changes
#Each line is a different tile of the pdc
#The order of the values for the timing must correspond to the order of the files in dflist
def avgTimingPlotter(dflist, nSpads = 64, nPdcs = 4, colorlist = ['r', 'b', 'g', 'orange', 'y', 'k']):
        #timing = [3.4,14.4, 25.4,36.4,47.4,58.6]
        timing = [14.4,140,1400,14000,18100]
        avgList = [[] for ipdc in range(nPdcs)]
        fs = 15.5 #Font size for plots
        plt.rcParams.update({'font.size': 12})

        for ipdc in range(nPdcs):
                
                for i in range(len(dflist)):
                        avgList[ipdc].append(sum(dflist[i][f'SPAD_TCR{ipdc}'])/len(dflist[i][f'SPAD_TCR{ipdc}']))
                plt.plot(timing, avgList[ipdc], 'o', color = colorlist[ipdc], label = f'Tile{ipdc}')
                plt.plot(timing, avgList[ipdc], '--', color = colorlist[ipdc])

        plt.legend()
        plt.xscale('log')
        plt.ylabel('Mean TCR', fontsize = fs)
        plt.xlabel('Hold-off Time [ns]', fontsize = fs)
        plt.title('PDC020')
        plt.tick_params(length = 10, width = 2)
        plt.tick_params(length = 10, width = 1, which = 'minor')
        plt.show()

def differencePlotter(dflist, colorlist = ['r', 'b', 'g', 'orange', 'y', 'k']):
    plt.rcParams.update({'font.size': 12})
    #ipdc = 2
    nSpads = 64
    averagelst = []
    differencelst = []

    for ipdc in range(4):
        averagelst = []
        differencelst = []
        for iSpad in range(nSpads):
            tcrlst = [dflist[run][f'SPAD_TCR{ipdc}'][iSpad] for run in range(len(dflist))]
            largest = max(tcrlst)
            smallest = min(tcrlst)
            difference = largest-smallest
            average = sum(tcrlst) / len(tcrlst)
            averagelst.append(average)
            differencelst.append(difference)
        plt.figure(0)
        plt.scatter(averagelst, differencelst, edgecolors = colorlist[ipdc], facecolors = 'none', label=f'PDC{ipdc}')
        #plt.figure(ipdc+1)
        #plt.scatter(averagelst, differencelst, edgecolors = colorlist[ipdc], facecolors = 'none', label=f'PDC{ipdc}')
        #plt.xscale('log')
        plt.xlabel('Average TCR')
        plt.ylabel('TCR Difference (min-max)')
        plt.legend()
        plt.tick_params(length = 8, width = 2)
        plt.tick_params(length = 6, width = 1, which = 'minor')
        plt.plot(np.linspace(200,30000), [8*i**0.5 for i in np.linspace(200,30000)])
    plt.title('PDC020 (y=8*x^0.5)')
    plt.show()


avgPlotter(dflist=dflist)