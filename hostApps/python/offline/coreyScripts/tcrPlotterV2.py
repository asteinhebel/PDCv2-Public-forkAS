#This Program is meant to be a compilation of all the working plots that I have to analyze the PDCs
#Corey Fox
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.optimize import curve_fit
from matplotlib import colors

nPdcs = 4
nSpads = 64
colorlist = ['r', 'b', 'g', 'orange', 'y', 'k']
#Use this dflist if there is not one specified in the function
dflist = [pd.read_csv("20250618_09h50m10_ZPP_PDC0_100ms.csv", sep = ';'),
          pd.read_csv("20250618_10h02m46_ZPP_PDC0_100ms.csv", sep = ';'),
          pd.read_csv("20250618_10h12m32_ZPP_PDC0_100ms.csv", sep = ';'),
          pd.read_csv("20250618_10h21m03_ZPP_PDC0_100ms.csv", sep = ';'),
          pd.read_csv("20250618_10h29m34_ZPP_PDC0_100ms.csv", sep = ';'),
          ]


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
def tcrPlotterMultipleRuns():
    dflist = [
        pd.read_csv( "20250604_14h43m19_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250604_14h52m39_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250604_15h03m36_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250604_15h07m18_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250604_15h09m21_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250604_15h25m41_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
    ]
    colorlist = ['r', 'b', 'g', 'orange', 'y', 'k']
    nmbRuns = len(dflist)
    nmbPdcs = 4
    legend = []
    colors = []

    for i in range(nmbRuns):
        legend.append(f'Run{i + 1}')
        legend.append(f'Mean{i + 1}')
        colors.append(colorlist[i % 6])

    for j in range(nmbPdcs): #j = pdc index
        for i in range(nmbRuns): #i = run index
            plt.figure(j)
            y_err = dflist[i][f'SPAD_TCR{j}'] ** 0.5
            plt.scatter(dflist[i][f'SPAD_idx{j}'], dflist[i][f'SPAD_TCR{j}'], facecolors = 'none', edgecolors=colors[i])
            #plt.errorbar(dflist[i][f'SPAD_idx{j}'], dflist[i][f'SPAD_TCR{j}'], yerr = y_err, fmt = 'none', ecolor = colors[i])
            plt.plot(dflist[i][f'SPAD_idx{j}'], [dflist[i][f'SPAD_TCR{j}'].mean()]*64, '--', color=colors[i])
            plt.title(f'PDC{j} Across {nmbRuns} Runs')
            plt.yscale('log')
        plt.legend(legend, bbox_to_anchor=(1.02, 1.15))
        plt.ylabel('TCR')
        plt.xlabel('SPAD Index')

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
def noScreamersPerPlotMultipleRuns(df1, index):
    nScreamers = []
    df = df1
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
    df = pd.read_csv('20250602_13h45m49_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv', sep=';')
    lst = []
    for ipdc in range(4):
        for i in range(64):
            lst.append(df[f'SPAD_TCR{ipdc}'][i])

    for i in range(53):
        lst.remove(max(lst))

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
def ccrPlotter():
    for i in range(len(dflist)):
        plt.scatter(dflist[i]['PDC0_SPAD_idx'], dflist[i]['PDC0_CCR (%)'], facecolors = 'none', edgecolors=colorlist[i])
        plt.title('PDC0 CCR')
    plt.show()

#Gives the number of spads that get a negative ccr value based on the holdoff time
#have to order dflist by holdoff time to get correct table
#returns the dataframe, if you want to view you have to print it
#Data must be in CCR format
def getNegativeCcrValues():
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

def getNegativeUcrValues():
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

def tcrUcrCcrPlotter():
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
def sumRunsGaussFit():
    nRemove = 10 #how many to remove. each increment of 1 removes 1 from each run, so it multiplies by the number of files.
    threshold = 380 #Any SPADs above this value will be removed
    ipdc = 0 #Which PDC to look at
    whichMethod = 2 # number method = 0; threshold method = 1; percent method = 2
    percent = 80
    dflist = [
        pd.read_csv( "20250604_14h43m19_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250604_14h52m39_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250604_15h03m36_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250604_15h07m18_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250604_15h09m21_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';'),
        pd.read_csv( "20250604_15h25m41_TCR_PDC0_PDC1_PDC2_PDC3_200ms.csv", sep=';')
    ]

    #Put all runs into a list to make a histogram from
    lst = []
    plt.figure(0)
    for i in range(len(dflist)): #scatter plot before screamers removed
        plt.subplot(1,2,1)
        plt.scatter(dflist[i][f'SPAD_idx{ipdc}'], dflist[i][f'SPAD_TCR{ipdc}'], facecolors = 'none', edgecolors = colorlist[i])
        plt.ylabel('TCR')
        plt.xlabel('SPAD index')
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
                    dflist[i] = dflist[i][dflist[i][f'SPAD_TCR{ipdc}'] != newdf[f'SPAD_distribution{ipdc}'][j]]
            print(f'{nScreamers} removed in total')   

    for i in range(len(dflist)): #loop over every run
        for j in dflist[i][f'SPAD_idx{ipdc}']: #loop over every SPAD
            lst.append(dflist[i][f'SPAD_TCR{ipdc}'][j]) #give lst all the values from the modified dataframes

    for i in range(len(dflist)): #Make a scatter plot after screamers are removed
        plt.figure(0)                          
        plt.subplot(1,2,2)
        plt.scatter(dflist[i][f'SPAD_idx{ipdc}'], dflist[i][f'SPAD_TCR{ipdc}'], facecolors = 'none', edgecolors = colorlist[i])
        plt.ylabel('TCR')
        plt.xlabel('SPAD index')
        plt.title('TCR, Screamers Removed')

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

    chiSquare = 0
    for i in range(len(midpoint)): #calculate chi squared
        chiSquare += ((model_func(midpoint, *popt)[i] - hist[0][i]) ** 2) / model_func(midpoint, *popt)[i]
    print('----------------------------')
    print(f'Chi-Squared = {chiSquare}')
    print(f'Sigma = {c_opt}')
    print(f'Number of bins = {len(midpoint)}')

    plt.plot(midpoint, model_func(midpoint, *popt), 'r-')
    plt.xlabel('TCR')
    plt.title(f'Histogram of PDC{ipdc} TCR')
    plt.show()

sumRunsGaussFit()