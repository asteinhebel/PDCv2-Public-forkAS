#----------------------------------------------------------------------------------
#-- Company: GRAMS
#-- Designer: Tommy Rossignol
#--
#-- Create Date: 2025-12-03
#-- Description:
#--     analysis functions used for energy spectrum analysis
#--
#-- Dependencies:
#-- Revision:
#-- Revision 1.0 - File Created
#-- Additional Comments:
#----------------------------------------------------------------------------------
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

from modules.pandasHelper import *

#----------------------------------------------------------------------------------
# analysis functions
#----------------------------------------------------------------------------------
def last_nonzero(arr, axis=1, invalid_val=-1):
    mask = arr!=0
    val = arr.shape[axis] - np.flip(mask, axis=axis).argmax(axis=axis) - 1
    return np.where(mask.any(axis=axis), val, invalid_val)

def first_nonzero(arr, axis=1, invalid_val=-1):
    mask = arr!=0
    return np.where(mask.any(axis=axis), mask.argmax(axis=axis), invalid_val)

def get_fwhm(X, Y):
    half_max=max(Y)/2.0
    d = np.sign(np.array(Y)-half_max)+1; # zero values are less than half_max
    left_idx = first_nonzero(d, axis=0)
    right_idx = last_nonzero(d, axis=0)
    return X[right_idx] - X[left_idx]


#----------------------------------------------------------------------------------
# Gaussian fit functions
#----------------------------------------------------------------------------------
def gaus(x, a, x0, sigma):
    """
    # single gaussian function
    # x: independant data for the fit
    # a: amplitude (gain) of the gaussian
    # x0: position of the gaussian peak
    # sigma: width of the gaussian
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

# NOTE: Since function is used in a fit to find optimal parameters,
#       fitY0 and nParamsPerGaus are global constants to this module.
# fitY0 = 0: no offset used in the fit,
# fitY0 = 1: offset is optimized during fit
fitY0 = 0
# nParamsPerGaus: number of parameters used per gaussian fit. Default is 3 (gain, center, sigma)
nParamsPerGaus = 3
def multiGaussian(x, *params):
    """
    # multiple gaussians function
    """
    # y0, [gain, center, sigma], [gain, center, sigma]
    if (len(params)-fitY0) % nParamsPerGaus != 0:
        if fitY0:
            raise TypeError("params must start with y0, then [gain, center, sigma] for each gaussian")
        else:
            raise TypeError(f"number of params must be divided by {nParamsPerGaus}. [gain, center, sigma] for each gaussian")
    nGauss = (len(params) - fitY0) // nParamsPerGaus
    if fitY0:
        y0 = params[0]
        y = y0 * np.ones_like(x)
    else:
        y = np.zeros_like(x)
    for i in range(nGauss):
        #gain = params[fitY0+3*i]
        #center = params[fitY0+3*i+1]
        #sigma = params[fitY0+3*i+2]
        #y += gaus(x, a=gain, x0=center, sigma=sigma)
        gaussParams = params[fitY0+nParamsPerGaus*i:fitY0+nParamsPerGaus*(i+1)]
        y += gaus(x, *gaussParams)
    return y


def gausParams2df(arr, columnNames):
    """
    # convert an array to a DataFrame with a column for each parameter.
    # columnNames: list of the name of the columns to generate
    """
    arr = np.array(arr)
    return pd.DataFrame(arr[fitY0:].reshape((len(arr)-fitY0)//nParamsPerGaus, nParamsPerGaus), columns=columnNames)


#----------------------------------------------------------------------------------
# Histogramming analysis functions
#----------------------------------------------------------------------------------

def rebinHisto(binCenters, binCounts, nBinPerBin):
    """
    # change the width of the bins of a histogram
    # binCenters: center of the bins to rebin
    # binCounts: count for each bin
    # nBinPerBin: number of original bins to regroup into the rebinned histogram
    """
    if binCounts.size % nBinPerBin != 0:
        # Truncate the arrays
        binCounts = np.array(binCounts[:-(binCounts.size % nBinPerBin)])
        binCenters = np.array(binCenters[:-(binCenters.size % nBinPerBin)])
    binCountsRebin = np.sum(binCounts.reshape(-1, nBinPerBin), axis=1)/nBinPerBin
    binCentersRebin = np.sum(binCenters.reshape(-1, nBinPerBin), axis=1)/nBinPerBin
    return binCentersRebin, binCountsRebin


#----------------------------------------------------------------------------------
# peak analysis functions
#----------------------------------------------------------------------------------
def cleanPeaks(dfPeaks, tolerance=50, nRecurs=5):
    """
    # Remove peaks from python "find_peaks" function
    # within tolerance by keeping the peak with the higher amplitude (gain).
    # dfPeaks must at least contain columns "gain", "center"
    """
    if not {"gain", "center"}.issubset(dfPeaks.columns):
        raise KeyError("dfPeak must contain columns 'gain' and 'center'")
    dfPeaks["tooClose"] = dfPeaks["center"].diff(periods=1).fillna(value=np.inf) <= tolerance
    
    if dfPeaks["tooClose"].any() and nRecurs>0:
        rowsToDrop = []
        lastIndex = -1
        for index, peakRow in dfPeaks.iterrows():
            if peakRow["tooClose"] and lastIndex != -1:
                if dfPeaks.loc[lastIndex, "gain"] > peakRow["gain"]:
                    # remove current index
                    rowsToDrop.append(index)
                else:
                    # remove previous index
                    rowsToDrop.append(lastIndex)
            lastIndex = index
        # remove duplicates
        rowsToDrop = list(set(rowsToDrop))
        #print(f"rowsToDrop:{rowsToDrop}")
        dfPeaks.drop(rowsToDrop, inplace=True, errors='ignore')
        # recursion
        dfPeaks = cleanPeaks(dfPeaks, tolerance=tolerance, nRecurs=nRecurs-1)
    dfPeaks.drop(columns="tooClose", axis=1, errors='ignore', inplace=True)
    
    return dfPeaks.reset_index(drop=True)














