#----------------------------------------------------------------------------------
#-- Company: GRAMS
#-- Designer: Tommy Rossignol
#--
#-- Create Date: 2026-03-04
#-- Description:
#--      Analyse and display csv results files from getDsumXt.py
#--
#-- Dependencies:
#-- Revision:
#-- Revision 1.0 - File Created
#--     Updated from Jupyter Notebook version created 2025-04-25
#-- Additional Comments:
#--     NOTE: recommended usage: Open an interactive python terminal
#--           and call the required functions
#----------------------------------------------------------------------------------
# python standard imports
import os, sys, re, ast
import inspect
import math
import random
import matplotlib.colors
import pandas as pd
import numpy as np
import time, datetime

import scipy.stats

try:
    from scipy.signal import savgol_filter, find_peaks, peak_prominences
    from scipy.optimize import curve_fit
    from scipy.signal import find_peaks
    from scipy import signal
except ImportError:
    print(f"ERROR: 'scipy' is not available. Please add it to your environment.")

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors

import numpy as np
import matplotlib.pyplot as plt
import scipy
import seaborn as sns

# custom modules
from modules.fgColors import fgColors
from modules.systemHelper import sectionPrint
import modules.pandasHelper as pdh
import modules.hexReadCsvParser as hrcp
import modules.energySpectrumAnalysisHelper as esah



# -----------------------------------------------
# --- list available functions
# -----------------------------------------------
def listAvailFunc(doPrint=False):
    """
    find all functions in this file
    """
    with open(__file__, 'r') as f:
        content = f.read()

    # ^def        : Matches 'def' at the absolute start of a line (column 0)
    # [^_]        : Ensures the function name doesn't start with an underscore
    # .*?         : Matches everything (including newlines) lazily...
    # (?=:\s*$)   : ...until it finds a colon at the end of a line
    pattern = r'^(def [^_].*?:(?=\s*$))'

    # re.MULTILINE: ^ matches start of lines
    # re.DOTALL: . matches newline characters
    matches = re.finditer(pattern, content, re.MULTILINE | re.DOTALL)
    funcList = [m.group(1).strip() for m in matches]

    if doPrint:
        for func in funcList:
            print(func)

    return funcList


# -----------------------------------------------
# --- Log execution time
# -----------------------------------------------
executionTime = {}
def logExecutionTime(label):
    """
    # keep the execution time of each step into a dictionary
    # label -> name of the key in the dictionary
    # E.g. "start", "import", "filter", ...
    """
    executionTime[label] = time.time()

def printExecutionTime():
    """
    # Take the dictionary with all time entries and
    # format it into a DataFrame with absolute time and each step time.
    """
    if "start" in executionTime:
        dfTime = pd.DataFrame.from_dict(executionTime,
                                        orient='index',
                                        columns=["time (sec)"])
        dfTime["time (sec)"] -= dfTime.loc["start"].values
        dfTime["dt (sec)"] = dfTime["time (sec)"].diff().fillna(value=0)
        first = dfTime.iloc[0, 0]
        last = dfTime.iloc[-1, 0]
        dfTime.loc["Total execution time"] = {"time (sec)": last, "dt (sec)": last-first}

        print()
        print('='*75)
        print(dfTime)
        print('='*75)


# -----------------------------------------------
# --- Import CSV data into a DataFrame
# -----------------------------------------------
def read_csv(datafile) -> pd.DataFrame:
    """
    # load CSV into a pandas dataframe
    """
    dfOut = hrcp.parseCsv(filename=datafile)

    # show the statistics of usage of the dataframe
    pdh.print_df_usage(dfOut)

    return dfOut


# -----------------------------------------------
# --- Extract Date and Time
# -----------------------------------------------
def getMeasureDateFromFileName(datafile) -> str:
    """
    # extract date from file path/name
    """
    date_pattern = r"\d{4}\d{2}\d{2}"
    dates = re.findall(date_pattern, datafile)
    if len(dates) == 0:
        raise ValueError
    else:
        measureDate = dates[0]
    return measureDate

def getMeasureTimeFromFileName(datafile) -> str:
    """
    # extract time from file path/name
    """
    time_pattern = r"\d{2}h\d{2}m\d{2}"
    times = re.findall(time_pattern, datafile)
    if len(times) == 0:
        raise ValueError
    else:
        measureTime = times[0]
    return measureTime


# -----------------------------------------------
# --- Extract number of events info
# -----------------------------------------------
def getNumberOfEvents(dfIn) -> int:
    """
    # get number of events
    """
    nEvents = dfIn['frameIdx'].nunique()
    print(f"nEvents = {nEvents}")
    return nEvents

def getEventList(dfIn) -> np.ndarray:
    """
    # get a list with the index of the events
    """
    eventList = dfIn['frameIdx'].unique()
    return eventList

def getFirstEventFromEventList(eventList) -> int:
    """
    # get the index of the first event in the list
    """
    return min(eventList)

def getLastEventFromEventList(eventList) -> int:
    """
    # get the index of the last event in the list
    """
    return max(eventList)


# -----------------------------------------------
# --- Extract PDC info
# -----------------------------------------------
def getPdcIdxInData(dfIn, doPrint=True) -> list:
    """
    # get which PDC indexes are found in DataFrame
    """
    pdcIdxFound = sorted(dfIn['pdcIdx'].unique())
    if doPrint:
        print(f"pdcIdxFound = {pdcIdxFound}")
    return pdcIdxFound


# -----------------------------------------------
# --- Extract the measurement duration
# --- and the event rate.
# -----------------------------------------------
def getMeasureDuration(dfIn, pdcIdxToUse:int=None, clkPrd=10e-9, doPlot=False):
    """
    # From timestamps of each digital sum bin, find the total acquisition time
    # pdcIdxToUse -> Using only a single PDC for the analysis.
    #                NOTE: The index selected must be present in the DataFrame.
    #                When set to None, using the index of the first PDC found.
    # clkPrd -> period of the system clock in seconds (use 10e-9 to specify ns)
    # doPlot -> When set to True, show the plot of the timestamps.
    """
    pdcList = getPdcIdxInData(dfIn, doPrint=False)
    if pdcIdxToUse is None:
        # Using the first PDC found in the DataFrame
        pdcIdxToUse = pdcList[0]
    elif not pdcIdxToUse in pdcList:
        # wrong parameter
        raise ValueError(f"PDC{pdcIdxToUse} is not in DataFrame ({pdcList})")

    # Extract the dataIdx of each digital sum bin
    # NOTE: Since the measurement has been saved with hexRead BIN_IDX_MODE="time",
    #       each "dataIdx" represents a clkPrd cycle.
    pdcTime = dfIn.loc[dfIn["pdcIdx"] == pdcIdxToUse, "dataIdx"].astype(np.int64)

    # Count the number of times the global counter overflowed (global counter is 32 bits)
    gblCntCycles = np.sum(np.diff(pdcTime) < -2**31)
    #print(f"gblCntCycles: {gblCntCycles}")

    # Start with the first timestamp, add the number of global counter cycles, and the number of cycles to the end.
    #             | from first timestamp | from gbl min to    | The number of global |
    #             | to gbl max value le  | the last timestamp | counter cycles       |
    acqClkCycles = 2**32-pdcTime.iloc[0] +  pdcTime.iloc[-1]  + 2**32*(gblCntCycles-1)

    # convert the number of clock cycles to a time
    acqTime = clkPrd*acqClkCycles
    #print(f"gblCntCycles={gblCntCycles}, acqClkCycles={acqClkCycles}, acqTime={acqTime}, acqTime={acqTime/60.0}")
    if doPlot:
        plt.figure()
        #plt.plot(gblCntCycles)
        plt.plot(pdcTime)
    print(f"Test duration:\n  {acqTime:.3f} seconds \n  {acqTime/60:.3f} min \n  {acqTime/3600:.3f} hours")
    return acqTime

def getMeasureEventRate(nEvents, acqTime):
    """
    # Calculate the rate from the number of events and the duration of the acquisition
    # nEvents -> number of events in the DataFrame. Calculated using getNumberOfEvents().
    # acqTime -> duration of the measure in the DataFrame. Calculated using getMeasureDuration(). 
    """
    nEventsPerSec = nEvents / acqTime
    print(f"{nEventsPerSec:0.01f} events per sec")
    return nEventsPerSec


# -----------------------------------------------
# --- Filter data from DataFrame
# -----------------------------------------------
def groupDfByFrameAndPdc(dfIn) -> pd.core.groupby.DataFrameGroupBy:
    """
    # group DataFrame by frameIdx and pdcIdx
    # The resulting group helps for further operations
    # on data from a same event frame from a PDC.
    # NOTE: Use this function to generate "dfGroup"
    #       required for the following functions.
    """
    return dfIn.groupby(["frameIdx", "pdcIdx"], as_index=False)

def filterFrameFromMaxBin(dfGroup, binMaxThGt) -> tuple:
    """
    # keep only frames in which a bin is higher than threshold.
    # binMaxThGt -> greater-than threshold value of the digital sum.
    # Within an event frame, at least one bin must exceed this threshold
    # to keep the frame for further analysis.
    # E.g. binMaxThGt=6: Frames with no bin higher than 6 are removed from analysis.
    """
    keepFramesOverNoise = dfGroup['dsum'].transform("max").gt(binMaxThGt)
    nKeep, nTotal, pctKeep, pctDrop = pdh.filterStats(keepFramesOverNoise,
                                                      description=inspect.currentframe().f_code.co_name)
    return keepFramesOverNoise, (nKeep, nTotal, pctKeep, pctDrop)

def filterFrameFromNumBin(dfGroup, numBinThGt) -> tuple:
    """
    # keep only frames on which more than N bins contain data
    # numBinThGt -> greater-than threshold of the number digital sum bin in a frame.
    # Within an event frame, the number of bins must exceed the threshold
    # to keep the frame for further analysis.
    # E.g. numBinThGt=2: Frames with 2 bins of less are removed from analysis.
    """
    keepFramesNumBins = dfGroup['dsum'].transform("size").gt(numBinThGt)
    nKeep, nTotal, pctKeep, pctDrop = pdh.filterStats(keepFramesNumBins,
                                                      description=inspect.currentframe().f_code.co_name)
    return keepFramesNumBins, (nKeep, nTotal, pctKeep, pctDrop)

def filterFrameFromTotalPhoton(dfGroup, numTotalPhThMin=0, numTotalPhThMax=4096) -> tuple:
    """
    # keep only frames based on the total number of detected photons in the frame.
    #   It applies a range to a preliminary estimation of the energy
    # numTotalPhThMin -> minimum value (inclusive) of the total number of detected photon a frame.
    # numTotalPhThMax -> maximum value (inclusize) of the total number of detected photon a frame.
    # Within an event frame, the sum of all digital sum bins must be in range
    # to keep the frame for further analysis.
    # E.g. numTotalPhThMin=100, numTotalPhThMax=2500: Frames with less than 100 or more than 2500
    # are removed from analysis.
    """
    prelimEnergy = dfGroup['dsum'].transform("sum")
    keepEnergyRngH = prelimEnergy.ge(numTotalPhThMin)
    keepEnergyRngL = prelimEnergy.le(numTotalPhThMax)
    keepEnergyRng = np.logical_and(keepEnergyRngH, keepEnergyRngL)
    nKeep, nTotal, pctKeep, pctDrop = pdh.filterStats(keepEnergyRng,
                                                      description=inspect.currentframe().f_code.co_name)
    return keepEnergyRng, (nKeep, nTotal, pctKeep, pctDrop)

def filterFramesByIndex(dfGroup, begin:int=None, end:int=None, eventList:list=None) -> tuple:
    """
    # keep only a portion of all the dataset
    # begin -> first index of the event to keep (None = begin with first)
    # end -> last index of the event to keep (None = end with last)
    # eventList -> specify the list of events (None = recalculate the list)
    # E.g. begin=None, end=None -> all data
    # E.g. begin=None, end=1000 -> begin of the data up to 1000th event
    # E.g. begin=-1000, end=None -> from 1000th event to the end
    """
    if eventList is None:
        # no list specified, recalculate
        eventList = getEventList(dfGroup.obj)

    if begin is None:
        # start from the beginning of the data
        idxBegin = getFirstEventFromEventList(eventList)
    else:
        # specified an integer index
        try:
            idxBegin = eventList[begin]
        except IndexError:
            printf(f"ERROR: begin ({begin}) is out of eventList range (0:{len(eventList)-1}), using all events.")
            idxBegin = getFirstEventFromEventList(eventList)
    #print(f"idxBegin: {idxBegin}")

    if end is None:
        # start from the beginning of the data
        idxEnd = getLastEventFromEventList(eventList)
    else:
        # specified an integer index
        try:
            # index from the end
            idxEnd = eventList[end]
        except IndexError:
            printf(f"ERROR: end ({end}) is out of eventList range (0:{len(eventList)-1}), using all events.")
            idxEnd = getLastEventFromEventList(eventList)
    #print(f"idxEnd: {idxEnd}")
    #print(f"rng: {idxEnd-idxBegin}")

    keepFrameIdxRngH = dfGroup['frameIdx'].transform("first").gt(idxBegin)
    keepFrameIdxRngL = dfGroup['frameIdx'].transform("last").lt(idxEnd)
    keepFrameIdxRng = np.logical_and(keepFrameIdxRngH, keepFrameIdxRngL)
    nKeep, nTotal, pctKeep, pctDrop = pdh.filterStats(keepFrameIdxRng,
                                                      description=inspect.currentframe().f_code.co_name)
    return keepFrameIdxRng, (nKeep, nTotal, pctKeep, pctDrop)

"""
# Examples of function call
filterFramesByIndex(dfPdcGroup, begin=None, end=None)
filterFramesByIndex(dfPdcGroup, begin=None, end=1000)
filterFramesByIndex(dfPdcGroup, begin=-1000, end=None)
"""

def filterFramesByPdcIndex(dfIn, pdcToKeep: list) -> tuple:
    """
    # keep frames based on PDC index
    # pdcToKeep -> list of PDCs in frame to keep.
    # E.g. pdcToKeep=[0,1,2,3] -> all PDCs from first head
    # E.g. pdcToKeep=[4,5,6,7] -> all PDCs from second head
    # E.g. pdcToKeep=[0] -> only PDC 0
    """
    keepPdcs = dfIn["pdcIdx"].isin(pdcToKeep)
    nKeep, nTotal, pctKeep, pctDrop = pdh.filterStats(keepPdcs,
                                                      description=inspect.currentframe().f_code.co_name)
    return keepPdcs, (nKeep, nTotal, pctKeep, pctDrop)

"""
# Examples of function call
filterFramesByPdcIndex(dfIn=df, pdcToKeep=[0,1,2,3])
filterFramesByPdcIndex(dfIn=df, pdcToKeep=[0,1])
"""

def combineFilters(filterList: list, operation="and"):
    """
    # With all the filters generated by the previous functions,
    # combine them into a single filter.
    # filterList -> List of filter masks generated with previous funtions.
    # operation -> Logical operation to apply. Supported values are "and", "or"
    # NOTE: The function works if the list contain a single item,
    #       but is useless, because the resulting filter is the same as the input.
    # E.g. filterList=[keepFrameIdxRng, keepEnergyRng, keepFramesNumBins]
    """
    if not isinstance(filterList, list):
        raise TypeError("'filterList' is expected to be a list of filters.\n"\
                        f"{' '*12}To specify a single item, use list operator '[]'.")
    supportedOperations = ["and", "or"]
    if not operation in supportedOperations:
        raise TypeError(f"Supported 'operation' options are {supportedOperations}")

    if operation == "and":
        logicalFunction = np.logical_and
    elif operation == "or":
        logicalFunction = np.logical_or

    combined = []
    for i, filt in enumerate(filterList):
        if i == 0:
            # use first item to initialize the filter
            combined = filt
        else:
            # apply 'logical and' on remaining items
            combined = logicalFunction(combined, filt)
    # return combined filters
    return combined

"""
# Examples of function call
combinedFilter = combineFilters([keepFrameIdxRng, keepFramesNumBins])
"""

def applyFilters(dfIn, dfFilter) -> pd.DataFrame:
    """
    # apply filters on a DataFrame to keep only items for which the filter is set to True
    # dfFilter -> Either one filter from the previous functions
    #             of a combined filter from "combineFilters".
    """
    return dfIn[dfFilter].copy()


# -----------------------------------------------
# --- Identify frames with more than one PDC
# --- and analyse for crosstalk
# -----------------------------------------------
def analyseNumPdcPerFrame(dfIn):
    """
    # For each frame, count the number of PDCs in this frame.
    # dfIn -> DataFrame on which to analyze the number of PDC per frame
    # NOTE: changes are directly applied to the DataFrame argument "dfIn"

    # Modify the DataFrame to add two new columns: "nPdcsInFrame" and "pdcsInFrame"
    # nPdcsInFrame: number of PDCs in the frame.
    # pdcsInFrame: bitwise value where each bit represents a PDC (PDC0 = LSB).
    # NOTE: This function must be used after properly filtering the input DataFrame.
    #       Some noise event frames (too few bins per frame, no bins larger than a given threshold)
    #       contain more than 1 PDC. They must be filtered first.

    """
    # for each frame, extract column "pdcIdx" and count the number of unique PDC index values
    # Each frame mostly 1 PDC per frame when using one scintillator per PDC
    dfIn.loc[:, "nPdcsInFrame"] = dfIn.groupby("frameIdx")["pdcIdx"].transform("nunique")

    # for most of the frames, since a single PDC is present, a vectorized operation works
    # transform pdcIdx to bitwise (PDC0=1, PDC1=2, PDC2=4, PDC3=8)
    # NOTE: casted as 8 bits unsigned integer to reduce memory usage.
    # For more a system with more than 8 PDCs, change the cast type.
    dtype = np.uint8
    dfIn.loc[:, "pdcsInFrame"] = 2**dfIn["pdcIdx"].astype(dtype)

    # find frames with more than 1 PDC
    nPdcMask = dfIn["nPdcsInFrame"] > 1

    # print statistics
    pdh.filterStats(nPdcMask, description=inspect.currentframe().f_code.co_name)

    # apply the proper function to combine the PDCs into a mask on 8 bits
    # E.g. if PDC0 and PDC1 are present: 0b00000001 + 0b00000010 = 0b00000011
    def maskPdc(x):
        return np.sum(x.unique(), dtype=dtype)
    dfIn.loc[nPdcMask, "pdcsInFrame"] = dfIn.loc[nPdcMask, :].groupby("frameIdx")["pdcsInFrame"].transform(maskPdc)

def analyseCrosstalkBetweenScintillators(dfIn, doPlot=True):
    """
    # Using which PDCs have data for each bin of data to evaluate crosstalk (light sharing)
    # between the scintillators (one scintillator per PDC, 4 scintillators)
    # doPlot -> True: show a plot, False: only print results in terminal
    NOTE: this function requires "analyseNumPdcPerFrame" to be run on the DataFrame
    """
    # verify "analyseNumPdcPerFrame" has been run on the DataFrame
    if not {"nPdcsInFrame", "pdcsInFrame"}.issubset(dfIn.columns):
        raise KeyError("First call 'analyseNumPdcPerFrame' on 'dfIn'"\
                       "to generate columns 'nPdcsInFrame' and 'pdcsInFrame'")

    from collections import Counter
    # variable incremented in sub function histogramXt,
    # to count the total number of bins
    global cumul
    cumul = 0
    if doPlot:
        fig, axes = plt.subplots(2, 2, figsize=(12,6))
        axes = axes.flatten()
    else:
        axes = [None]*4

    # function local to this function since they are not use anywhere else
    def bitwiseLabelToStr(label):
        labelStr = ""
        for iPdc in range(0,8):
            if label & (0x1<<iPdc):
                labelStr+=str(iPdc)
        return labelStr

    # function local to this function since they are not use anywhere else
    def allBtwiseLabelsToStr(labels):
        labelsOut = []
        for label in labels:
            labelsOut.append(bitwiseLabelToStr(label))
        return labelsOut

    def histogramXt(dfIn_, nPdc, ax):
        plural = "s" if nPdc > 1 else ""
        print(f"Event bins with {nPdc} PDC{plural}:")

        nTotal = dfIn_["frameIdx"].size
        # Keep only bins with the number of PDC specified (nPdc) in the frame
        xt = dfIn_.loc[dfIn_["nPdcsInFrame"] == nPdc, "pdcsInFrame"]
        if xt.size == 0:
            print(f"  No events detected.")
            return

        # Use Counter class from collection module to extract stats
        counts = Counter(xt)
        freqs  = 100.0*np.array(list(counts.values()))/nTotal
        xtTotal = np.sum(freqs)
        global cumul
        cumul += xtTotal

        # convert bisewise values of the keys into labels with index of PDCs
        # E.g. 7 -> 012 (PDC0, PDC1, PDC2)
        labels = allBtwiseLabelsToStr(counts.keys())

        # Print global statistics
        print(f"  Bins: {xt.size}/{nTotal} ({xtTotal:.6f} %), cumulative: {cumul:.6f} %")

        # sort by labels
        labels, freqs, values = zip(*sorted(zip(labels, freqs, counts.values())))
        for label, freq, val in zip(labels, freqs, values):
            print(f"    PDC{label}: {freq:.6f} % ({val})")

        if ax is not None:
            ax.bar(labels, freqs)
            ax.set_title(f"Event bins with {nPdc} PDC{plural}")
            ax.set_ylabel("% of total event bins")
            ax.set_xlabel("PDC groups")

    # print statistics for each number of PDC and plot a histogram if enabled by doPlot
    # NOTE: PDCs that are side to side have more coupling than the ones in diagonal.
    histogramXt(dfIn, nPdc=1, ax=axes[0])
    histogramXt(dfIn, nPdc=2, ax=axes[1])
    histogramXt(dfIn, nPdc=3, ax=axes[2])
    histogramXt(dfIn, nPdc=4, ax=axes[3])
    if doPlot:
        plt.tight_layout()


# -----------------------------------------------
# --- Preprocessing of digital sum data
# -----------------------------------------------
def getDsumDt(dfGroup):
    """
    # Calcul the time difference between each bin using dataIdx, since the zeros have been
    # suppressed in the original CSV file loaded into the DataFrame.
    # NOTE: changes are directly applied to the DataFrame refered by group "dfGroup"
    """
    dfGroup.obj.loc[:, "dt"] = dfGroup["dataIdx"].diff().fillna(value=1).astype(np.int16)

def getDsumBinIdx(dfGroup):
    """
    # For each frame, extract the index of the bin, based on dataIdx.
    # NOTE: Considering here the first data is in bin 0, which might not be true.
    # NOTE: changes are directly applied to the DataFrame refered by group "dfGroup"
    # Modify the DataFrame to add a new column: "binIdx"
    # binIdx: index of each bin in a frame (considering empty bins have been suppressed)
    """
    # _df is a reference to base Dataframe
    _df = dfGroup.obj

    # use "dataIdx" column to find the index of a bin in the frame
    _df.loc[:, "binIdx"] = (_df["dataIdx"].astype(np.int64) - dfGroup["dataIdx"].transform("min")).astype(np.int16)

def getDsumFramePeakInfo(dfGroup):
    """
    # For each frame, extract the position of the peak and its value.
    # NOTE: changes are directly applied to the DataFrame refered by group "dfGroup"
    # Modify the DataFrame to add new columns: "peakValue", "dataIdxPeak" and "dtPeak"
    # peakValue: value of the actual peak
    # dataIdxPeak: index of the peak
    # dtPeak: for each bin, the number of bins separating from the peak
    """
    # _df is a reference to base Dataframe
    _df = dfGroup.obj

    # Adding a column "peakValue" with the maximum value found in each PDC frame
    _df.loc[:, "peakValue"] = dfGroup["dsum"].transform("max").astype(np.uint16)

    # Adding a column "dataIdxPeak" with the index of the value at which the maximum occurs in the PDC frame
    _df.loc[:, "dataIdxPeak"] = _df.loc[dfGroup["dsum"].transform("idxmax").values, "dataIdx"].values

    # Adding a column "dtPeak" with the bin delta between the peak and the current bin
    _df.loc[:, "dtPeak"] = (dfGroup.obj["dataIdx"].astype(np.int64) - dfGroup.obj["dataIdxPeak"]).astype(np.int16)

def getDsumPrelimEnergy(dfGroup):
    """
    # For each frame, extract a preliminary energy value
    # NOTE: changes are directly applied to the DataFrame refered by group "dfGroup"
    # Modify the DataFrame to add a new column: "energy"
    # energy: for each frame, apply a summation of each digital sum bin
    """
    dfGroup.obj["energy"] = dfGroup["dsum"].transform("sum")


# -----------------------------------------------
# --- Sampling DataFrame
# -----------------------------------------------
def sampleDataFrame(dfIn, dfMaxFrames=200) -> pd.DataFrame:
    """
    From a DataFrame, return a sub DataFrame with a maximum of dfMaxFrames frames.
    dfMaxFrames -> maximum number of frames to keep
    """
    # if DataFrame is larger than the number of frames, only use a subset
    if dfIn.size > dfMaxFrames:
        print(f"Using a subset of the input DataFrame ({dfMaxFrames})")
        eventList = getEventList(dfIn)
        framesToKeep = eventList[np.linspace(0, len(eventList)-1, num=dfMaxFrames, dtype=int)]
        dfSampled = dfIn.loc[dfIn["frameIdx"].isin(framesToKeep), :]
    else:
        dfSampled = dfIn # view reference, not a copy
    return dfSampled


# -----------------------------------------------
# --- Using digital sum data to find hold-off
# -----------------------------------------------
def sortedDictFromPeaks(data, idx=None):
    """
    # From a list or array of data, find peaks, sort them (descending)
    # And return a dictionary
    # data -> data on which to find peaks
    # idx -> idx associated with the data.
    #        If idx==None, index are calculated as a continuous range
    """
    if idx is None:
        # set idx as a continuous range
        idx = range(len(data))
    # find peaks
    peaks, props = find_peaks(data)
    # building a dictionary
    peakValues = {}
    for peak in peaks:
        peakValues[idx[peak]] = data[peak]
    # sorting dictionary
    peakValues = dict(sorted(peakValues.items(), key=lambda x: x[1], reverse=True))
    return peakValues

def getHoldOffFromDsum(dfIn, dfMaxFrameToUse=20000,
                       method="all", tauClkCycle=4.0, doPlot=True) -> int:
    """
    # Estimate the hold-off time from the digital sum frames
    # dfIn -> DataFrame with digital sum to use to analyse the hold-off
    # dfMaxFrameToUse -> To speed up execution time, analysis is done on a subset of data.
    #                    Select here the maximum number of data to use for analysis.
    #                    The kept frames are equally spaced in the population.
    # method -> The method used to find the hold-off.
    #           "peak" : using the cumulated frame peaks
    #           "corr" : using cross-correlation of the cumulated frame to find hold-off
    #           "fit"  : using an exponential decay fit and looking for the maximum error value
    #           "all"  : using all methods
    # tauClkCycle -> when using method fit (fit or all), scintillation Tau used.
    #                It is represented in number of samples.
    #                It must be divided by the dsum clock period.
    # doPlot -> when set to True, will plot the selected method results
    # Returns the estimated hold-off value (in number of clock cycles).
    # If "all" method is used, the best value of the three methods is used.
    """
    supportedMethods = ["peak", "corr", "fit", "all"]
    if method not in supportedMethods:
        raise ValueError(f"Supported methods are: {supportedMethods}")

    # verify "getDsumFramePeakInfo" has been run on the DataFrame
    if not {"dtPeak"}.issubset(dfIn.columns):
        raise KeyError("First call 'getDsumFramePeakInfo' on 'dfIn'"\
                       "to generate column 'dtPeak'")

    # keep execution time
    start_time = time.time()

    ## if DataFrame is too large, only use a subset to reduce processing time
    dfSampled = sampleDataFrame(dfIn, dfMaxFrames=dfMaxFrameToUse)

    # group data
    dfSampledGroup = groupDfByFrameAndPdc(dfSampled)
    nGroups = dfSampledGroup.ngroups

    if doPlot:
        # get list of PDCs in DataFrame
        pdcIdxList = getPdcIdxInData(dfSampled, doPrint=False)
        if len(pdcIdxList) == 1:
            # single PDC in the DataFrame
            pdcLabel = f"PDC{pdcIdxList[0]}"
        else:
            pdcLabel = ""

    # the maximum number of bins in a PDC frame is 128
    maxDsumFrameSize = 128
    # since the peaks need to be aligned, use twice this size
    frameCumul = np.zeros((2*maxDsumFrameSize))
    idxCumul = np.arange(-maxDsumFrameSize, maxDsumFrameSize)
    for i, (groups, dsumFrame) in enumerate(dfSampledGroup):
        # printing progress
        if (i%100) == 0:
            print(f"\r{i}/{nGroups} ({100.0*i/nGroups:0.1f} %)", end="")
        #if len(dsumFrame) > 20:
        #if np.max(dsumFrame["dtPeak"]) > 20:
        frameCumul[dsumFrame["dtPeak"]+maxDsumFrameSize] += dsumFrame["dsum"]
    print(f"\r{i}/{nGroups} ({100.0*i/nGroups:0.1f} %)")

    # free memory, no longer required variables
    del dfSampled
    del dsumFrame

    # look for peaks only after the principal peak (decay)
    peak0 = np.argmax(frameCumul)
    xDecay = idxCumul[peak0+1:]
    yDecay = frameCumul[peak0+1:]

    # prepare figure and axes if required
    axPeak = None
    axCorr = None
    axFit = None
    if doPlot == True:
        fig = plt.figure(num=pdcLabel, figsize=(12,6))
        if method == "all":
            axPeak, axCorr, axFit = fig.subplots(1,3)
        elif method == "peak":
            axPeak = fig.subplots(1,1)
        elif method == "corr":
            axCorr = fig.subplots(1,1)
        elif method == "fit":
            axFit = fig.subplots(1,1)

    # analysis of the hold-off time using peaks in cumulative frame
    if method in ["peak", "all"]:
        print("Using 'peak' method:")
        # get peaks and sort them by value
        peakValues = sortedDictFromPeaks(data=yDecay, idx=xDecay)
        # get first item
        holdPeak = next(iter(peakValues.items()))
        HOLD_TIME_CLK_CYCLES = holdPeak[0]
        print(f"  HOLD_TIME_CLK_CYCLES: {holdPeak[0]}")

        if axPeak is not None:
            # update plot if required
            axPeak.plot(idxCumul, frameCumul, '.-', label="frameCumul")
            #axPeak.plot(idxCumul[peaks], frameCumul[peaks], 'x', color='red', label="peaks")
            axPeak.plot(peakValues.keys(), peakValues.values(), 'x', color='red', label="peaks")
            axPeak.plot([holdPeak[0], holdPeak[0]], [np.min(frameCumul), np.max(frameCumul)], '--', color="gray")
            axPeak.plot(holdPeak[0], holdPeak[1], 's', markerfacecolor='none', markeredgecolor="red", markersize=10, label=f"hold ({holdPeak[0]} clk cycles)")
            axPeak.set_yscale("log")
            axPeak.legend()
            axPeak.set_title("find peaks")


    # analysis of the hold-off time using cross-correlation of cumulative frame
    if method in ["corr", "all"]:
        print("Using 'corr' method:")
        # using autocorrelation to find repetitive patterns
        corr = np.correlate(yDecay, yDecay, mode="full")[len(yDecay)-1:]

        # get peaks and sort them by value
        peakValuesCorr = sortedDictFromPeaks(data=corr, idx=xDecay)
        holdCorrPeak = next(iter(peakValuesCorr.items()))
        HOLD_TIME_CLK_CYCLES = holdCorrPeak[0]
        print(f"  HOLD_TIME_CLK_CYCLES: {holdCorrPeak[0]}")

        if axCorr is not None:
            # update plot if required
            axCorr.plot(corr, '.-', label="autocorrelation")
            axCorr.plot(peakValuesCorr.keys(), peakValuesCorr.values(), 'x', color="red", label="peaks")
            axCorr.plot([holdCorrPeak[0], holdCorrPeak[0]], [np.min(corr), np.max(corr)], '--', color="gray")
            axCorr.plot(holdCorrPeak[0], holdCorrPeak[1], 's', markerfacecolor='none', markeredgecolor="red", markersize=10, label=f"hold ({holdCorrPeak[0]} clk cycles)")
            axCorr.set_yscale("log")
            axCorr.legend()
            axCorr.set_title("autocorrelation")

    if method in ["fit", "all"]:
        print("Using 'fit' method:")
        # Shift x to 0 for numerical stability during fitting
        xDecay0 = xDecay - xDecay[0]

        # Smart initial guesses: [amplitude, decay_rate]
        amplitude =  yDecay.max() - yDecay.min()
        p0 = [amplitude, tauClkCycle]

        # Find fit parameters
        try:
            popt, pcov = curve_fit(expDecay, xDecay0[:15], yDecay[:15], p0=p0)
            print(f"Tau (fit) = {popt[1]:.3f} samples")
        except RuntimeError as ex:
            print(f"    RuntimeError during curve_fit:\n{ex}")
            popt = p0

        # Apply fit and clip values
        clipMin = yDecay[yDecay > 0][-1]
        yDecayFit = np.clip(expDecay(xDecay0, *popt), a_min=clipMin, a_max=None)

        # find ratio between fit and data
        # using epsilon to prevent zero divisions and log(0)
        eps = 1e-12
        ratio = yDecay/(yDecayFit + eps)

        # fit flex point where data diverges from the fit
        gate = ratio > 1.5
        startPoint = xDecay[gate][0]
        xH = xDecay[startPoint-2:]
        yH = yDecay[startPoint-2:]

        peakValuesFit = sortedDictFromPeaks(data=yH, idx=xH)
        holdFitPeak = next(iter(peakValuesFit.items()))
        HOLD_TIME_CLK_CYCLES = holdFitPeak[0]
        print(f"  HOLD_TIME_CLK_CYCLES: {holdFitPeak[0]}")

        pltColors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        if axFit is not None:
            # update plot if required
            axFit.plot(idxCumul, frameCumul, '.-', color=pltColors[0], label="frameCumul")
            axFit.plot(xDecay, yDecayFit, '.-', color=pltColors[2], label="exponentialDecayFit")
            axFit.plot(xDecay, ratio, '.-', color=pltColors[1], label="ratio")
            axFit.plot([holdFitPeak[0], holdFitPeak[0]], [np.min(frameCumul), np.max(frameCumul)], '--', color="gray")
            axFit.plot(holdFitPeak[0], holdFitPeak[1], 'x', markeredgecolor="red", markersize=10)
            axFit.plot(holdFitPeak[0], holdFitPeak[1], 's', markerfacecolor='none', markeredgecolor="red", markersize=10, label=f"hold ({holdFitPeak[0]} clk cycles)")
            axFit.set_yscale("log")
            axFit.legend()
            axFit.set_title("exponential decay fit")

    if method == "all":
        print("With 'all' method: select the proper value")
        holdDict = {
            "peak": holdPeak[0],
            "corr": holdCorrPeak[0],
            "fit": holdFitPeak[0]
        }
        print(f"average : {np.round(np.average(list(holdDict.values()))):.0f}")
        # use the median to remove the outlier (if there is any)
        holdDict = sorted(holdDict.items(), key=lambda item: item[1])
        print(f"median  : {holdDict[1][1]} -> {holdDict[1][0]}")
        HOLD_TIME_CLK_CYCLES = holdDict[1][1]
        print(f"  HOLD_TIME_CLK_CYCLES: {HOLD_TIME_CLK_CYCLES}")

    # Print execution time
    stop_time = time.time()
    print(f"\nExecution time of {inspect.currentframe().f_code.co_name}(): {stop_time - start_time:.3f} sec")

    return HOLD_TIME_CLK_CYCLES

def expDecay(t, a, T):
    """
    # exponential decay formula
    """
    return a * np.exp(-t/T)

def reconstDsumLevelFromEdgeAcq(dfIn, holdOffClkCycle, dfColName=None, doPlot=False) -> None:
    """
    # Reconstruct the cumulative sum with the hold-off duration.
    # The measures were taken using "edge" synchronization of the pixels.
    # This reconstruction emulates a "level" synchronization of the pixels.
    # NOTE: changes are directly applied to the DataFrame argument "dfIn"
    # Modify the DataFrame to add a new column: "dsumLevel"
    # dsumLevel: dsum bin reconstructed as if pixels were sampled in 'level' mode.
    # NOTE: with all consecutive samples of dsum (no zero suppression), two methods are available:
    #         1- Convolution of each dsum frame with a kernel vector of length holdOffClkCycle
    #     ->  2- Cumulative sum of each dsum frame subtracted by the same cumulative sum, shifted of holdOffClkCycle
    # holdOffClkCycle -> hold-off duration represented in number of clock cycles.
    #                    holdOffClkCycle supported types:
    #                       int -> single value,
    #                       list/array -> one item for each PDC,
    #                       dictionary -> a key for each PDC and the value to use {0: 50, 1: 51}
    #                       pd.DataFrame -> with columns 'pdcIdx' and dfColName
    # dfColName -> only used if holdOffClkCycle is a DataFrame, the holdOffClkCycle DataFrame column to use
    # doPlot -> For large DataFrames, the execution time is non negligeable.
    #           With doPlot == True, the function plots a graph of the execution time.
    """
    # verify "getDsumBinIdx" has been run on the DataFrame
    if not {"binIdx"}.issubset(dfIn.columns):
        raise KeyError("First call 'getDsumBinIdx' on 'dfIn'"\
                       "to generate column 'binIdx'")
    start_time = time.time()

    # base reconstruction on dataIdxToUse ("binIdx")
    dataIdxToUse = "binIdx"

    # holdend is a temporary variable to keep the expected index of the end of the hold-off period for each bin
    if isinstance(holdOffClkCycle, (int, np.integer)):
        holdEnd = dfIn[dataIdxToUse].astype(np.int16) - holdOffClkCycle
        maxBinToApply = holdOffClkCycle
    elif isinstance(holdOffClkCycle, (list, np.ndarray)):
        pdcList = getPdcIdxInData(dfIn, doPrint=False)
        if len(holdOffClkCycle) != len(pdcList):
            raise ValueError(f"holdOffClkCycle specified as {type(holdOffClkCycle)}. Expecting {len(pdcList)} elements, but got {len(holdOffClkCycle)}")
        lookupDict = {pdcIdx: pdcHold for pdcIdx, pdcHold in enumerate(holdOffClkCycle)}
        seriesHoldOffClkCycle = dfIn["pdcIdx"].map(lookupDict).astype(np.int16)
        holdEnd = dfIn[dataIdxToUse].astype(np.int16) - seriesHoldOffClkCycle
        maxBinToApply = max(holdOffClkCycle)
    elif isinstance(holdOffClkCycle, dict):
        pdcList = getPdcIdxInData(dfIn, doPrint=False)
        if not all(pdcIdx in holdOffClkCycle for pdcIdx in pdcList):
            raise ValueError(f"holdOffClkCycle specified as a dictionary. Expecting a key for each PDC")
        seriesHoldOffClkCycle = dfIn["pdcIdx"].map(holdOffClkCycle).astype(np.int16)
        holdEnd = dfIn[dataIdxToUse].astype(np.int16) - seriesHoldOffClkCycle
        maxBinToApply = max(holdOffClkCycle.values())
    elif isinstance(holdOffClkCycle, pd.DataFrame):
        if dfColName is None:
            raise ValueError(f"When specifying prd as a DataFrame, please specify dfColName")
        if not {"pdcIdx", dfColName}.issubset(holdOffClkCycle.columns):
            raise KeyError(f"Make sure DataFrame holdOffClkCycle contains columns 'pdcIdx' and {holdClkCycles}")
        seriesHoldOffClkCycle = dfIn["pdcIdx"].map(holdOffClkCycle.set_index("pdcIdx")[dfColName]).astype(np.int16)
        holdEnd = dfIn[dataIdxToUse].astype(np.int16) - seriesHoldOffClkCycle
        maxBinToApply = np.max(holdOffClkCycle)
    else:
        raise TypeError(f"Unsupported type {type(holdOffClkCycle)} for holdOffClkCycle")

    # first bin no shift
    dfIn.loc[:, "dsumLevel"] = dfIn["dsum"]

    # init dataframe to shift (explicit copy, not a view or reference to original dfIn)
    dfShift = dfIn[["frameIdx", "pdcIdx", dataIdxToUse, "dsum"]].copy()

    # Dictionary to keep the execution time based on the number of bins to process
    exTime = {}
    exTime[0] = {"nBins":len(dfIn), "time":time.time()}

    # shift for each hold-off bin
    for iBinHold in range(maxBinToApply):
        # print progressing
        print(f"\r{iBinHold+1}/{maxBinToApply} bins", end="")

        # shift vector
        dfShift = dfShift.shift(1, fill_value=0)
        # NOTE: Instead of using grouped DataFrame functions, #
        # do the operation on the complete DataFrame to reduce execution time.
        # To execute the operation, the shifted and original DataFrames must
        # share the same frameIdx and the same pdcIdx.
        # The shifted value of column dataIdxToUse must also be less than the expected end of the hold-off.
        """
        mask = np.logical_and(np.logical_and(dfShift["frameIdx"] == dfIn["frameIdx"],
                                             dfShift["pdcIdx"] == dfIn["pdcIdx"]),
                              holdEnd < dfShift[dataIdxToUse])
        """
        mask = combineFilters([dfShift["frameIdx"] == dfIn["frameIdx"],
                               dfShift["pdcIdx"] == dfIn["pdcIdx"],
                               holdEnd < dfShift[dataIdxToUse]])
        # update the cumulative value of the digital sum for the bins specified by mask
        dfIn.loc[mask, "dsumLevel"] += dfShift.loc[mask, "dsum"]

        # log the execution time
        nBins = mask.sum()
        exTime[iBinHold] = {"nBins":nBins, "time":time.time()}
        if nBins == 0:
            break

    # clearing unrequired variable to same memory
    del mask
    del dfShift
    del holdEnd

    # print execution time
    stop_time = time.time()
    elapsed_time = stop_time - start_time
    print(f"\nExecution time of {inspect.currentframe().f_code.co_name}(): {elapsed_time:.3f} sec")

    # execution time of hold reconstruction, for statistics
    if doPlot:
        nBins, times = zip(*np.array([(nBinData["nBins"], nBinData["time"]) for nBinData in exTime.values()]))
        times = times-times[0]
        title = "execution time as a function of the number of bins to process"
        plt.close(title)
        plt.figure(title)
        plt.plot(nBins[1:]/nBins[0], np.diff(times), '.-')
        plt.xscale("log")
        plt.xlabel("number of bins to process")
        plt.ylabel("execution time (sec)")
        plt.title(title)


# -----------------------------------------------
# --- Plot a frame of digital sum and
# --- compare the different hold-off
# --- reconstruction methods
# -----------------------------------------------
def plotDsumFrameHoldOff(dfFrame, holdOffClkCycle, tauClkCycle):
    """
    # From the extracted frame (dfFrame), plot the digital sum with different
    # hold-off reconstruction methods (convolution and shifted cumulative sum)
    # holdOffClkCycle -> hold-off duration represented in number of clock cycles.
    # tauClkCycle -> scintillator time constant represented in number of clock cycles.
    """
    # Reconstruct the original frame with zero values
    # This cancels the zero suppression applied by hexRead
    binsZ = range(min(dfFrame["binIdx"]), max(dfFrame["binIdx"])+1)
    dsumZ = np.zeros((np.max(dfFrame["binIdx"])+1))
    dsumZ[dfFrame["binIdx"]] = dfFrame["dsum"]

    # use convolution to reproduce the hold-off effect (level acquisition)
    dsumZConv = np.convolve(dsumZ, np.ones((holdOffClkCycle)), mode="full")[:len(dsumZ)]

    # cumulative value of the bins. It corresponds to an infinite hold-off time.
    dsumZCumul = np.cumsum(dsumZ)

    # use shift of cumulative method to reproduce the hold-off effect (level acquisition)
    dsumZCumulH = np.copy(dsumZCumul)
    dsumZCumulH[holdOffClkCycle:] -= dsumZCumulH[:-holdOffClkCycle]

    # extract peak and its bin position
    peak, loc = dfFrame.loc[dfFrame["dtPeak"] == 0, ["dsum", "binIdx"]].values[0]
    print(f"peak:{peak}, loc:{loc}")

    # extract bins
    bins = dfFrame["binIdx"]

    # exponential with parameters found from data
    # NOTE: decay is not a fit, but an exponential decay using peak and loc as parameters
    #       To obtain better results, a fitting method could be used.
    decay = peak*(np.exp(-(bins-loc)/tauClkCycle))
    # remove values larger than peak and smaller than 1
    mask = np.logical_and(decay >= 1, decay <= peak)
    binsDecay = bins[mask]
    decay = decay[mask]
    decayCumul = np.cumsum(decay)

    # plot
    fig = plt.figure(figsize=(12,9))
    ax, axCumul, axLog = fig.subplots(3,1, sharex=True)
    pltColors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    iColorMax = len(pltColors)

    from enum import StrEnum
    class PltColors (StrEnum):
        df = pltColors[0]
        conv = pltColors[1]
        shift = pltColors[2]
        zero = pltColors[4]
        decay = "black"
        cumul = "gray"
        peak = "red"

    # === ax: Number of detected photons per bin === #
    axTitle = "Number of\ndetected photons\nper bin"
    ax.plot(bins, dfFrame["dsum"], 's', markersize=7, color=PltColors.df, label="dsum")
    ax.plot(binsZ, dsumZ, '.', color=PltColors.zero, label="dsum (zero-padded)")
    ax.plot(loc, peak, 'x', color=PltColors.peak, label="peak")
    ax.plot(binsDecay, decay, '-', color=PltColors.decay, label="decay")
    ax.set_title(axTitle.replace('\n', ' '))
    ax.set_ylabel(axTitle)

    # === axCumul: Cumulative number of detected photons per bin === #
    axCumulTitle = "Cumulative number of\ndetected photons\nper bin"
    axCumul.plot(bins, dfFrame["dsumLevel"], 's', markersize=7, color=PltColors.df, label="dsumLevel (df)")
    axCumul.plot(binsZ, dsumZCumul, '-', color=PltColors.cumul, label="dsum (cumulative)")
    axCumul.plot(binsZ, dsumZCumulH, '--', marker=".", markersize=9, color=PltColors.shift, label="dsumLevel (shifted cumulative)")
    axCumul.plot(binsZ, dsumZConv, '.', color=PltColors.conv, label="dsumLevel (convolution)")
    axCumul.plot(binsDecay, decayCumul, '-', color=PltColors.decay, label="decay (cumulative)")
    axCumul.set_title(axCumulTitle.replace('\n', ' '))
    axCumul.set_ylabel(axCumulTitle)

    # === axLog: Logarithm of the number of detected photons per bin === #
    axLogTitle = "Logarithm of the number of\ndetected photons\nper bin"
    dsumLog = np.log(dfFrame["dsum"])
    decayLog = np.log(peak)-(bins-loc)/tauClkCycle
    axLog.plot(bins, dsumLog, 's', markersize=7, color=PltColors.df, label="log(dsum)")
    axLog.plot(bins, np.maximum((decayLog).values, np.zeros_like(decayLog.values)), '.-', color=PltColors.decay, label="log(decay)")
    axLog.set_title(axLogTitle.replace('\n', ' '))
    axLog.set_ylabel(axLogTitle)

    # finalize the plot
    fig.supxlabel('Bin index')
    ax.legend()
    axCumul.legend()
    axLog.legend()
    plt.tight_layout()


# -----------------------------------------------
# --- Analyse the bin distribution
# --- of the digital sum
# -----------------------------------------------
def getNumFrameLongerThanPrd(dfGroup, prd, dfColName=None) -> tuple:
    """
    # extract frames with duration longer or equal to 'prd' period
    # prd -> Period to used to compare the duration of the frames.
    #        Specified in number of clock cycles.
    #        prd supported types:
    #            int -> single value,
    #            list/array -> one item for each PDC,
    #            dictionary -> a key for each PDC and the value to use {0: 50, 1: 51}
    #            pd.DataFrame -> with columns 'pdcIdx' and dfColName
    # dfColName -> only used if prd is a DataFrame, the prd DataFrame column to use
    """
    # verify "getDsumBinIdx" has been run on the DataFrame
    if not {"binIdx"}.issubset(dfGroup.obj.columns):
        raise KeyError("First call 'getDsumBinIdx' on 'dfGroup.obj'"\
                       "to generate column 'binIdx'")

    transformed = dfGroup["binIdx"].transform("last")
    if isinstance(prd, (int, np.integer)):
        frameLenGePrd = transformed.ge(prd)
    elif isinstance(prd, (list, np.ndarray)):
        pdcList = getPdcIdxInData(dfGroup.obj, doPrint=False)
        if len(prd) != len(pdcList):
            raise ValueError(f"prd specified as {type(prd)}. Expecting {len(pdcList)} elements, but got {len(prd)}")
        lookupDict = {pdcIdx: pdcPrd for pdcIdx, pdcPrd in enumerate(prd)}
        frameLenGePrd = transformed >= dfGroup.obj["pdcIdx"].map(lookupDict)
    elif isinstance(prd, dict):
        pdcList = getPdcIdxInData(dfIn, doPrint=False)
        if not all(pdcIdx in prd for pdcIdx in pdcList):
            raise ValueError(f"prd specified as a dictionary. Expecting a key for each PDC")
        frameLenGePrd = transformed >= dfGroup.obj["pdcIdx"].map(prd)
    elif isinstance(prd, pd.DataFrame):
        if dfColName is None:
            raise ValueError(f"When specifying prd as a DataFrame, please specify dfColName")
        if not {"pdcIdx", dfColName}.issubset(prd.columns):
            raise KeyError(f"Make sure DataFrame prd contains columns 'pdcIdx' and '{dfColName}'")
        seriesPrd = dfGroup.obj["pdcIdx"].map(prd.set_index("pdcIdx")[dfColName])
        frameLenGePrd = transformed >= seriesPrd
    else:
        raise TypeError(f"Unsupported type {type(prd)} for prd")

    nKeep, nTotal, pctKeep, pctDrop = pdh.filterStats(frameLenGePrd,
                                                      description=inspect.currentframe().f_code.co_name)
    return frameLenGePrd, (nKeep, nTotal, pctKeep, pctDrop)

def getNumFrameDtLargerThanTh(dfGroup, thGt=1, dfColName=None) -> tuple:
    """
    # count the number of frames with non consecutive bins
    # thGt -> Threshold of 'dt' to be greater than.
    #        Specified in number of clock cycles.
    #        thGt supported types:
    #            int -> single value,
    #            list/array -> one item for each PDC,
    #            dictionary -> a key for each PDC and the value to use {0: 50, 1: 51}
    #            pd.DataFrame -> with columns 'pdcIdx' and dfColName
    # dfColName -> only used if thGt is a DataFrame, the thGt DataFrame column to use
    """
    # verify "getDsumDt" has been run on the DataFrame
    if not {"dt"}.issubset(dfGroup.obj.columns):
        raise KeyError("First call 'getDsumDt' on 'dfGroup.obj'"\
                       "to generate column 'dt'")

    transformed = dfGroup["binIdx"].transform("last")
    if isinstance(thGt, (int, np.integer)):
        dtMaxFrame = transformed.gt(thGt)
    elif isinstance(thGt, (list, np.ndarray)):
        pdcList = getPdcIdxInData(dfGroup.obj, doPrint=False)
        if len(thGt) != len(pdcList):
            raise ValueError(f"thGt specified as {type(thGt)}. Expecting {len(pdcList)} elements, but got {len(thGt)}")
        lookupDict = {pdcIdx: pdcPrd for pdcIdx, pdcPrd in enumerate(thGt)}
        dtMaxFrame = transformed >= dfGroup.obj["pdcIdx"].map(lookupDict)
    elif isinstance(thGt, dict):
        pdcList = getPdcIdxInData(dfIn, doPrint=False)
        if not all(pdcIdx in thGt for pdcIdx in pdcList):
            raise ValueError(f"thGt specified as a dictionary. Expecting a key for each PDC")
        dtMaxFrame = transformed >= dfGroup.obj["pdcIdx"].map(thGt)
    elif isinstance(thGt, pd.DataFrame):
        if dfColName is None:
            raise ValueError(f"When specifying prd as a DataFrame, please specify dfColName")
        if not {"pdcIdx", dfColName}.issubset(thGt.columns):
            raise KeyError(f"Make sure DataFrame thGt contains columns 'pdcIdx' and '{dfColName}'")
        seriesPrd = dfGroup.obj["pdcIdx"].map(thGt.set_index("pdcIdx")[dfColName])
        dtMaxFrame = transformed >= seriesPrd
    else:
        raise TypeError(f"Unsupported type {type(thGt)} for thGt")

    nKeep, nTotal, pctKeep, pctDrop = pdh.filterStats(dtMaxFrame,
                                                      description=inspect.currentframe().f_code.co_name)
    return dtMaxFrame, (nKeep, nTotal, pctKeep, pctDrop)

def getNumDtPerFrame(dfGroup) -> None:
    """
    # for each frame, look at the number of different dt values.
    # print statistics on the number of frames with N different dt values.
    """
    # verify "getDsumDt" has been run on the DataFrame
    if not {"dt"}.issubset(dfGroup.obj.columns):
        raise KeyError("First call 'getDsumDt' on 'dfGroup.obj'"\
                       "to generate column 'dt'")

    # for each frame, extract the number of different dt values
    frameNDt = dfGroup["dt"].nunique()

    # for a small range of number of dt values, print the stats
    for thDt in range(1, 20):
        maskNDt = frameNDt["dt"] == thDt
        nKeep, _, _, _ = pdh.filterStats(maskNDt, description=f"frames with {thDt} dt")
        if nKeep == 0:
            # no more data, end loop
            break

def getFrameWithPeakNotAtFirstBin(dfGroup) -> tuple:
    """
    # frames on which the peak is not the first bin
    """
    # verify "getDsumFramePeakInfo" has been run on the DataFrame
    if not {"dtPeak"}.issubset(dfGroup.obj.columns):
        raise KeyError("First call 'getDsumFramePeakInfo' on 'dfGroup'"\
                       "to generate column 'dtPeak'")
    framesPeakNotFirstBin = dfGroup["dtPeak"].transform("min").le(-1)
    nKeep, nTotal, pctKeep, pctDrop = pdh.filterStats(framesPeakNotFirstBin, description=f"frames with peak not at first bin")
    return framesPeakNotFirstBin, (nKeep, nTotal, pctKeep, pctDrop)

def getMaxNumberOfPixelTriggered(dfIn, numPixEnabled:list):
    """
    # Extract the maximum number of pixels triggered at the same time.
    # Compare the results with the setting of the number of enabled pixels.
    # NOTE: changes are directly applied to the DataFrame argument "dfIn"
    # Modify the DataFrame to add new columns: "numPix" and "nAvail"
    # numPix: Total number of pixel per PDC.
    # nAvail: after each bin, calculate the number of available pixels,
    #         considering some pixels are in hold-off.
    # NOTE: for some frames, it is possible that the error on the hold-off estimation
    #       causes nAvail to go to zero even if dsum is not.
    # numPixEnabled -> Default setting of the number of enabled pixels.
    #                  The number of elements in numPixEnabled must match the number
    #                  of PDCs in the DataFrame (dfIn).
    """
    # verify "reconstDsumLevelFromEdgeAcq" has been run on the DataFrame
    if not {"dsumLevel"}.issubset(dfIn.columns):
        raise KeyError("First call 'reconstDsumLevelFromEdgeAcq' on 'dfIn'"\
                       "to generate column 'dsumLevel'")
    # Find maximum number of pixel that triggered in the DataFrame.
    # NOTE: Can be used as an estimate if number of enabled pixels is unknown.
    maxNumOfPixelTriggered = dfIn.groupby("pdcIdx")["dsumLevel"].max()

    for iPdc, (maxNumPixTrg, numPixEnPerPdc) in enumerate(zip(maxNumOfPixelTriggered, numPixEnabled)):
        print(f"PDC{iPdc}, max:{maxNumPixTrg}, en:{numPixEnPerPdc}, ({100.0*maxNumPixTrg/numPixEnPerPdc:.03f} %)")

    # for each sample, store the total number of enabled pixel
    dfIn["numPix"] = dfIn["pdcIdx"].map(pd.Series(numPixEnabled)).astype(np.uint16)

    # for each sample, store the number of remaining pixel after each bin
    # NOTE: casted numPix as an integer to prevent negative unsigned substraction to end up as large numbers
    # NOTE: added clip to keep value between 0 and numPix
    dfIn["nAvail"] = np.clip(dfIn["numPix"].astype(int)-dfIn["dsumLevel"]+dfIn["dsum"], 0, dfIn["numPix"]).astype(np.uint16)

    # Verify "nAvail" is valid
    criteria = 'dfIn["dsum"] > dfIn["nAvail"]'
    numBinCriteria = eval(criteria).sum()
    print(f'num of bins with {criteria}: {numBinCriteria}')


# -----------------------------------------------
# --- Linearity correction
# -----------------------------------------------
def applyDsumLinearity(dfIn) -> tuple:
    """
    # Apply linearity correction on "dsum" based on the number of pixel trigger
    # and the number of available pixels.
    # NOTE: changes are directly applied to the DataFrame argument "dfIn"
    # Modify the DataFrame to add a new column: "dsumLin"
    # dsumLin: Digital sum bins with applied linearity.
    """
    # verify "getMaxNumberOfPixelTriggered" has been run on the DataFrame
    if not {"nAvail"}.issubset(dfIn.columns):
        raise KeyError("First call 'getMaxNumberOfPixelTriggered' on 'dfIn'"\
                       "to generate column 'nAvail'")
    dfIn["dsumLin"] = dfIn["dsum"].astype(np.float32)

    # find where it is safe to apply linearity (no zero division, no log(0)
    maskLin = combineFilters([dfIn["nAvail"] > 0, dfIn["dsum"] <= dfIn["nAvail"]])


    N = dfIn.loc[maskLin, "nAvail"].max().astype(np.float32) # Number of SPADs

    logTerm = 1-dfIn.loc[maskLin, "dsum"].values/dfIn.loc[maskLin, "nAvail"].values
    dfIn.loc[maskLin, "dsumLin"] = -N*np.log(logTerm)

    # evaluate the number of bins affected by linearity
    maskLin = dfIn["dsum"] != dfIn["dsumLin"]
    nKeep, nTotal, pctKeep, pctDrop = pdh.filterStats(maskLin, description=f"bins with 'dsum' different after linearization")
    return maskLin, (nKeep, nTotal, pctKeep, pctDrop)


# -----------------------------------------------
# --- Pulse shape discrimination
# -----------------------------------------------
def applyPsd(dfIn, nPrompt_l=0, nPrompt_r=1, nTotal=-1, colToPsd="dsum") -> pd.DataFrame:
    """
    # Apply PSD on a DataFrame
    # NOTE: changes are directly applied to the DataFrame argument "dfIn"
    # Modify the DataFrame to add a new column: "psd"
    # psd: Result of the PSD calulation
    # dfIn -> DataFrame to apply PSD on.
    # nPrompt -> number of "prompt" bins to keep starting from peak value
    # nTotal -> number of "total" bins to keep starting from peak value.
    #           Set nTotal to -1 to use all data in the frame.
    # colToPsd -> column of the DataFrame to use for PSD calculation
    # Return post-processed CSV with Energy/PSD extracted for each event of each PDCs
    """
    # remove columns from a previous run
    dfIn.drop(columns=["psdPrompt", "psdTotal", "psd"], errors='ignore', inplace=True)

    # select which rows to use as prompt
    dfIn.loc[:, "psdPrompt"] = dfIn[colToPsd].where(dfIn["dtPeak"].between(nPrompt_l, nPrompt_r), 0)

    # select which rows to use as total
    if nTotal == -1:
        # using all data
        dfIn.loc[:, "psdTotal"] = dfIn[colToPsd]
    else:
        # using a subset of data
        dfIn.loc[:, "psdTotal"] = dfIn[colToPsd].where(dfIn["dtPeak"].between(nPrompt_l, nTotal), 0)

    # groupby data
    dfGroup = groupDfByFrameAndPdc(dfIn)

    # applying PSD on the selected data
    dfIn.loc[:, "psd"] = dfGroup["psdPrompt"].transform("sum")/dfGroup["psdTotal"].transform("sum")
    
    
    # Get the list of PDCs in the DataFrame
    pdcList = getPdcIdxInData(dfIn)
    df_psd = pd.DataFrame()
    for iAx, iPdc in enumerate(pdcList):
        dfPdc = dfIn.loc[dfIn["pdcIdx"] == iPdc, :]
        dfPdcGroup = dfPdc.groupby("frameIdx")
        
        x = dfPdcGroup["energy"].nth(0).values
        y = dfPdcGroup["psd"].nth(0).values
        
        df_psd = pd.concat([df_psd, pd.DataFrame({"Energy [a.u.]": x, "PSD":y, "iPDC":iPdc})], ignore_index=True)
        
    
    # remove temporary columns
    dfIn.drop(columns=["psdPrompt", "psdTotal"], errors='ignore', inplace=True)

    return df_psd


def plotPsdFctEnergy(dfIn):
    """
    # Generate a histogram of the PSD as a function of the energy
    # NOTE: this function requires "applyPsd" to be run on the DataFrame
    """
    # verify "getDsumPrelimEnergy" has been run on the DataFrame
    if not {"energy"}.issubset(dfIn.columns):
        raise KeyError("First call 'getDsumPrelimEnergy' on 'dfIn' to generate column 'energy'")
    # verify "applyPsd" has been run on the DataFrame
    if not {"psd"}.issubset(dfIn.columns):
        raise KeyError("First call 'applyPsd' on 'dfIn to generate column 'psd'")

    # Get the list of PDCs in the DataFrame
    pdcList = getPdcIdxInData(dfIn)

    # prepare a subplot for each PDC
    fig, axes = plt.subplots(nrows=len(pdcList), ncols=1, figsize=(12, 9),
                             squeeze=False, sharex=False)
    axes = axes.flatten()

    for iAx, iPdc in enumerate(pdcList):
        dfPdc = dfIn.loc[dfIn["pdcIdx"] == iPdc, :]
        dfPdcGroup = dfPdc.groupby("frameIdx")
        x = dfPdcGroup["energy"].nth(0).values
        y = dfPdcGroup["psd"].nth(0).values

        xmax = np.max(x)
        xmin = np.min(x)
        ymax = np.max(y)
        ymin = np.min(y)
        print(f"PDC{iPdc} -> xmin:{xmin:.03f}, xmax:{xmax:.03f}, ymin:{ymin:.03f}, ymax:{ymax:.03f}")

        xBins = np.arange(xmin*0.8, xmax+1)
        yBins = np.arange(ymin*0.8, ymax+0.01, 0.005)
        H, xedges, yedges = np.histogram2d( x=x, y=y,
                                            bins=(xBins, yBins))

        # show the histogram
        df_psd = pd.DataFrame({     "Energy [a.u.]": x, "PSD":y})
        g = sns.JointGrid(df_psd, x="Energy [a.u.]", y= "PSD", marginal_ticks=True)
        g.plot_joint(
            sns.histplot, discrete=(False, False),
            cmap="gist_rainbow", pmax=.8, cbar=False
        )
        g.plot_marginals(sns.histplot, element="step", color="#03012d")

        axes[iAx].imshow(H.T, cmap="gist_rainbow", origin="lower", aspect="auto",
                        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                        norm=colors.LogNorm(vmin=max(np.min(H), 1e-6), vmax=np.max(H)), zorder=0)
        axes[iAx].set_xlabel("energy (number of photons detected)")
        axes[iAx].set_ylabel(f"PDC{iPdc}\nPSD")

    
    plt.tight_layout()
    return fig, axes


# -----------------------------------------------
# --- Energy Spectrum
# -----------------------------------------------
def getDsumEnergy(dfIn, col="dsum", transform="sum", dtype=np.uint16):
    """
    # Calculate the energy based on the selected method (col, transform)
    # NOTE: changes are directly applied to the DataFrame argument "dfIn"
    # Modify the DataFrame to add a new column: "energy"
    # energy: calculated energy for each PDC and each frame.
    # col -> column to use to produce the energy spectra
    #        Recommended options: "dsum", "dsumLin", "dsumLevel", "energy"
    # transform -> DataFrame group transformation to apply on 'col'.
    #              Recommended options: "sum", "max", "mean"
    #              Can be a custom user function or a supported default tranform function.
    # E.g.:
    #     col="dsum", transform="sum" -> will sum all bins of the same group
    #     col="dsumLin", transform="sum" -> will sum all linearized bins of the same group
    #     col="dsumLevel", transform="max" -> will use the maximum value
    #                                         of the level reconstructed dsum
    #     col="energy", transform="mean" -> will use the average of the
    #                                       previously calculated energy column
    # NOTE: to keep multiple methods of energy calculation in the DataFrame,
    #       call this function, then rename the "energy" column,
    #       then call the function with different parameters.
    """
    dfIn.loc[:, "energy"] = groupDfByFrameAndPdc(dfIn)[col].transform(transform).astype(dtype)


def getEnergySpectra(dfIn, combinedFilter=None) -> pd.DataFrame:
    """
    # combinedFilter -> filter to select a specific subset of the DataFrame to use
    """
    # Get the list of PDCs in the DataFrame
    pdcList = getPdcIdxInData(dfIn, doPrint=False)

    #dfSpectra = pd.DataFrame(columns=["pdcIdx", "binCenters", "binCounts"])
    dfSpectra = pd.DataFrame()
    for iPdc in pdcList:
        try:
            pdcFilter = dfIn["pdcIdx"] == iPdc
            if combinedFilter is None:
                # only use PDC as filter
                pdcCombinedFilter = pdcFilter
            else:
                pdcCombinedFilter = combineFilters([pdcFilter, combinedFilter])

            label = f"PDC{iPdc}"
            _ = pdh.filterStats(pdcCombinedFilter, description=label)
            pdcData = dfIn.loc[pdcCombinedFilter].groupby("frameIdx")["energy"].transform("first")

        except KeyError as ex:
            print(f"Exception: {ex}")
            continue

        # generate the energy spectrum
        bins = np.arange(np.min(pdcData), np.max(pdcData))
        binCounts, binEdges = np.histogram(pdcData, bins=bins)
        binCenters = (binEdges[1:] + binEdges[:-1])/2.0
        # DataFrame with a single spectrum (one PDC)
        dfSpectrum = pd.DataFrame({"pdcIdx": [iPdc]*len(binCounts),
                                   "binCenters": binCenters,
                                   "binCounts": binCounts})
        # DataFrame with a spectrum for each PDC
        dfSpectra = pd.concat([dfSpectra, dfSpectrum], ignore_index=True)

    return dfSpectra

def extractPdcSpectrum(dfSpectra, iPdc) -> pd.DataFrame:
    """
    # get the spectrum data of a single PDC
    # dfSpectra: DataFrame with spectrum information from all PDCs.
    #            Generated using getEnergySpectra
    # iPdc -> index of the PDC to extract.
    """
    return dfSpectra.loc[dfSpectra["pdcIdx"] == iPdc, :].copy(deep=False)

def findPhotoPeaks(dfSpectrum, spectrumBinMin=0, spectrumBinMax=-1,
                   nBinPerBin=10,
                   relThreshold=0.08,
                   peaksTolerance=80) -> dict:
    """
    # fit peaks of an energy spectrum
    # NOTE: preliminary. Results are not garantied.
    # NOTE: Multiple parameters are hardcoded and not made available as a function parameter.
    #       User can tweak the values if required depending on the measurements.
    # dfSpectrum -> DataFrame with the spectrum data of a single PDC
    # spectrumBinMin -> miminum bin to consider as relevant data. Below is considered as noise.
    # spectrumBinMax -> maximum bin to consider as relevant data. Above is considered as noise.
    # nBinPerBin -> to help the algorithm to find peaks, histograms are rebinned.
    #               number of original bins to regroup into the rebinned histogram.
    # relThreshold -> relative threshold to the photoPeak (value from 0 to 1)
    # peakTolerance -> maximum distance between peaks to remove peaks
    #                  It prevents multiple peaks too close to each other.
    #
    # The function returns a dictionary with the following items:
    #   "photoPeak": the maximum value in the spectrum associated with the photodetection range.
    #   "noisePeak": the maximum value in the spectrum associated with the noise range (low energy).
    #   "centerReBin": the bin center of the rebinned energy spectrum
    #   "countsReBin": the bin counts of the rebinned energy spectrum
    #   "dfPeaks": DataFrame with the initial parameter for the fit ("gain", "center", "sigma")
    #   "dfParams": DataFrame with the information of the found peaks
    """
    outputDict = {}
    try:
        # print PDC idx
        print(f"=== PDC{dfSpectrum['pdcIdx'].unique()[0]} ===")

        # extract as two numpy variables
        binCenters = dfSpectrum["binCenters"].values
        binCounts = dfSpectrum["binCounts"].values.astype(np.float32)

         # keep noise peak bins in graph, but adjust max based on photopeak
        photoPeak = np.max(binCounts[spectrumBinMin:spectrumBinMax])
        noisePeak = np.max(binCounts[:spectrumBinMin])
        outputDict["photoPeak"] = photoPeak
        outputDict["noisePeak"] = noisePeak

        # rebin the histogram to remove noise, but keep tendency
        centerReBin, countsReBin = esah.rebinHisto(binCenters, binCounts, nBinPerBin=nBinPerBin)
        outputDict["centerReBin"] = centerReBin
        outputDict["countsReBin"] = countsReBin

        # find peaks (from scipy find_peaks) in re binned histogram
        peaksIdx, props = find_peaks(countsReBin, prominence=photoPeak*relThreshold*0.20, wlen=50)

        # create a DataFrame with the peaks
        dfPeaks = pd.DataFrame()
        dfPeaks["gain"] = countsReBin[peaksIdx]
        dfPeaks["center"] = centerReBin[peaksIdx]
        dfPeaks["sigma"] = 50.0 # initial value (ignored)

        # remove peaks that are too close (within peaksTolerance)
        dfPeaks = esah.cleanPeaks(dfPeaks, tolerance=peaksTolerance)
        print(f"dfPeaks:\n{dfPeaks}")

        # remove peaks with gain below or above thresholds
        peakLowGainThreshold = photoPeak*0.05
        peakHighGainThreshold = np.max(binCounts)
        gainKeepFilter = combineFilters([dfPeaks["gain"] > peakLowGainThreshold,
                                         dfPeaks["gain"] < peakHighGainThreshold])
        dfPeaks = dfPeaks.loc[gainKeepFilter, :].reset_index(drop=True)
        #print(f"dfPeaks:\n{dfPeaks}")

        outputDict["dfPeaks"] = dfPeaks

        # boundaries for the fit
        dfBoundsUp = dfPeaks.copy()
        dfBoundsDn = dfPeaks.copy()
        dfBoundsUp["sigma"] = 250.0
        dfBoundsDn["sigma"] = 0
        dfBoundsUp["center"] += peaksTolerance
        dfBoundsDn["center"] -= peaksTolerance
        dfBoundsDn.loc[dfBoundsDn["center"] < 1, "center"] = 1
        dfBoundsUp["gain"] = peakHighGainThreshold
        dfBoundsDn["gain"] = peakLowGainThreshold
        #print(f"up:\n{dfBoundsUp}")
        #print(f"dn:\n{dfBoundsDn}")

        p0 = dfPeaks.to_numpy().flatten()
        bUp = dfBoundsUp.to_numpy().flatten()
        bDn = dfBoundsDn.to_numpy().flatten()
        if esah.fitY0:
            p0 = [0] + p0
            bUp = [np.inf] + bUp
            bDn = [-np.inf] + bDn

        # apply fit (curve_fit from scipy)
        popt, pcov, infodict, errmsg, ier = curve_fit(esah.multiGaussian,
                                                      binCenters, binCounts,
                                                      p0=p0,
                                                      bounds=(bDn, bUp),
                                                      full_output=True)
        #print(f"infodict:\n{infodict}")
        #print(f"errmsg:\n{errmsg}")
        #print(f"ier:\n{ier}")

        # format fit results
        paramsName = ("gain", "center", "sigma")
        perr = np.sqrt(np.diag(pcov))
        dfParams = esah.gausParams2df(popt, paramsName)
        dfPerr = esah.gausParams2df(perr, [param+"_err" for param in paramsName])
        dfParams = dfParams.join(dfPerr)
        dfPerrPct = esah.gausParams2df(100.0*perr/popt, [param+"_errPct" for param in paramsName])
        dfParams = dfParams.join(dfPerrPct)

        # remove gaussian with contribution too small
        #dfParams = dfParams.loc[dfParams["gain"] > 0.05*dfParams["gain"].max(), :]

        # remove gaussian with error on gain larger than gain
        #dfParams = dfParams.loc[dfParams["gain_err"] < dfParams["gain"], :]

        # get each fit to extract the FWHM value
        for index, g in dfParams.iterrows():
            gfit = esah.gaus(x=binCenters, a=g["gain"], x0=g["center"], sigma=g["sigma"])
            dfParams.loc[index, "fwhm"] = esah.get_fwhm(binCenters, gfit)

        # overall fit to visualize the match
        yFit = esah.multiGaussian(binCenters, *dfParams[["gain", "center", "sigma"]].to_numpy().flatten())
        dfSpectrum.loc[:, "binCountsFit"] = yFit

        # energy resolution
        dfParams["ER"] = 100.0*dfParams["fwhm"]/dfParams["center"]

        # print results
        print(f"dfParams:\n{dfParams}")

        outputDict["dfParams"] = dfParams

    except RuntimeError as ex:
        print(f"    RuntimeError {ex}")

    return outputDict

def _plotOnPick(event):
    """
    # This callback function is used to update a descriptive text in the plot
    # to show peak parameters.
    # Used in plotEnergySpectra
    """
    # Retrieve data from the picked point
    artist = event.artist

    # Retrieve figure and ax
    fig = artist.figure
    ax = artist.axes#plt.gca()

    # Retrieve index of data
    ind = event.ind[0]

    # Retrieve clicked x and y
    x, y = np.array(artist.get_data())[:, ind]

    # Find the corresponding peak
    dx = (ax.dfParams["center"] - x).abs()
    peakIdx = dx.idxmin()
    peak = ax.dfParams.loc[peakIdx]

    # Update text and redraw the figure
    ax.textObj.set_text(f"Peak at {peak['center']:.01f} -> FWHM: {peak['fwhm']:.01f}, ER: {peak['ER']:.01f} %")
    fig.canvas.draw()


def plotEnergySpectra(dfSpectra, axes=None, peakInfos=None,
                      plotRebin=False, plotPeaks=True,
                      plotFit=True, plotEachPeakFit=True):
    """
    # plot a spectrum for each PDC in dfSpectra.
    # dfSpectra -> generated by function getEnergySpectra.
    #              Contains columns: "pdcIdx", "binCenters", "binCounts", ["binCountsFit"]
    # axes -> to plot to an existing plot, specify axes as a list (one item per PDC in dfSpectra)
    # peakInfos -> to add elements to the plots, specify peakInfos from function findPhotoPeaks.
    #              Can be a dictionary (for single PDC) or a list of dictionaries.
    # plotRebin -> only valid if peakInfos is specified. Add rebinned histogram.
    # plotPeaks ->`only valid if peakInfos is specified. Add peaks in plot.
    # plotFit -> only valid if peakInfos is specified. Add fit curve to plot
    # plotEachPeakFit -> only valid if peakInfos is specified. Add a fit curve for each peak.
    """
    pdcList = getPdcIdxInData(dfSpectra, doPrint=False)

    if axes is not None and len(pdcList) != len(axes):
        raise ValueError("Expected peakInfos size to match the number of PDCs in dfSpectra")

    if axes is None:
        # open a figure with an ax for each PDC
        fig, axes = plt.subplots(nrows=len(pdcList), ncols=1,
                                figsize=(12, 9), squeeze=False)
        axes = axes.flatten()
    else:
        fig = axes[0].get_figure()

    for idx, iPdc in enumerate(pdcList):
        dfSpectrum = extractPdcSpectrum(dfSpectra, iPdc)
        axes[idx].step(dfSpectrum["binCenters"], dfSpectrum["binCounts"],
                       label=f"PDC{iPdc}")

        if peakInfos is not None:
            try:
                if isinstance(peakInfos, list):
                    peakInfoDict = peakInfos[idx]
                elif isinstance(peakInfos, dict):
                    peakInfoDict = peakInfos
                else:
                    # will trigger an error
                    continue


                # if user wants to plot rebinned histogram
                if plotRebin:
                    axes[idx].step(peakInfoDict["centerReBin"],
                                   peakInfoDict["countsReBin"],
                                   label="rebinned")

                # if user wants to plot peaks
                if plotPeaks:
                    dfParams = peakInfoDict["dfParams"]
                    axes[idx].plot( dfParams["center"],
                                    dfParams["gain"],
                                    "s", color="red",
                                    zorder=100,
                                    label="peaks",
                                    picker=True)
                    axes[idx].textObj = axes[idx].text(0.10, 0.90, "", transform=axes[idx].transAxes)
                    axes[idx].dfParams = dfParams
                    if not hasattr(fig, "mpl_connect"):
                        fig.mpl_connect = fig.canvas.mpl_connect("pick_event", _plotOnPick)

                # if user wants to plot each peak as a gaussian
                if plotEachPeakFit:
                    for iPeak, peakDict in enumerate(dfParams.to_dict(orient="records")):
                        if iPeak == 0:
                            extraArgs = {"label": "single peak fit"}
                        else:
                            extraArgs = {}
                        peakFit = esah.gaus(x=dfSpectrum["binCenters"],
                                            a=peakDict["gain"],
                                            x0=peakDict["center"],
                                            sigma=peakDict["sigma"])
                        axes[idx].step( dfSpectrum["binCenters"], peakFit,
                                        "--", color="gray", **extraArgs)

                # if user wants to plot the fit
                if plotFit and "binCountsFit" in dfSpectrum.columns:
                    axes[idx].step( dfSpectrum["binCenters"], dfSpectrum["binCountsFit"],
                                    "-", color="black",
                                    label=f"fit")

                # adjust the y limits to prevent noise peak to set the range
                if "photoPeak" in peakInfoDict:
                    axes[idx].set_ylim((0, 1.2*peakInfoDict["photoPeak"]))

                # add legend
                axes[idx].legend()

            except (IndexError, KeyError) as ex:
                print(ex)
                continue


# -----------------------------------------------
# --- Methods to select a frame
# --- or a sub range of frames
# -----------------------------------------------
def getFramesByPdcIdxAndFrameIdx(dfIn, method:str, numFrames:int=None,
                                 pdcIdx:list=None, frameIdx=None) -> pd.DataFrame:
    """
    # Function to return a subset of a DataFrame (dfIn), based on PDC index, or frame index.
    # method -> method used to select the subset. Supported values are :
    #   "first": use the first "numFrames" items of the DataFrame.
    #   "last": use the last "numFrames" items of the DataFrame.
    #   "random": randomly select "numFrames" items of the DataFrame.
    #   "loc": use the .loc function (label) of the DataFrame to keep index specified with "frameIdx".
    #          WARNING: frameIdx indexes must be present in the DataFrame.
    #   "iloc": use the .iloc function (integer location) of the DataFrame to keep index specified with "frameIdx".
    #           NOTE: frameIdx indexes are the locations of the DataFrame.
    # numFrames -> with methods: first, last and random, specify the number of frame to keep.
    # pdcIdx -> specify the index of the PDCs to keep.
    #           If None, any PDC can be used.
    #           If "all", all PDC are used.
    # frameIdx -> with methods: loc and iloc, specify the frames labels or position to keep.
    """
    supportedMethods = ["first", "last", "random", "loc", "iloc"]
    if method not in supportedMethods:
        raise ValueError(f"Supported methods are: {supportedMethods}")
    if method in supportedMethods[-2:-1]:
        if frameIdx is None:
            raise ValueError(f"While using method '{method}', 'frameIdx' must be set")
        if numFrames is not None:
            print(f"{fgColors.yellow}WARNING: with method '{method}', 'numFrames' is ignored.{fgColors.endc}")
    else:
        if numFrames is None:
            print(f"{fgColors.yellow}WARNING: with method '{method}', 'numFrames' is required. Using 1 frame.{fgColors.endc}")
            numFrames = 1

    if pdcIdx is not None:
        if isinstance(pdcIdx, str) and pdcIdx.casefold() == "all":
            pdcIdx = getPdcIdxInData(dfIn, doPrint=False)
        elif not hasattr(pdcIdx, "__iter__"):
            pdcIdx = [pdcIdx]
    else:
        # resulting DataFrame will have frames with any PDC index
        pdcIdx = [-1]

    filters = []
    for iPdc in pdcIdx:
        if iPdc == -1:
            frames = dfIn[["pdcIdx", "frameIdx"]].drop_duplicates()
        else:
            frames = dfIn.loc[dfIn["pdcIdx"] == iPdc, ["pdcIdx", "frameIdx"]].drop_duplicates()

        if method == "first":
            framesToKeep = frames.iloc[:numFrames, :]
        elif method == "last":
            framesToKeep = frames.iloc[-numFrames:, :]
        elif method == "random":
            framesToKeep = frames.sample(n=numFrames)
        elif method == "loc":
            try:
                framesToKeep = frames.loc[frameIdx, :]
            except KeyError:
                print(f"{fgColors.red}Specified 'frameIdx' ({frameIdx}) is not found for the specified PDC.{fgColors.endc}")
                continue
        elif method == "iloc":
            framesToKeep = frames.iloc[frameIdx, :]

        for i, (pdc, frame) in framesToKeep.iterrows():
            filters.append({"pdcIdx": pdc, "frameIdx": frame})

    # transform list into DataFrame for next operation
    dfFilter = pd.DataFrame(filters, columns=["pdcIdx", "frameIdx"])

    # filter using fast "merge" method. Keeping only rows that are
    # in both DataFrames.
    dfOut = dfIn.merge(dfFilter, on=["pdcIdx", "frameIdx"], how='inner')
    return dfOut


# -----------------------------------------------
# --- Function to plot a histogram
# --- of a given column of the DataFrame
# -----------------------------------------------
def plotDfColHisto(dfIn, col) -> None:
    """
    # Plot a histogram from a DataFrame with the data at the specified column.
    # col -> name of the column to produce a histogram on.
    """
    print(f"Processing {col}")
    numberOfBins = abs(max(dfIn[col])-min(dfIn[col]))

    # limit the number of bins in the histogram
    maxNumBins = 256
    if numberOfBins > maxNumBins:
        print(f"  numberOfBins ({numberOfBins}) is too large, reducing to {maxNumBins}")
        numberOfBins = maxNumBins

    title = f"Histogram of {col}"
    plt.close(title)
    plt.figure(title)
    plt.hist(dfIn[col], bins=numberOfBins, histtype="bar",
             density=True, cumulative=False)
    plt.yscale("log")
    plt.xlabel(f"{col}")
    plt.ylabel("counts")
    plt.title(title)



# -----------------------------------------------
# --- Function to plot frames from a subset
# --- of a DataFrame
# -----------------------------------------------
def plotDsumFrames(dfFramesGroup,
                   xCol="binIdx",
                   yColsToPlot:list=["dsum"]) -> None:
    """
    # Plot a subset of a DataFrame.
    # xCol -> name of the column to use on the X axis
    # yColsToPlots -> list of names of the columns to use on the Y axis
    """
    # store the list of PDC present in the sub DataFrame
    pdcIdxFound = getPdcIdxInData(dfFramesGroup.obj, doPrint=False)

    # create a figure for each PDC
    figDict = {}
    for iPdc in pdcIdxFound:
        #fig, axes = plt.subplots(nrows=len(yColsToPlot), ncols=1, figsize=(12, 9), layout='constrained')
        fig, axes = plt.subplots(nrows=len(yColsToPlot), ncols=1, figsize=(12, 9), sharex=True)
        # In case there is only one PDC, axes will no be iterable. Make it iterable for generale purpose cases.
        if isinstance(axes, matplotlib.axes._axes.Axes):
            axes = [axes]
        figDict[iPdc] = {"fig":fig, "axes": axes}

        # Settings for each subplot
        for i, ax in enumerate(axes):
            ax.set_ylabel(yColsToPlot[i], fontsize=14)
            ax.tick_params(labelsize=12)
            ax.tick_params(labelsize=12)
            #ax.legend()
        fig.supxlabel(xCol)
        fig.suptitle(f"PDC{iPdc}")

    # Keep the index of frame for each PDC to alternate colors in each subPlot
    iPdcFrame = np.array([0]*len(pdcIdxFound))

    # loop for each frame of the group
    for i, (group, dfFrame) in enumerate(dfFramesGroup):
        #print(f"group:{group}")
        try:
            (iFrame, iPdc, peakBin) = group
        except Exception as ex:
            (iFrame, iPdc) = group

        # plot results
        iAx = pdcIdxFound.index(iPdc)

        # select the ax to plot to
        #ax = figDict[iPdc]["axes"][iAx]

        # index of the trace for each PDC
        index = iPdcFrame[iAx]

        for i, col in enumerate(yColsToPlot):
            ax = figDict[iPdc]["axes"][i]
            ax.plot(dfFrame[xCol], dfFrame[col],
                    marker="D", linestyle='None',
                    markersize=5,
                    label=col,
                    color=f"C{index%5}")

        # increment for next trace
        iPdcFrame[iAx] += 1

    for fig, axes in (pdcDict.values() for pdcDict in figDict.values()):
        fig.tight_layout()





# -----------------------------------------------
# --- testing functions
# -----------------------------------------------
if __name__ == "__main__":
    # NOTE: This is an example code using all the provided analysis functions.
    # User can run the module as is, or import it into a new file and use only
    # the required function to reduce execution time on large DataFrames.

    # Look for an environment variable pointing to a CSV file generated by hexRead
    datafile = os.environ.get("CSV_FILE", None)
    try:
        if datafile is None:
            datafile = input("Please specify here a '.csv' file to use for the analyse or\nuse 'CSV_FILE' environment variable:\n")

        # make sure the file exists
        if not os.path.isfile(datafile):
            print(f"{fgColors.red}ERROR: file does not exist:\n  {datafile}")
            sys.exit()

        # NOTE: Default behavior of the script is no plot, change here to enable them
        doPlot = os.environ.get("DO_PLOT", "False").lower() in ("enabled", "en", "yes", "y", "on", "true", "t", "1")
        if doPlot:
            # this example script will generate plots
            plt.ion()

        # keep time reference at the beginning
        # NOTE: call this function to add entry in the execution time dictionary
        logExecutionTime("start")

        # -----------------------------------------------
        # --- Show user list of available functions
        # -----------------------------------------------
        sectionPrint("Show user list of available functions")
        _ = listAvailFunc(doPrint=True)

        # -----------------------------------------------
        # --- Import CSV data into a DataFrame
        # -----------------------------------------------
        sectionPrint("Import CSV data into a DataFrame")

        # load file into a DataFrame
        df = read_csv(datafile=datafile)
        logExecutionTime("df imported")

        # -----------------------------------------------
        # --- Extract Date and Time
        # -----------------------------------------------
        sectionPrint("Extract Date and Time")
        # extract the date and time at which the measures have been taken
        measureDate = getMeasureDateFromFileName(datafile)
        measureTime = getMeasureTimeFromFileName(datafile)
        print(f"Measured acquired on {measureDate} at {measureTime}")

        # -----------------------------------------------
        # --- Some measurements parameters
        # -----------------------------------------------
        sectionPrint("Some measurements parameters")
        # number of pixels enabled for the measurements
        # NOTE: The number of enabled pixels per PDC must be changed
        # according to the measurement setup.
        # Here is an example for a scintillator of 4x4 mm² on each PDC of Head 0.
        # In this example, the total number of enabled pixels is 52 x 52, but drops when disabling pixels with
        # TCR above 100 cps.
        # 3239, 2431, 1525 or 770
        numPixEn = [1525, 2412, 2382, 2376, 0, 0, 0, 0] # scintillator over 52 x 52 pixels, TCR < 100 cps
        #numPixEn = [3239, 2412, 2382, 2376, 0, 0, 0, 0] # scintillator over 52 x 52 pixels, TCR < 100 cps
        print(f"number of enabled pixels per PDC (PDC0 to PDCX): {numPixEn}")

        # default settings for the hold-off period. It can later be extracted from the measurements.
        DSUM_CLK_PRD_NS = 10
        HOLD_TIME_NS = 250
        HOLD_TIME_CLK_CYCLES = int(round(HOLD_TIME_NS/DSUM_CLK_PRD_NS))
        print(f"default hold-off time set to {HOLD_TIME_NS} ns (e.g. {HOLD_TIME_CLK_CYCLES} clock cycles)")

        # setting the expected decay time of the scintillator
        # NOTE: in this example, LYSO has been used
        SCINT_TAU_NS = 40
        SCINT_TAU_CLK_CYCLES = int(round(SCINT_TAU_NS/DSUM_CLK_PRD_NS))
        print(f"Scintillator decay time (ns): {SCINT_TAU_NS} -> {SCINT_TAU_CLK_CYCLES} clock cycles")

        logExecutionTime("initial parameters")

        # -----------------------------------------------
        # --- Extract number of events info
        # -----------------------------------------------
        sectionPrint("Extract number of events info")
        # get the number of events and a list of events
        nEvents = getNumberOfEvents(df)
        eventList = getEventList(df)
        firstEvent = getFirstEventFromEventList(eventList)
        lastEvent = getLastEventFromEventList(eventList)
        print(f"firstEvent: {firstEvent}, lastEvent: {lastEvent}")
        logExecutionTime("get number of events")

        # -----------------------------------------------
        # --- Extract PDC info
        # -----------------------------------------------
        sectionPrint("Extract PDC info")
        # Get the list of PDCs included in the DataFrame
        pdcIdxFound = getPdcIdxInData(df)
        logExecutionTime("getPdcIdxInData()")

        # -----------------------------------------------
        # --- Extract the measurement duration
        # --- and the event rate.
        # -----------------------------------------------
        sectionPrint("Extract the measurement duration\n--- and event rate")
        # From timestamps of each digital sum bin, find the total acquisition time
        acqTime = getMeasureDuration(df, clkPrd=DSUM_CLK_PRD_NS*1e-9)
        logExecutionTime("getMeasureDuration()")

        # Calculate the rate from the number of events and the duration of the acquisition
        nEventsPerSec = getMeasureEventRate(nEvents, acqTime)
        logExecutionTime("getMeasureEventRate()")

        # -----------------------------------------------
        # --- Filter data from DataFrame
        # -----------------------------------------------
        sectionPrint("Filter data from DataFrame")
        # group DataFrame data per frame index and PDC index
        dfGroup = groupDfByFrameAndPdc(df)
        logExecutionTime("group df")

        # === Applying filters on data to remove noise/unwanted data === #
        # keep only frames in which a bin is higher than threshold
        # NOTE: Here binMaxThGt can be set to match the SUM_TH greater than used for the acquistion
        keepFramesOverNoise, _ = filterFrameFromMaxBin(dfGroup, binMaxThGt=10)
        logExecutionTime("filterFrameFromMaxBin()")

        # keep only frames on which more than N bins contain data
        keepFramesNumBins, _ = filterFrameFromNumBin(dfGroup, numBinThGt=2)
        logExecutionTime("filterFrameFromNumBin()")

        # keep only frames based on the total number of detected photons in the frame.
        #   It applies a range to a preliminary estimation of the energy
        # NOTE: keeping all energy ranges for this example
        keepEnergyRng, _ = filterFrameFromTotalPhoton(dfGroup,
                                                            numTotalPhThMin=0,
                                                            numTotalPhThMax=4096)
        logExecutionTime("filterFrameFromTotalPhoton()")

        # keep only a portion of all the dataset based on a range on the index of the frame
        # NOTE: here keeping everything (begin=None, end=None keeps all the frames)
        keepFrameIdxRng, _ = filterFramesByIndex(dfGroup,
                                                       begin=None,
                                                       end=None,
                                                       eventList=eventList)
        logExecutionTime("filterFramesByIndex()")

        # keep frames based on PDC index
        # NOTE: here keeping all PDCs on Head 0
        keepPdcs, _ = filterFramesByPdcIndex(df, pdcToKeep=[0, 1, 2, 3])
        logExecutionTime("filterFramesByPdcIndex()")

        # combine all previous data filters
        # NOTE: select here the filters to combine
        filterList = [
            keepFramesOverNoise,
            keepFramesNumBins,
            keepEnergyRng,
            keepFrameIdxRng,
            keepPdcs
        ]
        keepCombined = combineFilters(filterList=filterList)
        logExecutionTime("combineFilters()")

        # apply the combined filter
        dfKeep = applyFilters(df, dfFilter=keepCombined)
        logExecutionTime("applyFilters()")

        # group DataFrame data per frame index and PDC index
        dfKeepGroup = groupDfByFrameAndPdc(dfKeep)
        logExecutionTime("group dfKeep")

        # ---------------------------------------------------------------
        # NOTE: from now, using "dfKeep" for analysis instead of "df"
        # ---------------------------------------------------------------
        sectionPrint("using 'dfKeep' for analysis instead of 'df'")

        # -----------------------------------------------
        # --- Extract the measurement duration
        # --- and the event rate.
        # -----------------------------------------------
        sectionPrint("Extract the measurement duration\n--- and event rate for 'dfKeep'")
        # After filtering noise and undesired events, recalculate the event rate
        # get the number of events
        nEventsKeep = getNumberOfEvents(dfKeep)

        # From timestamps of each digital sum bin, find the total acquisition time
        acqTimeKeep = getMeasureDuration(dfKeep, clkPrd=DSUM_CLK_PRD_NS*1e-9)
        logExecutionTime("getMeasureDuration(dfKeep)")

        # Calculate the rate from the number of events and the duration of the acquisition
        nEventsPerSecKeep = getMeasureEventRate(nEventsKeep, acqTimeKeep)
        logExecutionTime("getMeasureEventRate(dfKeep)")

        # -----------------------------------------------
        # --- Identify frames with more than one PDC
        # --- and analyse for crosstalk
        # -----------------------------------------------
        sectionPrint("Identify frames with more than one PDC")
        # Once the base filters are applied, count the number of PDCs in each frame
        analyseNumPdcPerFrame(dfKeep)
        logExecutionTime("analyseNumPdcPerFrame()")

        # Analyse crosstalk (light leak) between scintillators
        analyseCrosstalkBetweenScintillators(dfKeep, doPlot=doPlot)
        logExecutionTime("analyseCrosstalkBetweenScintillators()")

        # -----------------------------------------------
        # --- Preprocessing of digital sum data
        # -----------------------------------------------
        sectionPrint("Preprocessing of digital sum data")
        # Calcul the time difference between each bin using dataIdx, since the zeros have been
        # suppressed in the original CSV file loaded into the DataFrame.
        getDsumDt(dfKeepGroup)
        logExecutionTime("getDsumDt()")

        # For each frame, extract the index of the bin, based on dataIdx.
        getDsumBinIdx(dfKeepGroup)
        logExecutionTime("getDsumBinIdx()")

        # For each frame, extract the position of the peak and its value.
        getDsumFramePeakInfo(dfKeepGroup)
        logExecutionTime("getDsumFramePeakInfo()")

        # For each frame, extract a preliminary energy value
        getDsumPrelimEnergy(dfKeepGroup)
        logExecutionTime("getDsumPrelimEnergy()")

        # print the memory usage of the DataFrame
        pdh.print_df_usage(dfKeep)

        # -----------------------------------------------
        # --- Using digital sum data to find hold-off
        # -----------------------------------------------
        sectionPrint("Using digital sum data to find hold-off")
        # Estimate the hold-off time from the digital sum frames
        """
        HOLD_TIME_CLK_CYCLES = getHoldOffFromDsum(dfKeep, method="all", doPlot=doPlot)
        logExecutionTime("getHoldOffFromDsum()")
        """
        holdOffClkPerPdc = []
        holdOffClkCyclesColName = "holdOffClkCycles"
        for pdcIdx in getPdcIdxInData(dfKeep):
            dfPdc = dfKeep.loc[dfKeep["pdcIdx"] == pdcIdx, :]
            holdOffClkCycles = getHoldOffFromDsum(dfPdc, method="all", doPlot=doPlot)
            holdOffClkPerPdc.append({"pdcIdx": pdcIdx, holdOffClkCyclesColName: holdOffClkCycles})
            logExecutionTime(f"getHoldOffFromDsum(PDC{pdcIdx})")

        # convert list of dictionaries into a DataFrame
        dfHold = pd.DataFrame(holdOffClkPerPdc)

        # Reconstruct the cumulative sum with the hold-off duration.
        reconstDsumLevelFromEdgeAcq(dfKeep,
                                    holdOffClkCycle=dfHold,
                                    dfColName=holdOffClkCyclesColName,
                                    doPlot=False)
        logExecutionTime(f"reconstDsumLevelFromEdgeAcq()")

        # Extract a single frame from DataFrame dfKeep
        # supported values for method: "first", "last", "random", "loc", "iloc"
        dfFrame = getFramesByPdcIdxAndFrameIdx(dfKeep, method="random",
                                               numFrames=1, pdcIdx=[0])
        logExecutionTime("getFramesByPdcIdxAndFrameIdx()")

        # -----------------------------------------------
        # --- Plot a frame of digital sum and
        # --- compare the different hold-off
        # --- reconstruction methods
        # -----------------------------------------------
        if doPlot:
            # From the extracted frame (dfFrame), plot the digital sum with different
            # hold-off reconstruction methods
            HOLD_TIME_CLK_CYCLES = dfHold.loc[dfFrame["pdcIdx"].unique(), holdOffClkCyclesColName].values[0]
            plotDsumFrameHoldOff(dfFrame,
                                 holdOffClkCycle=HOLD_TIME_CLK_CYCLES,
                                 tauClkCycle=SCINT_TAU_CLK_CYCLES)
            logExecutionTime("plotDsumFrameHoldOff()")

        # -----------------------------------------------
        # --- Analyse the bin distribution
        # --- of the digital sum
        # -----------------------------------------------
        sectionPrint("Analyse the bin distribution\n--- of the digital sum")
        # extract frames with duration longer or equal to 'prd' period
        #_ = getNumFrameLongerThanPrd(dfKeepGroup, prd=HOLD_TIME_CLK_CYCLES)
        #_ = getNumFrameLongerThanPrd(dfKeepGroup, prd=dfHold[holdOffClkCyclesColName].values)
        _ = getNumFrameLongerThanPrd(dfKeepGroup, prd=dfHold, dfColName=holdOffClkCyclesColName)
        logExecutionTime("getNumFrameLongerThanPrd()")

        # count the number of frames with non consecutive bins
        # Here, used to detect frames with non consecutive data
        _ = getNumFrameDtLargerThanTh(dfKeepGroup, thGt=1)
        logExecutionTime("getNumFrameDtLargerThanTh(1)")

        # Here, finding frames with a gap in binIdx larger than hold-off time
        #dtMaxFrame, _ = getNumFrameDtLargerThanTh(dfKeepGroup, thGt=HOLD_TIME_CLK_CYCLES)
        dtMaxFrame, _ = getNumFrameDtLargerThanTh(dfKeepGroup, thGt=dfHold, dfColName=holdOffClkCyclesColName)
        logExecutionTime("getNumFrameDtLargerThanTh(dfHold)")

        # Example of combining functions
        framesWithHigherCounts, _ = filterFrameFromMaxBin(dfKeepGroup, binMaxThGt=9)
        framesNoise = combineFilters([dtMaxFrame, ~framesWithHigherCounts])
        _ = pdh.filterStats(framesNoise,
                            description="max dt > than hold-off and max dsum < 10")

        # for each frame, look at the number of different dt values.
        getNumDtPerFrame(dfKeepGroup)
        logExecutionTime("getNumDtPerFrame()")

        # frames on which the peak is not the first bin
        getFrameWithPeakNotAtFirstBin(dfKeepGroup)
        logExecutionTime("getFrameWithPeakNotAtFirstBin()")

        # Extract the maximum number of pixels triggered at the same time.
        getMaxNumberOfPixelTriggered(dfKeep, numPixEn)
        logExecutionTime("getMaxNumberOfPixelTriggered()")


        # -----------------------------------------------
        # --- Filter for hold-off
        # -----------------------------------------------
        sectionPrint("Filter data to remove values up to hold-off + margin")
        filter_peakLen = dfKeep["dtPeak"].between(-1, HOLD_TIME_CLK_CYCLES-5)
        filterList = [filter_peakLen]
        combinedFilter = combineFilters(filterList)
        dfKeep = applyFilters(dfKeep, dfFilter=keepCombined)
        logExecutionTime("applyFilters()")

        # -----------------------------------------------
        # --- Linearity correction
        # -----------------------------------------------
        sectionPrint("Linearity correction")
        # Apply linearity correction on "dsum" based on the number of pixel trigger
        _ = applyDsumLinearity(dfKeep)
        logExecutionTime("applyDsumLinearity()")

        # -----------------------------------------------
        # --- Update the energy calculation
        # -----------------------------------------------
        sectionPrint("Update the energy calculation")
        # these are examples of energy calculation methods
        #getDsumEnergy(dfKeep, col="dsum", transform="sum")
        getDsumEnergy(dfKeep, col="dsumLin", transform="sum", dtype=np.float32)
        #getDsumEnergy(dfKeep, col="dsumLevel", transform="max")
        logExecutionTime("getDsumEnergy()")

        # -----------------------------------------------
        # --- Pulse shape discrimination
        # -----------------------------------------------
        sectionPrint("Pulse shape discrimination")

        df_psd = applyPsd(dfKeep, nPrompt_l=0, nPrompt_r=1, nTotal=HOLD_TIME_CLK_CYCLES-5, colToPsd="dsumLin")
        #df_psd.to_csv(f"{datafile.split('.csv')[0]}_psd.csv")
        logExecutionTime("applyPsd()")

        if doPlot:
            # show a plot of the PSD as a function of the energy
            fig_psd, axes_psd = plotPsdFctEnergy(dfKeep)

        # -----------------------------------------------
        # --- Energy Spectrum
        # -----------------------------------------------
        sectionPrint("Energy Spectrum")
        # Examples of filters to apply on data before applying energy analysis
        # keep only frames with only one PDC
        filter_nPdc = dfKeep["nPdcsInFrame"] == 1
        # for each frame, keep only data that are less than the hold-off period in bins from the peak (with margin)
        filter_peakLen = dfKeep["dtPeak"] < HOLD_TIME_CLK_CYCLES-5

        # keep only frames for which the peak is not the first bin
        filter_peakNotFirstBin = dfKeepGroup["dtPeak"].transform("min").lt(0)

        # keep only frames for which the PSD value is within a range
        # E.g. a scintillator with a 40 ns decay on which the calculated PSD is peak bin divided by sum of all bins whould give 0.22. Adding a range ±5 -> 0.17 <= PSD <= 0.27
        filter_psd = combineFilters([dfKeep["psd"] >= 0.17, dfKeep["psd"] <= 0.27])

        # keep only frames with energy in a given range
        # NOTE: to filter background, frames with higher energy can be removed.
        filter_energy = combineFilters([dfKeep["energy"] >= 0, dfKeep["energy"] <= 2500])

        # Select here the filters to apply.
        # NOTE: these are examples, other filter types can be used
        filterList = [
            filter_nPdc,
            filter_peakLen,
            filter_peakNotFirstBin,
            filter_psd,
            filter_energy
        ]
        # NOTE: keeping all data here
        filterList = []
        combinedFilter = combineFilters(filterList)
        logExecutionTime("combineFilters(filterList)")

        # get an updated list of the PDC
        pdcList = getPdcIdxInData(dfKeep, doPrint=False)

        # generate a DataFrame of the energy spectra (one per PDC in the input DataFrame)
        dfSpectra = getEnergySpectra(dfKeep,
                                     combinedFilter=None)
        logExecutionTime("getEnergySpectra()")

        if doPlot:
            # open a figure with an ax for each PDC
            fig, axes = plt.subplots(nrows=len(pdcList), ncols=1,
                                     figsize=(12, 9), squeeze=False)
            axes = axes.flatten()

        # extract peaks from the spectra
        for idx, iPdc in enumerate(pdcList):
            dfSpectrum = extractPdcSpectrum(dfSpectra, iPdc)
            infoDict = findPhotoPeaks(dfSpectrum,
                                      spectrumBinMin=30,
                                      spectrumBinMax=-1,
                                      nBinPerBin=10,
                                      relThreshold=0.08,
                                      peaksTolerance=80)
            logExecutionTime(f"findPhotoPeaks(PDC{iPdc})")

            if doPlot:
                plotEnergySpectra(dfSpectrum,
                                  axes=[axes[idx]],
                                  peakInfos=[infoDict],
                                  plotRebin=True, plotPeaks=True, plotFit=True)


        # -----------------------------------------------
        # --- Function to plot frames from a subset
        # --- of a DataFrame
        # -----------------------------------------------
        sectionPrint("Frame analysis")
        # Extract a few frames from DataFrame dfKeep
        # supported values for method: "first", "last", "random", "loc", "iloc"
        dfSubset = getFramesByPdcIdxAndFrameIdx(dfKeep, method="random",
                                                numFrames=5, pdcIdx="all")
        logExecutionTime("getFramesByPdcIdxAndFrameIdx(all, 5)")
        dfSubsetGroup = groupDfByFrameAndPdc(dfSubset)

        if doPlot:
            # for the
            plotDsumFrames(dfSubsetGroup,
                           xCol="binIdx",
                           yColsToPlot=["dsum", "dsumLevel", "dsumLin"])

        # if plots are enabled
        if doPlot:
            # -----------------------------------------------
            # --- Function to plot a histogram
            # --- of a given column of the DataFrame
            # -----------------------------------------------
            sectionPrint("Function to plot a histogram\n--- of a given column of the DataFrame")
            # NOTE: The histogram of some columns of the DataFrame make no sense.
            cols = ["dsum", "dt", "binIdx", "peakValue",
                    "dtPeak", "dsumLevel", "energy", "nAvail",
                    "dsumLin"]
            for col in cols:
            #for col in dfKeep.columns:
                plotDfColHisto(dfKeep, col)
            logExecutionTime("plotDfColHisto()")

        # If figures are open, prevent the end of the script
        if plt.get_fignums():
            input(f"Analysis Completed. Enter to continue")

    except KeyboardInterrupt:
        print(f"{fgColors.yellow}\nInterrupted.{fgColors.yellow}")
        sys.exit()

    finally:
        # table with all the execution time entries
        printExecutionTime()





