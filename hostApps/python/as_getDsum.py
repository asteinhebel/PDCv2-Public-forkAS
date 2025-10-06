#----------------------------------------------------------------------------------
#-- Company: GRAMS
#-- Designer: Tommy Rossignol
#--
#-- Create Date: 2024-05-21
#-- Description:
#--     Importing PDC digital sum (Dsum) from HDF5 file (from Controller)
#--     This script do not configure the PDCs nor the Controller.
#--     It only display the results from the digital sum.
#--     For data to be read user must do the following:
#--         1 - NFS must be configured on the Zynq board.
#--         2 - 'dataReader' app must run on the Zynq board.
#--         3 - an hex to HDF5 app must run on the NFS server (e.g. dma2h5, hexRead)
#--         4 - a script must be run to configure the Controller and PDCs for acquisition
#--
#-- Dependencies:
#-- Revision:
#-- Revision 1.0 - File Created
#-- Revision 2.0 - Updated for sharing on pdcv2-public
#-- Additional Comments:
#--     NOTE: Since this app is in Python, the execution speed is limited.
#--           Once an HDF5 file is read, it is deleted (deleteAfter parameter of h5Reader).
#--           Be carefull not to generate too much data too quickly.
#--           If the data generation speed is higher than the delete speed,
#--           you may fill your server disk, or get errors.
#--     NOTE: Possible options to specify the input HDF5 directory (HDF5_DATA_DIR):
#--           1- from HDF5_DATA_DIR environment variable
#--           2- from default directory (defaultHdf5Dir)
#--           3- from the hardcoded path (not recommended)
#--           HDF5_DATA_DIR variable is setup in module zynqHelper.py
#--     NOTE: Depending on the script generating your data, you may want to change
#--           the following user settings:
#--                SYS_CLK_PRD
#--                DSUM_SAMPLE_NCLK
#----------------------------------------------------------------------------------
import os, sys
from pathlib import Path
import time, datetime
import pandas as pd
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import matplotlib.patches as patches
import yaml

from enum import IntEnum

# to add buttons to menu bar
from matplotlib.backend_tools import ToolBase, ToolToggleBase
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.rcParams['toolbar'] = 'toolmanager'

# to save data
from itertools import zip_longest

import matplotlib # for get_backend()
import matplotlib.backends as backends

# custom modules
from modules.fgColors import fgColors
from modules.zynqEnvHelper import PROJECT_PATH, HOST_APPS_PATH, USER_DATA_DIR, HDF5_DATA_DIR
from modules.h5Reader import *

try:
    scriptName = os.path.basename(__file__)
    scriptAbsPath = os.path.dirname(os.path.abspath(__file__))
except NameError:
    scriptName = "fileNameNotFound.py"
    scriptAbsPath = os.path.abspath("./")

# -----------------------------------------------
# --- parse configuration values from YAML
# -----------------------------------------------
with open('config_all.yaml','r') as f:
    config_in = yaml.safe_load(f)
    pdcConfig_in = config_in['pdcConfig']
    instConfig_in = config_in['instrumentConfig']
    dataConfig_in = config_in['dataConfig']
    config_out = config_in.copy()


# -----------------------------------------------
# --- Global settings
# -----------------------------------------------
parsedH5 = []

def dbgPrint(message, enabled=False):
    if enabled: print(message)

# default path to open the menu to save data as a CSV
saveCsvInitialDir = os.path.join(USER_DATA_DIR, os.path.splitext(scriptName)[0])
Path(saveCsvInitialDir).mkdir(parents=True, exist_ok=True)

# -----------------------------------------------
# --- User settings - SYS_CLK_PRD
# -----------------------------------------------
# NOTE: Since the acquisition is triggered from another
#       script, this apps is not aware of the system clock period.
#       Using default value of 10 ns (100 MHz)
# system clock frequency of the Controller/PDCs
if os.environ.get("SYS_CLK_PRD") is not None:
    # using environment variable
    SYS_CLK_PRD = os.environ['SYS_CLK_PRD']
else:
    # default setting
    SYS_CLK_PRD=10.0e-9


# -----------------------------------------------
# --- User settings - DSUM_SAMPLE_NCLK
# -----------------------------------------------
# NOTE: This app is not designed to automatically get the sampling period
#       of the digital sum. It is hardcoded here for now.
#       DSUM_SAMPLE_NCLK =  1 -> one digital clock sample each SYS_CLK_PRD
#       DSUM_SAMPLE_NCLK = 10 -> one digital clock sample each 10 SYS_CLK_PRD
# number of clock cycles for each digital sum sample
if os.environ.get("DSUM_SAMPLE_NCLK") is not None:
    DSUM_SAMPLE_NCLK = os.environ['DSUM_SAMPLE_NCLK']
else:
    # default setting
    DSUM_SAMPLE_NCLK=1


# -----------------------------------------------
# --- User settings - HDF5_DATA_DIR
# -----------------------------------------------
# see line:
# from modules.zynqHelper import HDF5_DATA_DIR # to specify the HDF5 input path

# number of supported PDCs
N_PDC_MAX = 4 #8

# database of the data
dp = None
fname = "_"+dataConfig_in['name'] if len(dataConfig_in['name'])>0 else ""

# -----------------------------------------------
# --- Class for the data
# -----------------------------------------------
class dsumPlotter:
    def __init__(self, figName, nPdcMax):
        """
        create empty object with no data, but figure properly formatted
        """
        self.run = True
        self.figName = figName
        self.nPdcMax = nPdcMax
        self.label = "DSUM"
        self.nEvents = [0]*self.nPdcMax
        self.DSUM_time = [0]*self.nPdcMax
        self.DSUM_data = [0]*self.nPdcMax
        self.DSUM_PLOT_time = [0]*self.nPdcMax
        self.DSUM_PLOT_data = [0]*self.nPdcMax

    def getAllPdcData(self, db):
        """
        read from a H5 data base object (see settings class)
        """

        # init empty data
        self.DSUM_time = [0]*self.nPdcMax
        self.DSUM_data = [0]*self.nPdcMax

        # open hdf5 file
        db.h5Open()

        # get content
        for iPdc in range(self.nPdcMax):
            [self.DSUM_time[iPdc], self.DSUM_data[iPdc]] = db.getPdcDsum(iPdc=iPdc)

        # close hdf5 file
        db.h5Close()

        # process the data only once the HDF5 file is closed
        for iPdc in range(self.nPdcMax):
            if (not np.shape(self.DSUM_data[iPdc]) == ()):
                self.nEvents[iPdc] += 1

        #equivalent of update method
        try:
            self.DSUM_PLOT_data += np.array(self.DSUM_data)
        except ValueError: #first entry
            self.DSUM_PLOT_data = np.array(self.DSUM_data)
        self.DSUM_PLOT_time = self.DSUM_time


    def saveCsv(self, extraname:str=""):
        """
        save data to a csv file
        """
        # init empty DataFrame
        df = pd.DataFrame()

        # for each PDC, check if plot data is available
        for iPdc in range(dp.nPdcMax):
            if dp.nEvents[iPdc] > 0:
                df.insert(loc=len(df.columns), column=f"time{iPdc}", value=self.DSUM_PLOT_time[iPdc])
                df.insert(loc=len(df.columns), column=f"data{iPdc}", value=self.DSUM_PLOT_data[iPdc])


        dateStr=datetime.datetime.now().strftime("%Y%m%d_%Hh%Mm%S")
        ename = "_"+extraname if len(extraname)>0 else ""
        filename = f"{dateStr}_DSUM{ename}{fname}.csv"
        datafile = os.path.join(Path(saveCsvInitialDir), filename)
        if not os.path.isdir(saveCsvInitialDir):
            print(f"{fgColors.red}ERROR: Specified path '{datafile}' is not valid. No data exported.{fgColors.endc}")
            return ""

        # if there are data to export
        if df.size > 0:
            print(f"{fgColors.green}Saving data to file {datafile}{fgColors.endc}")
            df.to_csv(datafile, sep=',', index=False, float_format="%.3E")
        else:
            print(f"{fgColors.bYellow}WARNING: No data to save.{fgColors.endc}")


# ---------------------------------------
# --- Script main execution
# ---------------------------------------
try:
    #   loop until user Ctrl+C
    dp = dsumPlotter(figName="DSUM PLOTTER",
                     nPdcMax=N_PDC_MAX)
    print(f"{fgColors.bBlue}Collecting DSUM data.{fgColors.endc}")

    countedEvents=0
    maxEvents = dataConfig_in['maxEvents'] if dataConfig_in['maxEvents']>0 else 1e9
    t_start = datetime.datetime.now()
    t_elapse = datetime.timedelta(seconds=dataConfig_in['maxRunTime']) if dataConfig_in['maxRunTime']>0 else datetime.timedelta(seconds=1e9)
    
    #remove existing h5 files
    for f in os.listdir(HDF5_DATA_DIR):
        os.remove(os.path.join(HDF5_DATA_DIR, f))

    while (countedEvents<maxEvents) and (datetime.datetime.now()-t_start<t_elapse):

        #begin data collection 
        db = h5Reader(deleteAfter=True,
                      #hfRelPath="HDF5",
                      hfAbsPath=HDF5_DATA_DIR,
                      sysClkPrd=SYS_CLK_PRD,
                      dsumPrd=DSUM_SAMPLE_NCLK,
                      hfFile="")

        if db.newFileReady():
            # -----------------------------------------------
            # --- Open HDF5 file to get Controller Data
            # -----------------------------------------------
            # get all PDC data for a given event in db
            dp.getAllPdcData(db=db)
            countedEvents+=1
        else: # nothing new to plot
            time.sleep(0.1)

except (KeyboardInterrupt, SystemExit) as ex:
    if "dp" in locals():
        dp.run = False
    if isinstance(ex, SystemExit):
        print(f"\n{fgColors.yellow}Program interrupted: exit program{fgColors.endc}")
    else:
        print(f"\n{fgColors.yellow}Keyboard Interrupt: exit program{fgColors.endc}")
finally:
    #automatically save data
    t_end = datetime.datetime.now()
    str_duration = f"{(t_end-t_start).seconds}s_{countedEvents}evts"
    dp.saveCsv(extraname=str_duration)
    sys.exit()

