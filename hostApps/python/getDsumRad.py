#----------------------------------------------------------------------------------
#-- Company: GRAMS
#-- Designer: Tommy Rossignol
#--
#-- Create Date: 2025-09-24
#-- Description:
#--      Using python ssh libraries to send remote commands to the ZCU102
#--      This script prepare the Controller and the PDCs for an acquisition.
#--      Based on the number of photons detected, an acquisition is started.
#--      Results are stored to a CSV file, to be analysed with another script.
#--
#-- Dependencies:
#-- Revision:
#-- Revision 1.0 - File Created
#-- Revision 2.0 - Updated for sharing on pdcv2-public
#-- Additional Comments:
#--
#----------------------------------------------------------------------------------
import sys, os
import numpy as np
import random
import time
import datetime
from itertools import chain
import statistics
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# custom modules
from pdcv2_modules.fgColors import fgColors
from pdcv2_modules.zynqEnvHelper import PROJECT_PATH, HOST_APPS_PATH, USER_DATA_DIR, HDF5_DATA_DIR
import pdcv2_modules.sshClientHelper as sshClientHelper
import pdcv2_modules.systemHelper as systemHelper
import pdcv2_modules.pixMap as pixMap
from pdcv2_modules.zynqCtlPdcRoutines import initCtlPdcFromClient, packetBank
from pdcv2_modules.zynqDataTransfer import zynqDataTransfer
from pdcv2_modules.systemHelper import sectionPrint
from pdcv2_modules.pdcHelper import *
#from modules.zynqHelper import *
from pdcv2_modules.h5Reader import *

import pdcv2_modules.pdcSpadFunctions as pdcSpadFunctions

try:
    scriptName = os.path.basename(__file__)
except NameError:
    scriptName = "fileNameNotFound.py"

# -----------------------------------------------
# --- Global vars
# -----------------------------------------------
# NOTE: user can set hold, recharge and flag time using the following environment variables:
#       HOLD_TIME_NS, RECH_TIME_NS, FLAG_TIME_NS


# NOTE: set environment variable SPAD_BIAS_V to store it in data file name
spadBiasStr=""
if os.environ.get("SPAD_BIAS_V") is not None:
    spadBias = os.environ['SPAD_BIAS_V']
    if '.' in spadBias:
        spadBiasStr = "_" + spadBias.replace('.', "V")
    elif ',' in spadBias:
        spadBiasStr = "_" + spadBias.replace(',', "V")
    else:
        spadBiasStr = "_" + spadBias + "V"
    print(f"SPAD bias voltage set to {spadBias} V")


#NOTE: set environment variable HEAD_ID to store it in data file name
#      only set the integer value (e.g. 48)
headStr = ""
if os.environ.get("HEAD_ID") is not None:
    headId = os.environ['HEAD_ID']
    headStr = f"H{headId}_"
    print(f"Using head {headId}")


# NOTE: set environment variable RAD_SOURCE to specify
#       the radiation source used for the measurement
#radSource = "Co57"
#radSource = "Cs137"
#radSource = "Ge68"
#radSource = "Am241"
#radSource = "xray"
radSource = "background"
if os.environ.get("RAD_SOURCE") is not None:
    radSource = os.environ['RAD_SOURCE']
print(f"Radiation source: {radSource}")


# NOTE: when set to analogOnly, no csv file is generated.
#       Flag output is set as flag, instead of dsum threshold
analogOnly = os.environ.get("ANALOG_ONLY", default="False").lower().strip()
if analogOnly in ('true', '1'):
    analogOnly = True
else:
    analogOnly = False
print(f"analogOnly: {analogOnly}")


# NOTE: specify a TCR file (from getSpadTcrUSingFlag.py) to set which pixels to enable
tcrFile = None # default value will throw an error
if os.environ.get("TCR_FILE") is not None:
    tcrFile = os.environ["TCR_FILE"]
    # check if full path is specified
    if os.path.isfile(tcrFile):
        # user specified with full path
        print(f"TCR_FILE specified as absolute path:\n  {tcrFile}")
    else:
        # check if user only specified file name, look into default path
        tcrFile = os.path.join(USER_DATA_DIR, "TCR", "getSpadTcrUsingFlag", tcrFile)
        if os.path.isfile(tcrFile):
            # user specified file name from default path
            print(f"TCR_FILE specified as file name in default path:\n  {tcrFile}")
        else:
            print(f"{fgColors.red}ERROR: could not find specified file {tcrFile}{fgColors.endc}")
            sys.exit()

if tcrFile is None:
    print(f"{fgColors.red}ERROR: no TCR file specified (use TCR_FILE=filename in system call){fgColors.endc}")
    sys.exit()

# -----------------------------------------------
# --- open a connection with the ZCU102 board
# -----------------------------------------------
sectionPrint("open a connection with the ZCU102 board")
# parameters of the ZCU102 board
# open a client based on its name in the ssh config file
client = sshClientHelper.sshClientFromCfg(hostCfgName="zcudev")

# -----------------------------------------------
# --- prepare Zynq platform
# -----------------------------------------------
sectionPrint("prepare Zynq platform")
zynq = zynqDataTransfer(sshClientZynq=client)
#zynq.hexAppName = "hexReadMax"
#zynq.init() # a custom setting is required here, not using init()
zynq.initNfs()

# make sure no file remains from previous test run
zynq.cleanDataPath()

# start dataReader app on ZCU102
zynq.initDataReader(dataReaderLaunch=True)

# dsum module settings
#AS - must create this dir if doesn't exist
CSV_DATA_DIR = os.path.join(USER_DATA_DIR, f"DSUM_CSV_3D_{radSource}")
DATE_STR = datetime.datetime.now().strftime("%Y%m%d_%Hh%Mm%S")
# DATA_TYPE: "all", "NZ", "NZKF", "dt", "max"
DATA_TYPE = "NZKF"
BIN_IDX_MODE = "time" # "continuous", "frame", "time"
DATA_FILE_NAME = f"{DATE_STR}_{os.path.splitext(scriptName)[0]}_{headStr}{DATA_TYPE}_{BIN_IDX_MODE}{spadBiasStr}.csv"
dsumCsvFile = os.path.join(CSV_DATA_DIR, DATA_FILE_NAME)

# zpp module settings
CSV_ZPP_DIR = os.path.join(USER_DATA_DIR, f"ZPP_CSV_3D_{radSource}")
#zppOptions = f"-M zpp " \
#             f"-o {CSV_ZPP_DIR} " \
#             f"-f getDsumOnAnyFlagZpp.csv"
zppOptions = "" # disabled zpp readout

# hexRead options (not saving to HDF5 file)
if not analogOnly:
    zynq.initHex(autoStart=True,
                archive=False,
                printParsed=False,
                exportH5=False,
                hexReadExtraArgs=f"-M dsum " \
                                f"-o {CSV_DATA_DIR} " \
                                f"-f {DATA_FILE_NAME} " \
                                f"-b {BIN_IDX_MODE} -t {DATA_TYPE} " \
                                f"{zppOptions} ")



# -----------------------------------------------
# --- prepare controller for acquisition
# -----------------------------------------------
# NOTE: select the PDC to use:
#       pdcEn=0x1 -> PDC0
#       pdcEn=0x2 -> PDC1
#       pdcEn=0x4 -> PDC2
#       pdcEn=0x8 -> PDC3
#       pdcEn=0xF -> PDC0, PDC1, PDC2, PDC3
# NOTE: set environment variable PDC_EN tp set which PDCs to use
pdcEn = int(os.environ.get("PDC_EN", default="0xF"), 0)
icp = initCtlPdcFromClient(client=client, sysClkPrd=10e-9, pdcEn=pdcEn)

# -----------------------------------------------
# --- set system clock period
# -----------------------------------------------
icp.setSysClkPrd()

# -----------------------------------------------
# --- reset of the controller
# -----------------------------------------------
icp.resetCtl()

# -----------------------------------------------
# --- configure controller packet
# -----------------------------------------------
# NOTE always set SCSA register first to store other configuration registers in HDF5
# configure CFG_STATUS_A
    # 0x8000 = PDC_CFG
    # 0x4000 = CTL_CFG
    # 0x2000 = PDC_STATUS
    # 0x1000 = PDC_STATUS_ALL
    # 0x0007 = ALL CTL_STATUS
SCSA = 0x1000
# configure CTL_DATA_A
    # 0x0001 = GBL_CTL_TDC
SCDA = 0x0001
# configure PDC_DATA_A
    # 0x0100 = DSUM
    # 0x00F7 = ZPP
SPDA = 0x0100
icp.setCtlPacket(bank=packetBank.BANKA, SCS=SCSA, SCD=SCDA, SPD=SPDA)

# -----------------------------------------------
# --- configure the Controller FSM
# -----------------------------------------------
sectionPrint("configure the Controller FSM")
# NOTE: Always start by configuring FSM registers with the higher register addresses first.
#       This way, the FSM will not start while being configured.
# NOTE: All registers with a CANCEL_MODE selection (FACQ, FSM_ACQ_CANCEL_MODE), use the following bits:
#       bit 0 - NO_COINCIDENCE          NOCN
#       bit 1 - SUM THRESHOLD AND       STHA
#       bit 2 - SUM THRESHOLD ANY (OR)  STHO
#       bit 3 - EXTERNAL 1              EXT1
#       bit 4 - EXTERNAL 2              EXT2
#       bit 5 - PDC_FUNC_ALL            PDCA
#       bit 6 - PDC_FUNC_ANY (OR)       PDCO


# Number of samples to acquire
# NOTE: nAcqSamples is the sum of nAcqSamplesFast and nAcqSamplesSlow.
#       This value should not be larger than 128.
#       nAcqSamplesFast = 128, means 128 samples each Controller clock cycle (e.g. default is 10 ns)
#       nAcqSamplesSlow = 128, means 128 samples, but with a spacing defined by Controller register SLW0 and SLW1
#       User can decide to split between fast and slow (e.g. 28 fast and 100 slow).
#       nAcqSamplesHistory is the number of samples before the trigger of the acquisition.
#       To get nAcqSamplesHistory, acquisition must run continuously to get samples before the trigger.
nAcqSamplesFast = 0
nAcqSamplesSlow = 128
#nAcqSamplesSlow = 64
#nAcqSamplesSlow = 20
nAcqSamples = nAcqSamplesFast + nAcqSamplesSlow
nAcqSamplesHistory = 28
#nAcqSamplesHistory = 14
#nAcqSamplesHistory = 0
print(f"nAcqSamplesFast     = {nAcqSamplesFast}")
print(f"nAcqSamplesSlow     = {nAcqSamplesSlow}")
print(f"nAcqSamples         = {nAcqSamples}")
print(f"nAcqSamplesHistory  = {nAcqSamplesHistory}")
print("")
# depending if ACQ is fast or slow or both, this is the trigger to start the acquistion
#FSM_ACQ_TRG_MODE = 2    # NOTE CHANGE HERE
FSM_ACQ_TRG_MODE = 5    # NOTE CHANGE HERE -----------
#FSM_ACQ_TRG_MODE = 6    # NOTE CHANGE HERE
#FSM_ACQ_TRG_MODE = 1    # NOTE CHANGE HERE
                        # 0 = FLAG,
                        # 1 = COINCIDENCE_OK
                        # 2 = EXTERNAL 1
                        # 3 = EXTERNAL 2
                        # 4 = SUM_THRESHOLD_AND
                        # 5 = SUM_THRESHOLD_OR
                        # 6 = COUNTER
                        # 7 = BYPASS


# configure FSM timeout register (TOUT)
FSM_TIMEOUT_AUTO = 0    # 1 = timeout period is automatically calculated
FSM_TIMEOUT_A1O0 = 1    # 1 = All PDC must timeout to sent a timeout, 0 = only 1 PDC to sent a timeout
FSM_TIMEOUT_RETRY = 1   # number of timeout retries before sending a timeout
FSM_TIMEOUT_PRD = 48    # period in clock cycles of the timeout
toutReg = ((FSM_TIMEOUT_AUTO&0x1)<<15) | ((FSM_TIMEOUT_A1O0&0x1)<<14) | \
          ((FSM_TIMEOUT_RETRY&0x3)<<12) | (FSM_TIMEOUT_PRD&0x3FF)
client.runPrint(f"ctlCfg -a TOUT -r 0x{toutReg:04x} -g")


# configure FSM END register (FEND)
FSM_END_EN = 1      # 1 = enable end state of the FSM
FSM_END_MODE = 2    # 0 = END_DONE_CNT, 1 = END_DONE_DAQ, 2 = END_DONE_PCK
FSM_END_DELAY = 0   # delay to use to leave the FSM end state. (used with END_DONE_CNT)
fendReg = ((FSM_END_EN&0x1)<<15) | ((FSM_END_MODE&0x3)<<8) | (FSM_END_DELAY&0xFF)
client.runPrint(f"ctlCfg -a FEND -r 0x{fendReg:04x} -g")


# configure FSM transmit register 1 (FTX1)
FSM_TX_CNL_SEND_RSTN = 0    # 1 = when cancelling the transmit, sent a RSTN command to the PDCs
FSM_TX_CANCEL_MODE = 0      # bitwise setting to set the conditions on which to sent a RSTN command to the PDCs
FSM_TX_N = nAcqSamples      # number of samples of the digital sum to transmit from the PDC to the Controller
ftx1Reg = ((FSM_TX_CNL_SEND_RSTN&0x1)<<15) | ((FSM_TX_CANCEL_MODE&0x7F)<<8) | (FSM_TX_N&0xFF)
client.runPrint(f"ctlCfg -a FTX1 -r 0x{ftx1Reg:04x} -g")

# configure FSM transmit register 0 (FTX0)
FSM_TX_EN = 1           # 1 = enable transmit state of the Controller FSM
FSM_TX_CMD_MODE = 1     # 0 = done with a counter after the command, 1 = done when TX command is sent
FSM_TX_MODE = 2         # 0 = DONE_PDC_FUNC_ALL,
                        # 1 = DONE_PDC_FUNC_ANY,
                        # 2 = DONE_CNT, (completion from counter FSM_TX_N)
                        # 3 = DONE_AUTO (completion auto calculated from given settings)
FSM_TX_CMD_DELAY = 0    # delay in clock cycle to wait after sending the command
ftx0Reg = ((FSM_TX_EN&0x1)<<15) | ((FSM_TX_CMD_MODE&0x1)<<10) | ((FSM_TX_MODE&0x3)<<8) | (FSM_TX_CMD_DELAY&0xFF)
client.runPrint(f"ctlCfg -a FTX0 -r 0x{ftx0Reg:04x} -g")


# configure FSM acquisition-transmission register 1 (ATX1)
FSM_ATX_CNL_SEND_RSTN = 0    # 1 = when cancelling the acquisition-transmission, sent a RSTN command to the PDCs
FSM_ATX_CANCEL_MODE = 0      # bitwise setting to set the conditions on which to sent a RSTN command to the PDCs
FSM_ATX_N = 0                # number of samples of the digital sum to acquire and transmit from the PDC
atx1Reg = ((FSM_ATX_CNL_SEND_RSTN&0x1)<<15) | ((FSM_ATX_CANCEL_MODE&0x7F)<<8) | (FSM_ATX_N&0xFF)
client.runPrint(f"ctlCfg -a ATX1 -r 0x{atx1Reg:04x} -g")

# configure FSM acquisition-transmission register 0 (ATX0)
FSM_ATX_EN = 0          # 1 = enable acquisition-transmission state of the Controller FSM
FSM_ATX_MODE = 0        # 0 = DONE_PDC_FUNC_ALL,
                        # 1 = DONE_PDC_FUNC_ANY,
                        # 2 = DONE_CNT, (completion from counter FSM_TX_N)
                        # 3 = DONE_AUTO (completion auto calculated from given settings)
FSM_ATX_CMD_MODE = 0    # 0 = FLAG,
                        # 1 = COINCIDENCE_OK
                        # 2 = EXTERNAL 1
                        # 3 = EXTERNAL 2
                        # 4 = SUM_THRESHOLD_AND
                        # 5 = SUM_THRESHOLD_OR
                        # 6 = COUNTER
                        # 7 = BYPASS
FSM_ATX_CMD_DELAY = 0   # delay in clock cycle to wait after sending the command
atx0Reg = ((FSM_ATX_EN&0x1)<<15) | ((FSM_ATX_MODE&0x3)<<11) | \
          ((FSM_ATX_CMD_MODE&0x7)<<8) |  (FSM_ATX_CMD_DELAY&0xFF)
client.runPrint(f"ctlCfg -a ATX0 -r 0x{atx0Reg:04x} -g")


# configure FSM slow acquisition register 1 (SLW1)
FSM_ACQ_SLOW_PRD = 0; # NOTE: change PRD setting here, any setting between 0 and 255 is valid
FSM_ACQ_SLOW_N = nAcqSamplesSlow
slw1Reg = ((FSM_ACQ_SLOW_PRD&0xFF)<<8) | (FSM_ACQ_SLOW_N&0xFF)
client.runPrint(f"ctlCfg -a SLW1 -r 0x{slw1Reg:04x} -g")

# configure FSM slow acquisition register 0 (SLW0)
FSM_ACQ_SLOW_EN = 1 if nAcqSamplesSlow > 0 else 0   # enable slow acquisition state the the Controller FSM
FSM_ACQ_SLOW_MODE = 0   # 0 = FLAG,
                        # 1 = COINCIDENCE_OK
                        # 2 = EXTERNAL 1
                        # 3 = EXTERNAL 2
                        # 4 = SUM_THRESHOLD_AND
                        # 5 = SUM_THRESHOLD_OR
                        # 6 = COUNTER
                        # 7 = BYPASS
FSM_ACQ_SLOW_DELAY = 0  # delay in clock cycle to wait at the end of the state
if nAcqSamplesFast == 0:
    FSM_ACQ_SLOW_MODE = FSM_ACQ_TRG_MODE
    FSM_ACQ_SLOW_DELAY = nAcqSamplesSlow - nAcqSamplesHistory
slw0Reg = ((FSM_ACQ_SLOW_EN&0x1)<<15) | ((FSM_ACQ_SLOW_MODE&0x7)<<8) | (FSM_ACQ_SLOW_DELAY&0xFF)
client.runPrint(f"ctlCfg -a SLW0 -r 0x{slw0Reg:04x} -g")


# configure FSM fast acquisition register 1 (FST1)
FSM_ACQ_FAST_N = nAcqSamplesFast
fst1Reg = (FSM_ACQ_FAST_N&0xFF)
client.runPrint(f"ctlCfg -a FST1 -r 0x{fst1Reg:04x} -g")

# configure FSM fast acquisition register 0 (FST0)
FSM_ACQ_FAST_EN = 1 if nAcqSamplesFast > 0 else 0   # enable fast acquisition state the the Controller FSM
FSM_ACQ_FAST_MODE = FSM_ACQ_TRG_MODE
                        # 0 = FLAG,
                        # 1 = COINCIDENCE_OK
                        # 2 = EXTERNAL 1
                        # 3 = EXTERNAL 2
                        # 4 = SUM_THRESHOLD_AND
                        # 5 = SUM_THRESHOLD_OR
                        # 6 = COUNTER
                        # 7 = BYPASS
FSM_ACQ_FAST_DELAY = nAcqSamplesFast - nAcqSamplesHistory # delay in clock cycle to wait at the end of the state
fst0Reg = ((FSM_ACQ_FAST_EN&0x1)<<15) | ((FSM_ACQ_FAST_MODE&0x7)<<8) | (FSM_ACQ_FAST_DELAY&0xFF)
client.runPrint(f"ctlCfg -a FST0 -r 0x{fst0Reg:04x} -g")


# configure FSM acquisition register (FACQ)
FSM_ACQ_CNL_SEND_RSTN = 0       # 1 = when cancelling the acquisition, sent a RSTN command to the PDCs
FSM_ACQ_CANCEL_MODE = 0         # bitwise setting to set the conditions on which to sent a RSTN command to the PDCs
FSM_ACQ_OVERLAP = 0             # number of samples to overlap while sending the next command
FSM_ACQ_CMD_MODE = 7            # 0 = FLAG # NOTE change here
#FSM_ACQ_CMD_MODE = 2            # 0 = FLAG # NOTE change here
                                # 1 = COINCIDENCE_OK
                                # 2 = EXTERNAL_1
                                # 3 = EXTERNAL_2
                                # 7 = BYPASS
facqReg = ((FSM_ACQ_CNL_SEND_RSTN&0x1)<<15) | ((FSM_ACQ_CANCEL_MODE&0x7F)<<8) | \
          ((FSM_ACQ_OVERLAP&0xF)<<4) | (FSM_ACQ_CMD_MODE&0x7)
client.runPrint(f"ctlCfg -a FACQ -r 0x{facqReg:04x} -g")


# configure FSM misc register (FSMM)
# NOTE: keep FSM disabled for now, will be enabled when everything else is ready
FSM_TX_ALL_PDC = 1      # force to transmit from all PDCs, not only those with data
FSM_ACQ_ALL_PDC = 1     # force to acquire from all PDCs, not only those with events
FSM_SEQ_START_MODE = 0  # 0 = DISABLED
                        # 1 = COMMAND only
                        # 2 = TIMER
                        # 3 = BYPASS (as soon as possible)
fsmmReg = ((FSM_TX_ALL_PDC&0x1)<<8) | ((FSM_ACQ_ALL_PDC&0x1)<<4) | (FSM_SEQ_START_MODE&0x3)
client.runPrint(f"ctlCfg -a FSMM -r 0x{fsmmReg:04x} -g")


# -----------------------------------------------
# --- set delay of CFG_DATA pins
# -----------------------------------------------
icp.setDelay(signal="CFG_DATA", delay=300)

# -----------------------------------------------
# --- check for power good
# -----------------------------------------------
icp.checkPowerGood()

# -----------------------------------------------
# --- enable CFG_RTN_EN
# -----------------------------------------------
icp.setCfgRtnEn()

# -----------------------------------------------
# --- prepare PDC for configuration
# -----------------------------------------------
icp.preparePDC()

# -----------------------------------------------
# --- Using all of the 4096 pixels,
# --- set the number of SPADS to 4096
# -----------------------------------------------
# NOTE: Default icp.nSpad is 64 (2D CMOS SPAD)
#       To user all of the 4096 pixels,
#       uncomment the following line
icp.nSpad = 4096

# --------------------------
# --- configure the PDCs ---
# --------------------------
sectionPrint("configure the PDCs")
PDC_SETTING = pdc_setting()
client.runPrint("ctlCmd -c MODE_CFG")  # set PDCs to configuration mode

# === PIXL REGISTER ===
print("\n=== PIXL REGISTER ===")
# active quenching of the front-end
ACTIVE_QC_EN = 1; # 0=disabled/passive, 1=enabled/active
# trigger using QC front-end (FE) or digital only (DGTL)
TRG_DGTL_FEN = 0; # 0=FE, 1=DGTL
# enable flag output of the pixel
FLAG_EN = 1; # 1=enabled, 0=disabled
# EDGE_LVLN and DIS_MEM on synchronizer
EDGE_LVLN = 1
DIS_MEM = 1
PIXL = ((DIS_MEM<<13) + (EDGE_LVLN<<12) + (FLAG_EN<<8) + (TRG_DGTL_FEN<<4) + (ACTIVE_QC_EN<<1))
client.runPrint(f"pdcCfg -a PIXL -r 0x{PIXL:04x} -g")  # configure pixel register
PDC_SETTING.PIXL = PIXL

# === TIME REGISTER ===
print("\n=== TIME REGISTER ===")
HOLD_TIME = float(os.environ.get("HOLD_TIME_NS", default=250.0))
RECH_TIME = float(os.environ.get("RECH_TIME_NS", default=10.0))
FLAG_TIME = float(os.environ.get("FLAG_TIME_NS", default=2.0))
client.runPrint(f"pdcTime --hold {HOLD_TIME} --rech {RECH_TIME} --flag {FLAG_TIME} -g")
PDC_SETTING.TIME = client.runReturnSplitInt('pdcTime -g')

# === ANLG REGISTER ===
print("\n=== ANLG REGISTER ===")
#ANLG = 0x0000; # disabled
ANLG = 0x0001; # full amplitude (~30 µA)
client.runPrint(f"pdcCfg -a ANLG -r 0x{ANLG:04x} -g")  # set analog monitor
PDC_SETTING.ANLG = ANLG

# === STHH REGISTER ===
#print("\n=== STHL REGISTER ===")
#STHH = 0x1FFF
#client.runPrint(f"pdcCfg -a STHH -r 0x{STHH:04x} -g")  # set sum threshold high

# === STHL REGISTER ===
print("\n=== STHL REGISTER ===")
EN_SUM_TH = 1   # enables both sum threshold registers STHH and STHL
# NOTE: change SUM_THL to optimize the trigger rate of events. A higher value filters the DCR, but might remove valid events
SUM_THL = 6 # number of photon to exceed in digital sum to raise SUM_GT
STHL = ((EN_SUM_TH&0x1)<<15) | (SUM_THL&0x1FFF)
client.runPrint(f"pdcCfg -a STHL -r 0x{STHL:04x} -g")  # set sum threshold low
PDC_SETTING.STHL = STHL
# NOTE: 100 TCR max per pixel on 2500 pixels = 250 000 cps = 1 count each 4 us.
#       With a hold off of 250 us, STHL must be > 16

# === XXXX REGISTER ===
# skipping registers ACQA to DBGC

# === FIFO REGISTER ===
print("\n=== FIFO REGISTER ===")
OVERWRITE_DIS = 0   # disable overwriting of the FIFO during acquisition
FIFO_EN = 1         # to enable the FIFO
SEL_FIFO_IN = 0     # 0 = digital sum
                    # 1 = debug counter
                    # 2 = debug register
                    # 3 = unused (zeros)
#FIFO_DEPTH = 127    # the size of the FIFO - 1. (e.g. FIFO_DEPTH=127 -> FIFO has 128 bins)
if nAcqSamples > 0:
    FIFO_DEPTH = nAcqSamples-1    # the size of the FIFO - 1. (e.g. FIFO_DEPTH=127 -> FIFO has 128 bins)
else:
    FIFO_DEPTH = 127

FIFO = ((OVERWRITE_DIS&0x1)<<13) | ((FIFO_EN & 0x1)<<12) | ((SEL_FIFO_IN&0x3)<<8) | (FIFO_DEPTH&0x7F)
client.runPrint(f"pdcCfg -a FIFO -r 0x{FIFO:04x} -g")
PDC_SETTING.FIFO = FIFO

# === OUTD REGISTER ===
print("\n=== OUTD REGISTER ===")
#DATA_FUNC = OUT_MUX.FLAG
#DATA_FUNC = OUT_MUX.TRG
DATA_FUNC = OUT_MUX.DATA
#DATA_FUNC = OUT_MUX.SUM_GT
#DATA_FUNC = OUT_MUX.VSS
#DATA_FUNC = OUT_MUX.VDD
OUTD = (DATA_FUNC & 0x1F) + ((DATA_FUNC & 0x1F)<<6)
client.runPrint(f"pdcCfg -a OUTD -r 0x{OUTD:04x} -g")
PDC_SETTING.OUTD = OUTD

# === OUTF REGISTER ===
print("\n=== OUTF REGISTER ===")
if analogOnly:
    # in analog only mode, use FLAG pin as FLAG fonction
    FLAG_FUNC = OUT_MUX.FLAG
else:
    # in digital mode, use FLAG pin as SUM_TH_GT function
    #FLAG_FUNC = OUT_MUX.FLAG
    #FLAG_FUNC = OUT_MUX.TRG
    FLAG_FUNC = OUT_MUX.SUM_GT # NOTE: use this setting to set flag output as sum greater than
    #FLAG_FUNC = OUT_MUX.VSS
    #FLAG_FUNC = OUT_MUX.VDD
OUTF = (FLAG_FUNC & 0x1F) + ((FLAG_FUNC & 0x1F)<<6)
client.runPrint(f"pdcCfg -a OUTF -r 0x{OUTF:04x} -g")
PDC_SETTING.OUTF = OUTF

# === TRGC REGISTER ===
print("\n=== TRGC REGISTER ===")
TRGC = 0x0000
client.runPrint(f"pdcCfg -a TRGC -r 0x{TRGC:04x} -g")  # disable trigger command
PDC_SETTING.TRGC = TRGC

# === DISABLE ALL THE PIXELS ===
print("\n=== DISABLE ALL THE PIXELS ===")
    # NOTE: pdcPix returns the PDC to acquisition mode,
    #       if mode is not specified
client.runPrint("pdcPix --dis --mode NONE")

# === VALIDATE CONFIGURATIONS ===
print("\n=== VALIDATE CONFIGURATIONS ===")
icp.validPdcCfg()

# === OUTC REGISTER ===
print("\n=== OUTC REGISTER ===")
    # disable configuration output last once configuration was validated
FLAG_CFG_FUNC = OUT_MUX.CLK_CS     # default function
#FLAG_CFG_FUNC = OUT_MUX.VSS        # disabled
DATA_CFG_FUNC = OUT_MUX.CFG_VALID  # default function
#DATA_CFG_FUNC = OUT_MUX.VSS        # disabled
OUTC = (DATA_CFG_FUNC & 0x1F) + ((FLAG_CFG_FUNC & 0x1F)<<6)
client.runPrint(f"pdcCfg -a OUTC -r 0x{OUTC:04x} -g")
PDC_SETTING.OUTC = OUTC

# print the settings of all the PDCs
print("\n=== PDC SETTINGS ===")
PDC_SETTING.print()


# ------------------------
# --- enable PDC SPADs ---
# ------------------------
sectionPrint("enable PDC SPADs")

# generate a pattern to select which SPADs to enable
pixEnMask = np.zeros((pixMap.TOP_NX_PIX, pixMap.TOP_NY_PIX), dtype=int) # init an empty mask to set later
maskCX = 34 # center position in X axis (from wirebond 1 to wirebond 32)
maskCY = 35 # center position in Y axis (from CMOS pads to 2D SPADs)
maskX = 41 # width in X axis (used to match scintillator size)
maskY = 41 # width in Y axis (used to match scintillator size)


# load TCR CSV file into a pandas dataframe
dfTcr = pd.read_csv(tcrFile, header=0, sep=';')

# print both to terminal and a file the enabled pixels statistics
pixelStatsFileName = os.path.splitext(dsumCsvFile)[0]+".txt"
os.makedirs(os.path.dirname(pixelStatsFileName), exist_ok=True)
pixelStatsFile = open(pixelStatsFileName, 'w')
defaultStdout = sys.stdout
sys.stdout = systemHelper.Tee(defaultStdout, pixelStatsFile)

# loop for each PDC
for iPdc in range(icp.nPdcMax):
    if ((icp.pdcEnUser >> iPdc) & 0x1) == 0x1:
        print(f"{fgColors.blue}PDC{iPdc}:{fgColors.endc}")
        # PDC is enabled, configure it
        try:
            # NOTE: specify here a pattern to place on the PDCs
            # square/rectangle is the default
            pixEnMask[maskCX-maskX//2:maskCX+maskX//2, maskCY-maskY//2:maskCY+maskY//2] = 1

            # checker pattern can help for crosstalk analysis
            checkerPattern=False
            if checkerPattern:
                pixEnMask[maskCX-maskX//2:maskCX+maskX//2, maskCY-maskY//2:maskCY+maskY//2] = 1
                pixEnChMask = np.zeros((pixMap.TOP_NX_PIX, pixMap.TOP_NY_PIX), dtype=int)
                pitch = 2
                pixEnChMask[pitch//2::pitch, ::pitch] = 1
                pixEnChMask[::pitch, pitch//2::pitch] = 1
                pixEnMask = np.logical_and(pixEnMask, pixEnChMask)

            # NOTE: here is an example to manually disable some pixels specific for each PDC
            # NOTE: User can use different masks for each PDC here
            """
            if iPdc == 0:
                #pixEnMask[maskCX-maskX//2:maskCX+maskX//2, maskCY-maskY//2:maskCY+maskY//2] = 0
                pixEnMask[43, 41] = 0

            if iPdc == 1:
                #pixEnMask[maskCX-maskX//2:maskCX+maskX//2, maskCY-maskY//2:maskCY+maskY//2] = 0
                pixEnMask[22, 53] = 0

            if iPdc == 2:
                #pixEnMask[maskCX-maskX//2:maskCX+maskX//2, maskCY-maskY//2:maskCY+maskY//2] = 0
                pixEnMask[32, 32] = 0

            if iPdc == 3:
                #pixEnMask[maskCX-maskX//2:maskCX+maskX//2, maskCY-maskY//2:maskCY+maskY//2] = 0
                pixEnMask[54, 27] = 0
            """
            # NOTE: convertPixArrayToReg function from python module pdcSpadFunctions
            # supported methods: constant, average, percent, medianFactor, medianToMin
            # Here, keeping only SPADs with TCR below 100 cps (thConst)
            regs = pdcSpadFunctions.convertPixArrayToReg(
                        pixArray=dfTcr[f"SPAD_TCR{iPdc}"],
                        thMethod=pdcSpadFunctions.ThreshMethod.percent,
                        thConst=100.0,
                        thPct=95,
                        thOp=pdcSpadFunctions.ThreshOp.le,
                        pixEnMask = pixEnMask,
                        returnAnalysis=False,
                        log=False,  # set log to True to print the values of the registers to program
                        plot=False) # set plot to True to show tha map of the enabled pixels

            # set plotMask to True to see the pixel mask (pixEnMask) not considering the TCR
            plotMask = False
            if plotMask:
                # Create a custom colormap from green to red
                # You can define the colors at specific points along the colormap
                colors = [(0, 'black'), (1, 'white')]
                cmap = mcolors.LinearSegmentedColormap.from_list("BlackWhite", colors)
                plt.ion()
                plt.figure()
                plt.imshow(pixEnMask.T, cmap=cmap, origin='lower', aspect="equal")
                plt.show()

            # enable the proper SPADs
            icp.cfgAllPixRegs(regs, iPdc)
        except KeyError:
            print(f"{fgColors.red}ERROR: no data found for this PDC.{fgColors.endc}")
            # disable all pixels when no data found for a PDC
            cmd = f"pdcSpad --dis --spdc {iPdc} --mode NONE"
            client.runPrint(cmd, printCmd=False)

# Restore original stdout and close the log file when done
sys.stdout = defaultStdout
pixelStatsFile.close()

# done with all the PDCs, make sure to return the configuration to all PDCs
client.runPrint(f"ctlCfg -a CFGS -r 0x0000 -g") # disable single configuration


# ---------------------------------------
# --- return PDCs to acquisition mode ---
# ---------------------------------------
sectionPrint("return PDCs to acquisition mode")
client.runPrint("ctlCmd -c MODE_ACQ")

#
# ---------------------------------------------------------
# --- configure streamManager to maximize the bandwidth ---
# ---------------------------------------------------------
sectionPrint("configure streamManager to maximize the bandwidth")
TLAST_THRESHOLD = 500
    # number of packets (each ended with a tlast) to store
    # before transfering from the PL to the PS
TIMEOUT_STR = "-T 1 --sec"
    # this sets the timeout to 1 second.
    # After 1 second, if the threshold is not reached,
    # a transfer will occur from the PL to the PS.
client.runPrint(f"streamManager --thresh {TLAST_THRESHOLD} {TIMEOUT_STR}")

# ---------------------------------------------------
# ---  configure the Controller sum TH interface  ---
# ---------------------------------------------------
sectionPrint("configure the Controller sum TH interface")
client.runPrint(f"ctlCfg -a ITF2 -r 0b00000001 -g")

# ---------------------------------------------------
# --- configure the Controller coincidence module ---
# ---------------------------------------------------
sectionPrint("configure the Controller coincidence module")
# configure AUX0 out (OUT1 on PCB) to COINC_OK
#client.runPrint(f"auxOut --channel 0 --func 5 -g")

# select which PDCs to include into the coincidence
client.runPrint(f"ctlCfg -a COI0 -r 0x{icp.pdcEnUser&0xFFFF:04x} -g")
client.runPrint(f"ctlCfg -a COI1 -r 0x{(icp.pdcEnUser>>16)&0xFFFF:04x} -g")

# set the width of the coincidence windows (in clock cycle period)
COIN_WLEN = 1   # coincidence windows of X clock cycles # NOTE match with slow period ?
client.runPrint(f"ctlCfg -a COIW -r 0x{COIN_WLEN&0x03FF:04x} -g")

# set the coincidence thresholds
NCH_TH = 1      # number of PDC channel for a coincidence
NUM_BANK = 4096    # number of hits per PDC
cothReg = ((NCH_TH&0x7F)<<8) | (NUM_BANK&0x7)
client.runPrint(f"ctlCfg -a COTH -r 0x{cothReg:04x} -g")


# ---------------------------------------
# --- configure strobe timer
# ---------------------------------------
sectionPrint("configure strobe timer")
client.runPrint(f"set-ctl-tmr-prd FREQ=10e3")

# ---------------------------------------
# --- configure auxiliary ios
# ---------------------------------------
sectionPrint("configure auxiliary ios")
#client.runPrint("auxOut --ch 1 --func 0x5 -g") # COINC_OK
#client.runPrint("auxOut --ch 1 --func 0x8 -g") # SUM_TH_OR
##client.runPrint("auxOut --ch 1 --func 0x3 -n 0 -g") # AUX_IN0
client.runPrint("auxOut --ch 0 --func 23 -n 0 -g") # STRB_TMR
client.runPrint("auxOut --ch 1 --func 23 -n 0 -g") # STRB_TMR
#client.runPrint("auxIn --func 0 --channel 4 -g") # FSM_EXT1 = AUXI0
client.runPrint("auxIn --func 0 --channel 2 -g") # FSM_EXT1 = STRB_TMR

# ---------------------------------------
# --- notify user of manual steps
# ---------------------------------------
try:
    print(f"{fgColors.bYellow}Apply HV here{fgColors.endc}")
    input("Press [enter] key to continue")
except KeyboardInterrupt:
    print("\nKeyboard Interrupt: exit program")
    sys.exit()


# ------------------------------------------------
# --- start Controller FSM acquisition
# ------------------------------------------------
sectionPrint("start Controller FSM acquisition")
if not analogOnly:
    client.runPrint(f"ctlCfg -a FSMM -r 0x{fsmmReg|0x3:04x} -g"); # starts the FSM


# get the total execution time of the test
test_start_time = time.time()

# ------------------------
# --- ready to operate ---
# ------------------------
print("\n=== READY TO OPERATE ===")
# NOTE: Implement here a specific routine
try:
    input("Press [enter] key to exit")
except KeyboardInterrupt:
    print("\nKeyboard Interrupt: exit program")

finally:
    client.runPrint("stop")
    # total execution time
    test_stop_time = time.time()
    test_duration_sec = test_stop_time-test_start_time
    print(f"{fgColors.bBlue}Test duration:\n  {test_duration_sec:.3f} seconds \n  {test_duration_sec/60:.3f} min \n  {test_duration_sec/3600:.3f} hours{fgColors.endc}")
    # WARNING remove empty file at the end of the execution
    if "dsumCsvFile" in locals() and os.path.exists(dsumCsvFile) and os.path.getsize(dsumCsvFile) == 0:
        os.remove(dsumCsvFile)
        if "pixelStatsFileName" in locals() and os.path.exists(pixelStatsFileName):
            # remove pixel stats if file dsumCsvFile is empty
            os.remove(pixelStatsFileName)
    zynq.close()
    sys.exit()





