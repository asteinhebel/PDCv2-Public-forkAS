#----------------------------------------------------------------------------------
#-- Company: GRAMS
#-- Designer: Tommy Rossignol
#--
#-- Create Date: 2025-02-28
#-- Description:
#--      Using python ssh libraries to send remote commands to the ZCU102
#--      This script prepare the PDCs (general setup).
#--      It is a starting point or a template to create new scripts.
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

# custom modules
from modules.fgColors import fgColors
from modules.zynqEnvHelper import PROJECT_PATH, HOST_APPS_PATH, USER_DATA_DIR, HDF5_DATA_DIR
import modules.sshClientHelper as sshClientHelper
import modules.systemHelper as systemHelper
from modules.zynqCtlPdcRoutines import initCtlPdcFromClient, packetBank
from modules.zynqDataTransfer import zynqDataTransfer
from modules.systemHelper import sectionPrint
from modules.pdcHelper import *
#from modules.zynqHelper import *
from modules.h5Reader import *

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
zynq.init()

# -----------------------------------------------
# --- prepare controller for acquisition
# -----------------------------------------------
# NOTE: select here the PDC to use:
#       pdcEn=0x1 -> PDC0
#       pdcEn=0x2 -> PDC1
#       pdcEn=0x4 -> PDC2
#       pdcEn=0x8 -> PDC3
#       pdcEn=0xF -> PDC0, PDC1, PDC2, PDC3
icp = initCtlPdcFromClient(client=client, sysClkPrd=10e-9, pdcEn=0xF)

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
SCSA = 0x0000
# configure CTL_DATA_A
SCDA = 0x0000
# configure PDC_DATA_A
    # 0x0100 = DSUM
    # 0x00F7 = ZPP
SPDA = 0x00F7
icp.setCtlPacket(bank=packetBank.BANKA, SCS=SCSA, SCD=SCDA, SPD=SPDA)

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
#icp.nSpad = 4096

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
EDGE_LVLN = 0
DIS_MEM = 0
PIXL = ((DIS_MEM<<13) + (EDGE_LVLN<<12) + (FLAG_EN<<8) + (TRG_DGTL_FEN<<4) + (ACTIVE_QC_EN<<1))
client.runPrint(f"pdcCfg -a PIXL -r {PIXL} -g")  # configure pixel register
PDC_SETTING.PIXL = PIXL

# === TIME REGISTER ===
print("\n=== TIME REGISTER ===")
HOLD_TIME = 150.0
RECH_TIME = 10.0
FLAG_TIME = 10.0
client.runPrint(f"pdcTime --hold {HOLD_TIME} --rech {RECH_TIME} --flag {FLAG_TIME} -g")
PDC_SETTING.TIME = client.runReturnSplitInt('pdcTime -g')


# === ANLG REGISTER ===
print("\n=== ANLG REGISTER ===")
ANLG = 0x0000; # disabled
#ANLG = 0x001F; # full amplitude (~30 µA)
client.runPrint(f"pdcCfg -a ANLG -r {ANLG} -g")  # set analog monitor
PDC_SETTING.ANLG = ANLG

# === XXXX REGISTER ===
# skipping registers STHH to DTXC

# === OUTD REGISTER ===
print("\n=== OUTD REGISTER ===")
#DATA_FUNC = OUT_MUX.FLAG
#DATA_FUNC = OUT_MUX.TRG
DATA_FUNC = OUT_MUX.PIX_QC
#DATA_FUNC = OUT_MUX.VSS
#DATA_FUNC = OUT_MUX.VDD
OUTD = (DATA_FUNC & 0x1F) + ((DATA_FUNC & 0x1F)<<6)
client.runPrint(f"pdcCfg -a OUTD -r 0x{OUTD:04x} -g")
PDC_SETTING.OUTD = OUTD

# === OUTF REGISTER ===
print("\n=== OUTF REGISTER ===")
FLAG_FUNC = OUT_MUX.FLAG
#FLAG_FUNC = OUT_MUX.TRG
#FLAG_FUNC = OUT_MUX.VSS
#FLAG_FUNC = OUT_MUX.VDD
OUTF = (FLAG_FUNC & 0x1F) + ((FLAG_FUNC & 0x1F)<<6)
client.runPrint(f"pdcCfg -a OUTF -r 0x{OUTF:04x} -g")
PDC_SETTING.OUTF = OUTF

# === TRGC REGISTER ===
print("\n=== TRGC REGISTER ===")
TRGC = 0x0000
client.runPrint(f"pdcCfg -a TRGC -r {TRGC} -g")  # disable trigger command
PDC_SETTING.TRGC = TRGC

# === DISABLE ALL THE PIXELS ===
print("\n=== DISABLE ALL THE PIXELS ===")
    # NOTE: pdcPix returns the PDC to acquisition mode,
    #       if mode is not specified
client.runPrint("pdcPix --dis --mode NONE")

# Enable some pixels
print("\n=== ENABLE SOME PIXELS ===")
# NOTE: see pdcPix app help to see available options
pixIndex = 0
client.runPrint(f"pdcPix --index {pixIndex} --mode NONE")

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

# ---------------------------------------
# --- return PDCs to acquisition mode ---
# ---------------------------------------
sectionPrint("return PDCs to acquisition mode")
client.runPrint(f"ctlCmd -c MODE_ACQ")

# ready to operate
print("\n=== READY TO OPERATE ===")
# NOTE: Implement here a specific routine
try:
    input("Press [enter] key to exit")
except KeyboardInterrupt:
    print("\nKeyboard Interrupt: exit program")
    sys.exit()





