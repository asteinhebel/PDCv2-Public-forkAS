import os
import modules.sshClientHelper as sshClientHelper 
from modules.zynqEnvHelper import PROJECT_PATH, HOST_APPS_PATH, USER_DATA_DIR, HDF5_DATA_DIR
from modules.zynqDataTransfer import zynqDataTransfer  
from modules.pdcHelper import *


def openZynqConnection(obj, hostCfgName="zcudev"):
    # open a client based on its name in the ssh config file
    obj.client = sshClientHelper.sshClientFromCfg(hostCfgName=hostCfgName) 

def prepZynq(obj):
    zynq = zynqDataTransfer(sshClientZynq=obj.client)  
    zynq.hexAppOutPathDefault = os.path.join(HDF5_DATA_DIR, os.path.splitext(obj.scriptName)[0]) 
    zynq.init()  

def setConfigMode(clientObj):
    # set PDCs to configuration mode
    clientObj.runPrint("ctlCmd -c MODE_CFG")  

def configAllRegisters(setupObj):
    configRegister_pixl(setupObj)
    configRegister_time(setupObj)
    configRegister_anlg(setupObj)
    configRegister_outd(setupObj)
    configRegister_outf(setupObj)
    configRegister_trgc(setupObj)

def configRegister_pixl(setupObj):
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
    setupObj.client.runPrint(f"pdcCfg -a PIXL -r {PIXL} -g")  # configure pixel register  
    setupObj.pdcSetting.PIXL = PIXL  

def configRegister_time(setupObj):
    # === TIME REGISTER ===
    print("\n=== TIME REGISTER ===")
    HOLD_TIME = float(os.environ.get("HOLD_TIME_NS", default=150.0))
    RECH_TIME = float(os.environ.get("RECH_TIME_NS", default=10.0))
    FLAG_TIME = float(os.environ.get("FLAG_TIME_NS", default=10.0))
    setupObj.client.runPrint(f"pdcTime --hold {HOLD_TIME} --rech {RECH_TIME} --flag {FLAG_TIME} -g") 
    setupObj.pdcSetting.TIME = setupObj.client.runReturnSplitInt('pdcTime -g') 

def configRegister_anlg(setupObj):
    # === ANLG REGISTER ===
    print("\n=== ANLG REGISTER ===")
    #ANLG = 0x0000; # disabled
    ANLG = 0x001F; # full amplitude (~30 µA)
    setupObj.client.runPrint(f"pdcCfg -a ANLG -r {ANLG} -g")  # set analog monitor
    setupObj.pdcSetting.ANLG = ANLG

def configRegister_outd(setupObj):
    # === OUTD REGISTER ===
    print("\n=== OUTD REGISTER ===")
    #DATA_FUNC = OUT_MUX.FLAG
    #DATA_FUNC = OUT_MUX.TRG
    #DATA_FUNC = OUT_MUX.PIX_QC
    DATA_FUNC = OUT_MUX.VSS
    #DATA_FUNC = OUT_MUX.VDD
    OUTD = (DATA_FUNC & 0x1F) + ((DATA_FUNC & 0x1F)<<6)
    setupObj.client.runPrint(f"pdcCfg -a OUTD -r 0x{OUTD:04x} -g")
    setupObj.pdcSetting.OUTD = OUTD

def configRegister_outf(setupObj):
    # === OUTF REGISTER ===
    print("\n=== OUTF REGISTER ===")
    FLAG_FUNC = OUT_MUX.FLAG
    #FLAG_FUNC = OUT_MUX.TRG
    #FLAG_FUNC = OUT_MUX.VSS
    #FLAG_FUNC = OUT_MUX.VDD
    OUTF = (FLAG_FUNC & 0x1F) + ((FLAG_FUNC & 0x1F)<<6)
    setupObj.client.runPrint(f"pdcCfg -a OUTF -r 0x{OUTF:04x} -g")
    setupObj.pdcSetting.OUTF = OUTF

def configRegister_trgc(setupObj):
    # === TRGC REGISTER ===
    print("\n=== TRGC REGISTER ===")
    TRGC = 0x0000
    setupObj.client.runPrint(f"pdcCfg -a TRGC -r {TRGC} -g")  # disable trigger command
    setupObj.pdcSetting.TRGC = TRGC

def configRegister_outc(setupObj):
    # === OUTC REGISTER ===
    print("\n=== OUTC REGISTER ===")
        # disable configuration output last once configuration was validated
    #FLAG_CFG_FUNC = OUT_MUX.CLK_CS     # default function
    FLAG_CFG_FUNC = OUT_MUX.VSS        # disabled
    #DATA_CFG_FUNC = OUT_MUX.CFG_VALID  # default function
    DATA_CFG_FUNC = OUT_MUX.VSS        # disabled
    OUTC = (DATA_CFG_FUNC & 0x1F) + ((FLAG_CFG_FUNC & 0x1F)<<6)
    setupObj.client.runPrint(f"pdcCfg -a OUTC -r 0x{OUTC:04x} -g")
    setupObj.pdcSetting.OUTC = OUTC

def disableAllPixels(clientObj):
    # === DISABLE ALL THE PIXELS ===
    print("\n=== DISABLE ALL THE PIXELS ===")
        # NOTE pdcPix returns the PDC to acquisition mode NOTE
    clientObj.runPrint("pdcPix --dis")

def configZppMethod(setupObj):
    class ZppModuleSetMethod(IntEnum):
        app=0,
        registers=1

    method = ZppModuleSetMethod.app
    if method == ZppModuleSetMethod.app:
        # new application available from 20250509 image
        setupObj.client.runPrint(f"set-ctl-zpp-prd {setupObj.measTime}")
    elif method == ZppModuleSetMethod.registers:
        print("  Configure ZPP Timer High Period")
        ZPP_HIGH_PRD=setupObj.measTime # seconds
        ZPP_HIGH_REG=int(ZPP_HIGH_PRD/setupObj.sysClkPrd)
        setupObj.client.runPrint(f"ctlCfg -a ZPH0 -r 0x{ZPP_HIGH_REG&0xFFFF:04x} -g")
        setupObj.client.runPrint(f"ctlCfg -a ZPH1 -r 0x{(ZPP_HIGH_REG>>16)&0xFFFF:04x} -g")

        print("  Configure ZPP Timer Low Period")
        ZPP_LOW_PRD=setupObj.sysClkPrd
        ZPP_LOW_REG=int(ZPP_LOW_PRD/setupObj.sysClkPrdd)
        setupObj.client.runPrint(f"ctlCfg -a ZPL0 -r 0x{ZPP_LOW_REG&0xFFFF:04x} -g")
        setupObj.client.runPrint(f"ctlCfg -a ZPL1 -r 0x{(ZPP_LOW_REG>>16)&0xFFFF|0x8000:04x} -g")  # |0x8000 to enable ZPP

def prepFSM(contObj):
    contObj.run(f"ctlCfg -a FSMM -r 0x0101 -g"); # triggered by a COMMAND