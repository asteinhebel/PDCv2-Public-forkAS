import os
from modules.zynqCtlPdcRoutines import initCtlPdcFromClient, packetBank 


def setPDCController(setupObj, sysClkPrd=10e-9):
    setupObj.sysClkPrd = sysClkPrd
    # NOTE: select the PDC to use:
    #       pdcEn=0x1 -> PDC0
    #       pdcEn=0x2 -> PDC1
    #       pdcEn=0x4 -> PDC2
    #       pdcEn=0x8 -> PDC3
    #       pdcEn=0xF -> PDC0, PDC1, PDC2, PDC3
    # NOTE: set environment variable PDC_EN tp set which PDCs to use
    setupObj.pdcEn = int(os.environ.get("PDC_EN", default="0xF"), 0)
    setupObj.icp = initCtlPdcFromClient(client=setupObj.client, sysClkPrd=setupObj.sysClkPrd, pdcEn=setupObj.pdcEn) 

def resetController(contObj):
    contObj.resetCtl()

def setSysClkPrd(contObj):
    contObj.setSysClkPrd()

def configControllerPacket(contObj):
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
    contObj.setCtlPacket(bank=packetBank.BANKA, SCS=SCSA, SCD=SCDA, SPD=SPDA)

def setDelay(contObj):
    contObj.setDelay(signal="CFG_DATA", delay=300)

def checkPowerGood(contObj):
    contObj.checkPowerGood()

def setCfgRtnEn(contObj):
    contObj.setCfgRtnEn()

def prepPDC(contObj):
    contObj.preparePDC()

def prepController(contObj):
    setDelay(contObj)
    checkPowerGood(contObj)
    setCfgRtnEn(contObj)
    prepPDC(contObj)

def setPixNmb(setupObj):
    setupObj.nspad = int(os.environ.get("N_SPAD", default=4096))#64))
    setupObj.icp.nSpad = setupObj.nspad # AS - TEST IN LAB

def validateConfig(contObj):
    # === VALIDATE CONFIGURATIONS ===
    print("\n=== VALIDATE CONFIGURATIONS ===")
    contObj.validPdcCfg()