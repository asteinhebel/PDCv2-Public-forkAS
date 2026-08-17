

from runModules.codeSetup import *
from runModules.zynqSetup import *
from runModules.icpSetup import *
#from runModules.hvSetup import * # AS - TEST IN LAB
from modules.systemHelper import sectionPrint

# -----------------------------------------------
# --- 
# ----------------------------------------------- 
setupObj = setupValues()

setName(setupObj,os.path.basename(__file__))

# -----------------------------------------------
# --- Global vars
# -----------------------------------------------
# database of the data
tp = None

setattr(setupObj, "measTime", float(os.environ.get("MEAS_TIME", default=0.2))) #second
setBias(setupObj)
setHeadID(setupObj)
setattr(setupObj, "fname", str(os.environ.get("FNAME")))


# -----------------------------------------------
# --- open a connection with the ZCU102 board
# -----------------------------------------------
sectionPrint("open a connection with the ZCU102 board")
openZynqConnection(setupObj)
# -----------------------------------------------
# --- prepare Zynq platform
# -----------------------------------------------
sectionPrint("prepare Zynq platform")
prepZynq(setupObj)
# -----------------------------------------------
# --- prepare controller for acquisition
# -----------------------------------------------
setPDCController(setupObj)


# -----------------------------------------------
# --- set system clock period
# -----------------------------------------------
#setSysClkPrd(setupObj.icp) # AS - TEST IN LAB
# -----------------------------------------------
# --- reset of the controller
# -----------------------------------------------
#resetZynq(setupObj.icp) # AS - TEST IN LAB


# -----------------------------------------------
# --- configure controller packet
# -----------------------------------------------
configControllerPacket(setupObj.icp)


# -----------------------------------------------
# --- set delay of CFG_DATA pins
# -----------------------------------------------
# -----------------------------------------------
# --- check for power good
# -----------------------------------------------
# -----------------------------------------------
# --- enable CFG_RTN_EN
# -----------------------------------------------
# -----------------------------------------------
# --- prepare PDC for configuration
# -----------------------------------------------
#prepController(setupObj.icp) # AS - TEST IN LAB

# -----------------------------------------------
# --- Testing all the pixels,
# --- 4096 for a complete 3D SPAD array
# --- 64 for the embedded 2D CMOS SPADs
# -----------------------------------------------
setPixNmb(setupObj)

# --------------------------
# --- configure the PDCs ---
# --------------------------
sectionPrint("configure the PDCs")
createPDCSetting(setupObj)
#setConfigMode(setupObj.client) # AS - TEST IN LAB

#configAllRegisters?
"""configRegister_pixl(setupObj)
configRegister_time(setupObj)
configRegister_anlg(setupObj)
configRegister_outd(setupObj)
configRegister_outf(setupObj)
configRegister_trgc(setupObj)"""

disableAllPixels(setupObj.client)
validateConfig(setupObj.icp)

#configRegister_outc(setupObj) # AS - TEST IN LAB
#setupObj.pdcSetting.print() # AS - TEST IN LAB

# ---------------------------------------
# --- configure Controller ZPP module ---
# ---------------------------------------
sectionPrint("configure Controller ZPP module")
#configZppMethod(setupObj) # AS - TEST IN LAB

# ---------------------------------------
# --- ORNL SPECIFIC - Set up power supply settings
# ---------------------------------------
sectionPrint("configure HV")
#setupHV(setupObj) # AS - TEST IN LAB

# ---------------------------------------
# --- Notify user of manual steps
# ---------------------------------------
#deliverBias(setupObj) # AS - TEST IN LAB

# ------------------------------------------------
# --- Prepare Controller FSM for the acquisition
# ------------------------------------------------
#prepFSM(setupObj.client)

#check object
print(', '.join("%s: %s" % item for item in vars(setupObj).items()))