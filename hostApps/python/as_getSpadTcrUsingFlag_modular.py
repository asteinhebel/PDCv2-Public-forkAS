

from runModules.codeSetup import *
from runModules.zynqSetup import *
from runModules.icpSetup import *
from runModules.hvSetup import * 
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
setSysClkPrd(setupObj.icp) 

# -----------------------------------------------
# --- reset of the controller
# -----------------------------------------------
resetController(setupObj.icp) 

# -----------------------------------------------
# --- configure controller packet
# -----------------------------------------------
configControllerPacket(setupObj.icp) #add option for ZPP or DSUM

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
prepController(setupObj.icp) 

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
setConfigMode(setupObj.client)

configAllRegisters(setupObj) #add option to change timing variables

disableAllPixels(setupObj.client)
validateConfig(setupObj.icp)

configRegister_outc(setupObj) 
setupObj.pdcSetting.print() 

# ---------------------------------------
# --- configure Controller ZPP module ---
# ---------------------------------------
sectionPrint("configure Controller ZPP module")
configZppMethod(setupObj) 

# ---------------------------------------
# --- ORNL SPECIFIC - Set up power supply settings
# ---------------------------------------
sectionPrint("configure HV")
inst_dcps, inst_bkps = setupHV(setupObj) 

# ---------------------------------------
# --- Notify user of manual steps
# ---------------------------------------
deliverBias(setupObj) 

# ------------------------------------------------
# --- Prepare Controller FSM for the acquisition
# ------------------------------------------------
prepFSM(setupObj.client)

#check object
print(', '.join("%s: %s" % item for item in vars(setupObj).items()))

powerRampDown(inst_dcps, inst_bkps)