

from runModules.codeSetup import *

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




#check object
print(', '.join("%s: %s" % item for item in vars(setupObj).items()))