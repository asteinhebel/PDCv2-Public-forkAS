#----------------------------------------------------------------------------------
#-- Company: GRAMS
#-- Designer: Tommy Rossignol
#--
#-- Create Date: 2023-03-22
#-- Description:
#--     module related to the ZYNQ settings and functions
#--
#-- Dependencies:
#-- Revision:
#-- Revision 1.0 - File Created
#-- Additional Comments:
#----------------------------------------------------------------------------------
import os, sys
from enum import IntEnum

# custom modules
from modules.fgColors import fgColors
import modules.systemHelper as systemHelper
import modules.sshClientHelper as sshClientHelper
#from modules.zynqHelper import *

try:
    moduleName = os.path.basename(__file__)
except NameError:
    moduleName = "moduleNameNotFound.py"


# -----------------------------------------------
# --- Convert board_info into a dict for export
# -----------------------------------------------
def boardInfo2Dict(client):
    """
    Running boardInfo app on zynq, list all parameters
    and returning them as a python dictionary.
    """
    biStr = client.runReturnStr('boardInfo -l', printCmd=False)
    biDict=dict((key.strip(), value.strip())
        for key, value in (line.split(':')
        for line in biStr))
    return biDict

# -----------------------------------------------
# --- Specifications on hardware used
# -----------------------------------------------
class hardware_info:
    def __init__(self, headID, adaptID, cableLen, description=""):
        """
        Used to store informations about a test setup
        """
        self.headID=headID
        self.adaptID=adaptID
        self.cableLen=cableLen
        self.description=description

# -----------------------------------------------
# --- testing functions and modules
# -----------------------------------------------
if __name__ == "__main__":
    # testing h5helper class only when module is not called via 'import'
    try:
        # open a client
        client = sshClientHelper.sshClient(host="zcudev")

        # test boardInfo2Dict
        BI = boardInfo2Dict(client)

    except BaseException as ex:
        ## Get current system exception
        systemHelper.printException(ex)

    finally:
        # executed no matter what
        del client









