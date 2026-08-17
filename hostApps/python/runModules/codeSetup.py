import os
import time
from modules.pdcHelper import *

class Descriptor():
    def __init__(self, typeIn:type):
        self.typeIn = typeIn
    def __set_name__(self, owner, name):
        self._name = name

    def __get__(self, instance, owner):
        return instance.__dict__[self._name]

    def __set__(self, instance, value):
        if self.typeIn==str and value==None:
            instance.__dict__[self._name] = ""
            return
        elif self.typeIn==float and value==None:
            instance.__dict__[self._name] = float('nan')
            return
            
        try:
            instance.__dict__[self._name] = self.typeIn(value)
        except ValueError:
            raise ValueError(f'"{self._name}" must be a {self.typeIn}') from None


class setupValues:
    """Class to carry input values"""
    scriptName = Descriptor(str)
    startTime = Descriptor(float)
    measTime = Descriptor(float)
    spadBias = Descriptor(float)
    spadBiasStr = Descriptor(str)
    headId = Descriptor(float)
    headIdStr = Descriptor(str)
    fname = Descriptor(str)
    pdcEn = Descriptor(int)
    sysClkPrd = Descriptor(float)
    nspad = Descriptor(int)

def setBias(setupObj):
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
    else:
        spadBias = None
        spadBiasStr = ""
    setupObj.spadBias = spadBias
    setupObj.spadBiasStr = spadBiasStr

def setHeadID(setupObj):
    #NOTE: set environment variable HEAD_ID to store it in data file name only set the integer value (e.g. 48)
    headStr = ""
    if os.environ.get("HEAD_ID") is not None:
        headId = os.environ['HEAD_ID']
        headStr = f"H{headId}_"
        print(f"Using head {headId}")
    else:
        headId = None
        headStr = ""
    setupObj.headId = headId
    setupObj.headStr = headStr
        
def setName(setupObj, pathIn):
    try:
        setupObj.scriptName = pathIn
    except NameError:
        setupObj.scriptName = "fileNameNotFound.py"

    # get the total execution time of the test
    setupObj.startTime = time.time()

def createPDCSetting(setupObj):
    setupObj.pdcSetting = pdc_setting()