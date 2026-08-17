import datetime, sys
from modules.fgColors import fgColors
import modules.pixMap as pixMap

# --------------------------------------------------
# --- Function to get ZPP of each PDC from h5 file
# --------------------------------------------------
def waitForH5File(timeOutSec=10):
    """
    function to wait for a new HDF5 file
    """
    t0 = datetime.datetime.now()
    while 1:
        db = h5Reader(deleteAfter=True,
                      hfAbsPath=zynq.h5Path,
                      hfFile="")

        if db.newFileReady():
            return db
        else:
            if datetime.datetime.now()-t0 > datetime.timedelta(seconds=timeOutSec):
                print(f"{fgColors.red}ERROR: Timeout while waiting for HDF5 data ({timeOutSec} seconds){fgColors.endc}")
                sys.exit()

def measCntRate(measTime, numPdc,
                spadRow=None, spadCol=None,
                spadIndex=None):
    """
    measCntRate: send cmd and cfg to Controller and PDC to get the SPAD count rate
    1- enable spads based on 64 bits pattern given and return to acquisition mode
    2- reset ZPP module
    3- wait for measTime for stats to build up
    4- send a Controller data packet with ZPP data
    5- wait to receive the file, fetch the ZPP data and close the file
    """
    if type(spadRow) != type(None) and type(spadRow) != type(None):
        # pixel/SPAD specified by row and column
        client.runPrint(f"pdcPix --dis --row {spadRow} --col {spadCol} --mode NONE")

    elif type(spadIndex) != type(None):
        # pixel/SPAD specified by index
        client.runPrint(f"pdcPix --dis --index {spadIndex} --mode NONE")

    else:
        # do not change pixels settings
        print("pixel/SPAD not specified")


    # send commands to the Controller/PDC
    client.runPrint("ctlCmd -c MODE_ACQ; " \
                    "ctlCmd -c RSTN_ZPP; " \
                    f"sleep {measTime:.06f}; " \
                    "ctlCmd -c PACK_TRG_A; ")

    # wait for the HDF5 result file
    db = waitForH5File()
    db.h5Open()

    # get ZPP results
    AVG = [-1]*numPdc
    for iPdc in range(0, numPdc):
        ZPP = db.getPdcZPP(iPdc=iPdc, zppSingle=PDC_ZPP_ITEM.AVG)
        if (ZPP != None) and (ZPP.AVG != -1):
            # use ZPP value and normalize count rate to 1 sec to have cps
            AVG[iPdc] = ZPP.AVG / measTime
            print(f"  PDC {iPdc} TCR = {AVG[iPdc]}")

    db.h5Close()

    return AVG


# ---------------------------------------
# --- data acquisition as a thread
# ---------------------------------------
class LoopingMethod(IntEnum):
    rows_cols=1,
    index=2,

def test_all_pixels(tp: tcrPlotter, update=False, numPdc=icp.nPdcMax):
    if type(tp) == type(None):
        print(f"{fgColors.red}tcrPlotter object must be created first{fgColors.endc}")
        sys.exit()

    if tp.nPdcMax != numPdc:
        numPdc = tp.nPdcMax

    # acquire data
    if icp.nSpad == 4096:
        # testing a full array
        rows = range(0, 64)
        cols = range(0, 64)
        #rows = range(63, -1, -1) # starts from the end
        #cols = range(63, -1, -1) # starts from the end
        lMethod = LoopingMethod.rows_cols

    elif icp.nSpad == 64:
        # testing only 2D CMOS SPADs
        rows = [63]
        cols = range(0, icp.nSpad)
        lMethod = LoopingMethod.rows_cols
    else:
        # testing based on index (fallback case, should not be used)
        rows = [0]
        cols = range(0, icp.nSpad)
        lMethod = LoopingMethod.index

    for iRow in rows:
        for iCol in cols:
            if not tp.run:
                break
            if lMethod == LoopingMethod.index:
                spadIndex = iCol
                spadRow = None
                spadCol = None
                iSpad = spadIndex
            elif lMethod == LoopingMethod.rows_cols:
                spadIndex = None
                spadRow = iRow
                spadCol = iCol
                if icp.nSpad == 4096:
                    iSpad = pixMap.idx_map(y=spadRow, x=spadCol)
                else:
                    iSpad = iCol
            # enable the proper pixels and measure the total countrate
            AVG_TCR = measCntRate(spadIndex=spadIndex,
                                  spadRow=spadRow,
                                  spadCol=spadCol,
                                  measTime=measTime,
                                  numPdc=numPdc)

            # add new data for each PDC
            for iPdc in range(0, numPdc):
                # put new data into data object
                tp.newData(iPdc=iPdc,
                           iSpad=iSpad,
                           avg=AVG_TCR[iPdc])

            tp.current_pixel_index = iSpad

            if update:
                tp.updatePlot()

    # indicate all tests are completed
    tp.done_test_all_pixels = True