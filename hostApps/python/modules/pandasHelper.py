#----------------------------------------------------------------------------------
#-- Company: GRAMS
#-- Designer: Tommy Rossignol
#--
#-- Create Date: 2025-12-10
#-- Description:
#--     helpers with pandas analysis
#--
#-- Dependencies:
#-- Revision:
#-- Revision 1.0 - File Created
#-- Additional Comments:
#----------------------------------------------------------------------------------
import pandas as pd

#----------------------------------------------------------------------------------
# pandas DataFrame functions
#----------------------------------------------------------------------------------
# print available columns and memory usage
def print_df_usage(dataFrame):
    #print(list(dataFrame))
    dtypeColumns = dataFrame.dtypes
    dtypeColumns = pd.concat([pd.Series([dataFrame.index.dtype], index=["Index"]), dtypeColumns])
    memUsageColumns = dataFrame.memory_usage(index=True)/(1024.0*1024.0) # MB
    dfSpecs = pd.concat([dtypeColumns, memUsageColumns], axis=1, keys=["dtype", "memory (MB)"])
    print(dfSpecs)
    print(f"{memUsageColumns.sum()/1024.0:0.3f} GB")
    del dtypeColumns, memUsageColumns, dfSpecs
    #print(dtypeColumns)

# using a mask, show the statistics of the number of kept, total, %kept, %dropped data
def filterStats(keepMask, description="", doPrint=True):
    nKeep = keepMask.sum()
    nTotal = len(keepMask)
    pctKeep = 100.0*nKeep/nTotal
    pctDrop = 100.0-pctKeep
    if doPrint:
        print(f"Keeping {nKeep} values over {nTotal} ({pctKeep:.1f} %). Droping {pctDrop:.1f} %", end="")
        if description != "":
            print(f" ({description})")
        else:
            print("")
    return nKeep, nTotal, pctKeep, pctDrop
