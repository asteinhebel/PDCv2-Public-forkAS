#----------------------------------------------------------------------------------
#-- Company: GRAMS
#-- Designer: Tommy Rossignol
#-- 
#-- Create Date: 2023-03-22
#-- Description:
#--     module related to the system and error management
#--
#-- Dependencies: 
#-- Revision:
#-- Revision 1.0 - File Created
#-- Additional Comments:
#----------------------------------------------------------------------------------
import subprocess
import sys
import traceback
from modules.fgColors import fgColors

def sectionPrint(msg):
    print('')
    print('-'*50)
    print(f"--- {msg}")
    print('-'*50)

def printException(ex):
    # Get current system exception
    # https://docs.python.org/3/library/traceback.html
    ex_type, ex_value, ex_traceback = sys.exc_info()

    # Extract unformatter stack traces as tuples
    trace_back = traceback.extract_tb(ex_traceback)

    # print the message
    print(f"{fgColors.red}{'-'*75}{fgColors.endc}")
    print(f"{fgColors.red}{ex_type.__name__}{fgColors.endc}{' '*33}Traceback (most recent call last)")
    #print(f"{stack_trace}")
    for trace in trace_back:
        tfile=trace[0]
        tline=trace[1]
        tfunc=trace[2]
        tmess=trace[3]
        print(f"File {fgColors.bGreen}{tfile}:{tline}{fgColors.endc}")
        print(f"    Func.Name : {tfunc}")
        print(f"    Message : {tmess}")
    print(f"{fgColors.red}{ex_type.__name__}{fgColors.endc}: {ex_value}")

def get_gitVersion():
    try:
        version = '0x0000000'
        process = subprocess.Popen(['git', 'rev-parse', '--short', 'HEAD'],
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if not stderr:
            # no error
            version = '0x' + stdout.decode().strip('\n')
            tmp = int(version, base=16)
        return version
    except BaseException as ex:
        return version

