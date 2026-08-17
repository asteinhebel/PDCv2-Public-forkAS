import pyvisa, sys, glob, os, time
from modules.fgColors import fgColors

def checkSystHealth(healthbytes, instrument):

    if healthbytes=="000000" or healthbytes=="000001":
        #Require all systems off, including over-current and over-voltage protections
        print('All good')
        return 
    elif healthbytes=="000009": #output is on
        print("Output is on - turn off to begin")
        instrument.write("OUT OFF")
        return
    elif healthbytes=="0000A1": #ovp and ovc still set
        print("Over-current and -voltage protections are set - turn off to begin")
        instrument.write("PROT:OCP OFF")
        instrument.write("PROT:OVP OFF")
        return
    elif healthbytes=="0000A9": #output, ovp and ovc still set
        print("Settings remain from previous run - turn off to begin")
        instrument.write("OUT OFF")
        instrument.write("PROT:OCP OFF")
        instrument.write("PROT:OVP OFF")
        return
    else:
        print(f"Failed system health: {healthbytes}")
        sys.exit()

def powerRampDown(inst_dcps, inst_bkps):
    # ORNL SPECIFIC - turn off power 
    print("Turn off power supplies")
    print("Keysight ramping down....")
    for vDown in range(25,-1,-1):
        time.sleep(0.5)
        inst_dcps.write(f"APPL {float(vDown)}, 0.1")
    print("Keysight HV is turned off")

    print(f"Run over - turn off Keysight output")
    inst_dcps.write("APPL 0.0, 0.0")
    inst_dcps.write("OUTP OFF")
    inst_dcps.close()
    print("BK ramping down....")
    for vDown in range(250,-1,-10):
        time.sleep(0.5)
        inst_bkps.write(f"SOUR:VOLT {float(vDown)}")
    print("BK HV is turned off")

    print(f"Run over - turn off BK output")
    inst_bkps.write("SOUR:VOLT 0.0")
    inst_bkps.write("OUT OFF")
    inst_bkps.write("PROT:OCP OFF")
    inst_bkps.write("PROT:OVP OFF")
    inst_bkps.close()
    return

def setupHV(setupObj):
    ###Identify supplies
    rm = pyvisa.ResourceManager('@py')
    serialID_bkps = glob.glob('/dev/serial/by-id/*CP2102*')[0]
    serialID_dcps = glob.glob('/dev/serial/by-id/*Prolific*')[0]
    resource_bkps = f"ASRL/dev/{os.readlink(serialID_bkps)[-7:]}::INSTR"
    resource_dcps = f"ASRL/dev/{os.readlink(serialID_dcps)[-7:]}::INSTR"

    foundResources = rm.list_resources()
    if resource_bkps not in foundResources:
        print('Could not find BK device in expected location - try unplugging')
        sys.exit()
    if resource_dcps not in foundResources:
        print('Could not find Keysight device in expected location - try unplugging')
        sys.exit()

    ###Set up BK Precision to supply substrate voltage
    print('opening resource: ' + resource_bkps)
    inst_bkps = rm.open_resource(resource_bkps, write_termination = '\n',read_termination='\r\n',baud_rate=57600)
    print(f'BKPS Device: {inst_bkps.query("*IDN?")}')

    #clear status and errors
    systHealth = inst_bkps.query("STATUS?")[:6]
    checkSystHealth(systHealth, inst_bkps)

    #Set over current and over voltage protections
    inst_bkps.write("PRO:OCP:LEV 0.00005") #50 uA current limit
    inst_bkps.write("PRO:OVP:LEV 260.0") #260 V voltage limit
    print('Set over-current and over-voltage protection limits')
    inst_bkps.write("PROT:OCP ON")
    inst_bkps.write("PROT:OVP ON")

    ###Set up Keysight to supply bias voltage
    print('opening resource: ' + resource_dcps)
    inst_dcps = rm.open_resource(resource_dcps, write_termination = '\n',read_termination='\n',baud_rate=9600)
    print(f'DCPS Device: {inst_dcps.query("*IDN?")}')
    systHealth = inst_dcps.query("*TST?")
    if '0' in systHealth:
        print('All good')
        inst_dcps.write("SYST:BEEP")
    else:
        print(f"Failed system health: {systHealth}")
        sys.exit()

    #Ensure that the output is not on while setting proper values
    if int(inst_dcps.query("OUTP?")) == 1:
        print("output is on - turn it off")
        inst_dcps.write("OUTP OFF")
    #Set over current and over voltage protections
    inst_dcps.write(f"VOLT:PROT {setupObj.spadBias+0.5}")


    setupObj.inst_bkps = inst_bkps
    setupObj.inst_dcps = inst_dcps

    return inst_dcps, inst_bkps

def deliverBias(setupObj):
    try:
        print(f"{fgColors.bYellow}Apply HV here{fgColors.endc}")
        if os.environ.get("BATCH_MODE") is None:
            #ORNL SPECIFIC - turn on power supplies
            input("Apply substrate bias?")
            #Ramp output voltage to 250 in 10V steps
            setupObj.inst_bkps.write("OUT ON")

            print("BK HV Ramping up....")
            for vUp in range(10, 250+1, 10):
                time.sleep(0.5)
                if vUp==10: #begin ramp slowly
                    for vUpFine in range(10):
                        setupObj.inst_bkps.write(f"SOUR:VOLT {float(vUpFine)}")
                        time.sleep(0.5)
                setupObj.inst_bkps.write(f"SOUR:VOLT {float(vUp)}")
            print("BK HV is supplied")

            input("Apply SPAD bias?")
            #Ramp output voltage to 25 in 1V steps
            setupObj.inst_dcps.write("APPL 0.0, 0.0")
            setupObj.inst_dcps.write("OUTP ON")

            print("Keysight HV Ramping up....")
            for vUp in range(int(setupObj.spadBias)+1):
                time.sleep(0.5)
                setupObj.inst_dcps.write(f"APPL {float(vUp)}, 0.1")
            print("Keysight HV is supplied")

            input("Press [enter] key to continue")
    except KeyboardInterrupt:
        print("\nKeyboard Interrupt: exit program")
        powerRampDown()
        sys.exit()