#!/bin/bash

#OUTPUT FILE/PLOT PARAMS
export SPAD_BIAS_V=25
export HEAD_ID=62
export ANALOG_ONLY=False
export FNAME=""
export SAVE_PLOT=True
export MEAS_TIME=0.2

#PDC SETTING PARAMS
export PDC_EN="0xF"
export N_SPAD=4096
export DATA_TYPE="ZPP" #DSUM, ZPP

#TIMING PARAMS
export HOLD_TIME_NS=150.0
export RECH_TIME=10.0
export FLAG_TIME=2.0

#COINCIDENCE PARAMS
export COIN_WLEN=1
export COIN_NCH_TH=1
export COIN_NUM_BANK=1

#SCREAMER ID PARAMS
export SCREAMER_METHOD=average #average, threshold, percent, medianFactor, medianToMin
export SCREAMER_THRESHOLD=10000.0 #for "threshold" OR for getDsumRad (only uses threshold method)
export SCREAMER_PERCENT=90.0 #for "percent"
export SCREAMER_FACTOR=1.5 #for "medianFactor"

#SCRIPT SPECIFIC PARAMS
export RAD_SOURCE_NAME="test"
export TCR_FILE="/home/i8x/PDCv2-data/TCR/as_getSpadTcrUsingFlag/20260612_14h58m17_TCR_H54_200ms_25V_PDC0_PDC1_PDC2_PDC3.csv"
#export TCR_FILE="/home/i8x/PDCv2-data/TCR/getSpadTcrUsingFlag/20260226_15h43m23_TCR_H55_500ms_25V_PDC0_PDC1_PDC2_PDC3.csv"
