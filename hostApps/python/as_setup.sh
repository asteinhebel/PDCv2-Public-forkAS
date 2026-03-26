#!/bin/bash

export SPAD_BIAS_V=25
export HEAD_ID=54
export ANALOG_ONLY=False
export RAD_SOURCE="test"
#export TCR_FILE="/home/i8x/PDCv2-data/TCR/getSpadTcrUsingFlag/20260302_17h42m27_TCR_H54_500ms_25V_PDC0_PDC1_PDC2_PDC3.csv"
export TCR_FILE="/home/i8x/PDCv2-data/TCR/getSpadTcrUsingFlag/20260226_15h43m23_TCR_H55_500ms_25V_PDC0_PDC1_PDC2_PDC3.csv"
export PDC_EN="0xF"
export HOLD_TIME_NS=5000.0
export RECH_TIME=10.0
export FLAG_TIME=2.0
export N_SPAD=4096
export SCREAMER_THRESHOLD=400.0
export FNAME="test_new"
export SAVE_PLOT=False
export DATA_TYPE="DSUM"
export COIN_WLEN=1
export COIN_NCH_TH=1
export COIN_NUM_BANK=1
export MEAS_TIME=0.5