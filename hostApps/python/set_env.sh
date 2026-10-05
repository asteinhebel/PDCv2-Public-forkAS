#!bin/bash

export SHOW_PLOT=False
export DO_PLOT=True
export RUN_VERBOSE=False
export SAVE_DATA=False
export ENERGY=dsumLin #dsum, dsumLevel, dsumLin

read -p "Input PLOT_DIR    " plot_dir
export PLOT_DIR=$plot_dir
echo $PLOT_DIR

read -p "Input FILE_IN    " -i "/" -e file_in
export FILE_IN=$file_in