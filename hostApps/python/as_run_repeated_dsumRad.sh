#!bin/bash
# run like : as_run_repeated_dsumRad.sh <INTEGER_NUMBER_OF_REPETITIONS>

if [[ -n "$1" ]]
then
	REPEATED_MEASUREMENTS="$1"
else
	read -p "Enter the number of repititions: " REPEATED_MEASUREMENTS
fi


for i in $( seq 1 $REPEATED_MEASUREMENTS ); do
	python as_getDsumRad_noInputs.py;
done