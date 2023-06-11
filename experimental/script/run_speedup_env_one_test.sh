#!bin/bash

## Define here what you need
SOURCE_FILE="base-one-test"
SIMULATION_NAME="ex6"
MIN_PROC=64
MAX_PROC=64
#MAX_PROC=4096

## work

F()
{
	for (( solv = 1; solv <= 9; solv++ ))
	do
		for (( pre = 0; pre <= 6; pre++ ))
		do
			cp $SOURCE_FILE ${SOURCE_FILE}-${1}-${solv}-${pre}.bash
			FILE=${SOURCE_FILE}-${1}-${solv}-${pre}.bash
			ccc_msub ${FILE}
		done
	done
}

for (( k = $MIN_PROC ; k <= $MAX_PROC; k*=2 ))
do
	F $k
done
