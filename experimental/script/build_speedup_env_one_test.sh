#!bin/bash

## Define here what you need
SOURCE_FILE="base-one-test"
SIMULATION_NAME="ex6"
MAX_PROC=32
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
			sed -i -e "s/PROC/${1}/g" $FILE
			sed -i -e "s/SOLVER/${solv}/g" $FILE
			sed -i -e "s/PRECONDITIONER/${pre}/g" $FILE
		done
	done
}

for (( k = 32; k <= $MAX_PROC; k*=2 ))
do
	F $k
done
