#!bin/bash

## Define here what you need
SOURCE_FILE="base"
SIMULATION_NAME="ex6"
MAX_PROC=4096
MULT=8

## work

F()
{
	cp $SOURCE_FILE ${SOURCE_FILE}-${1}.bash
	FILE=${SOURCE_FILE}-${1}.bash
	sed -i -e "s/PROC/${1}/g" $FILE
}

for (( k = 1; k <= $MAX_PROC; k*=2 ))
do
	F $k
done
