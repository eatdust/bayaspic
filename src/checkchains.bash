#!/bin/bash

if [ $# -eq 0 ]; then
    echo 'Usage: checkchains.bash listofmodel.dat chaindir'
    exit
fi

modelfile=$1
chaindir=$2

listnames=`cat $modelfile | awk '{print $1$2}'`

for name in $listnames
do
    isthere=`ls $chaindir/*.ini | grep $name`
    if [[ -z "$isthere" ]]; then
	echo 'Missing model '$name
    else
	echo 'Model '$name
	echo 'Found file: '$isthere
    fi
    echo
done
	
