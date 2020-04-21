#!/bin/bash
#set absorption length in mm
lambdaMin=200
lambdaMax=1000
lambdaStep=100
file="run.mac"
runprefix="run001"
g4simpleLoc="/mnt/geda00/patrick/software/g4simple_build/g4simple"

while [[ $lambdaMin -le $lambdaMax ]]
do
                sed -i "/lArAbsLength/c\\/optics\/lArAbsLength $lambdaMin mm" $file
        ofile="$runprefix-$lambdaMin.root"
        sed -i "/write\/filename/c\\/write\/filename $ofile" $file
        echo "Start simulation with absoprtion length of $lambdaMin mm and writting the results to $ofile"
	$g4simpleLoc $file > "$ofile.log" &
        sleep 1
        lambdaMin=$(( $lambdaMin + $lambdaStep ))
 done
