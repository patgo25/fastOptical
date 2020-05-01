#!/bin/bash
#wlsr radius in mm
lambdaMin=300
lambdaMax=1000
lambdaStep=100
file="run.mac"
runprefix="run0011"
g4simpleLoc="/mnt/geda00/patrick/software/g4simple_build/g4simple"

while [[ $lambdaMin -le $lambdaMax ]]
do
        sed -i "/wlsr\/radius/c\\/geometry\/wlsr\/radius $lambdaMin mm" $file
        sed -i "/generator\/SetRadiusMax/c\\/generator\/SetRadiusMax $lambdaMin mm" $file
	ofile="$runprefix-$lambdaMin.root"
        sed -i "/write\/filename/c\\/write\/filename $ofile" $file
        echo "Start simulation with WLSR radius of $lambdaMin mm and save the results into $ofile";
	$g4simpleLoc $file > "$ofile.log" &
        sleep 1
        lambdaMin=$(( $lambdaMin + $lambdaStep ))
 done
