#!/bin/bash
rm *eps *png
mkdir Plots031023
rm -r Plots031023/Between_Models_Comparison 
mkdir Plots031023/Between_Models_Comparison

for energystat in 250 both three
do
    for pid in dNdx
    do
	for error in Stat 
	do
	    root -q test_models.C\(\"${energystat}\",\"${error}\",\"${pid}\",\"paper\"\)
        done
    done
done

mv *png Plots031023/Between_Models_Comparison/.
mv *eps Plots031023/Between_Models_Comparison/.
