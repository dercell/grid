#!/bin/sh

m=0
T=1.79999
while [ $m -lt "20" ]
do 	
	T=$(echo "$T + 0.000001" | bc)
	let m++
	nohup ./ed $T>out.txt>&1 &

done


