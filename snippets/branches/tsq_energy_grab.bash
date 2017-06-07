#!/bin/bash

list=`ls *tsq*`
for item in $list
do
	energy=`sed -n 2p $item`
	echo "$item $energy" >> tsqfile.txt
	#echo $item
done

#list2=`ls *.log`
#for item2 in $list2
#d#o
#	#echo $item2
#	DFTgrab=`grep 'xyz read in from DFT' $item2`
#	echo $item2 $DFTgrab \n >> test.log
#	i=1
# 	for item3 in $DFTgrab
#	do
#	  if [ $i -eq 9 ] 
#		echo $item3
#	  fi
#	let i=i+1
#	done
#	#echo $DFTgrab 
#done
