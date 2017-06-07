#!/bin/bash

 list2=`ls *.inp.out`
  for item2 in $list2
   do
   dH=`grep 'Total Enthalpy' $item2`
   dS=`grep "Total Entropy" $item2`
   echo $item-$item2 $dH $dS "\n" >> Freqy.txt
   done
