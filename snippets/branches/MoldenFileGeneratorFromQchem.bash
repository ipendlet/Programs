#!/bin/bash

list=`ls *.out`

for item in $list
do
a=`grep -n "MOLDEN-FORMATTED INPUT FILE FOLLOWS" $item`
b=${a%:*}
c=`grep -n "END OF MOLDEN-FORMATTED INPUT FILE" $item`
d=${c%:*}
num=$((b+1))
echo $num
csplit $item $num $d
mv xx01 $item.molden
done


#a='hello:world:of:tomorrow'

#echo "${a%:*}"
#echo "${a%%:*}"
#echo "${a#*:}"

