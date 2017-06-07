#!/bin/bash

list="q003*.out"

for item in $list
do 
top=`grep -ni "Localized Orbital Bonding Analysis" $item | cut -d':' -s -f1`
bot=`grep -ni "Oxidation" $item | cut -d':' -s -f1`
csplit $item $top $bot
mv xx01 $item.orbitals.txt
done
