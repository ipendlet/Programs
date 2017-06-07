#! /bin/bash

name=`ls *.out`

#GASE=$(grep -i "INTERNAL ENERGY IN GAS" q001.3a-227-244.xyz_converted.inp.out)
for item in $name
do
SOLVE=$(grep -i 'TOTAL FREE ENERGY IN SOLVENT' $item)
#GASEonly=$(sed -n '/INTERNAL ENERGY IN GAS = /,/A.U./p' $name.out)
#echo $SOLVE
echo $item $SOLVE "\n" >> calculatedsolventenergies.txt
done
