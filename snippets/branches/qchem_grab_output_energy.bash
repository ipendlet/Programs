#! /bin/bash

list=`ls *.inp.out`
for item in $list
do
energy=`grep 'Final energy is' $item`
# >> optimized_energies_output_fromqchem.txt
echo $item $energy >> optimized_energies_output_fromqchem.txt
done


