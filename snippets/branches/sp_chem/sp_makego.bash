#!/bin/bash

echo "setting up go script "
file=sp_go.qsh

echo "#PBS -t " >> $file
echo "#PBS -l nodes=1:ppn=4 -l pmem=2000MB -l walltime=90:00:00" >> $file
echo "#PBS -e `pwd` -o `pwd`" >> $file
echo "#PBS -q zimmerman -N" "opt"${PWD##*/}".qsh"  >> $file
echo " " >> $file
echo 'ID=`printf "%0*d\n" 3 ${PBS_ARRAYID}`' >> $file
echo " " >> $file
echo "source /export/zimmerman/paulzim/qchem/qchemjan42013c/paul.set.local" >> $file
echo " " >> $file
echo 'cd $PBS_O_WORKDIR' >> $file
echo 'name=`ls q$ID*.inp`' >> $file
echo 'qchem -np 4 $name $name.out' >> $file
echo " " >> $file
echo 'rm $QCSCRATCH' >>	$file

