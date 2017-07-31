#!/bin/bash
file=runorca_pbs.qsh

echo "#PBS -t " >> $file
echo "#PBS -l nodes=1:ppn=1 -l walltime=30:00:00 " >> $file
echo "#PBS -e `pwd` -o `pwd`" >> $file
echo "#PBS -q zimmerman -N orca_run.qsh" >> $file
echo " " >> $file
echo 'ID=`printf "%0*d\n" 4 ${PBS_ARRAYID}`' >> $file
echo " " >> $file
echo 'cd $PBS_O_WORKDIR' >> $file
echo 'name=`ls q$ID.*.inp`' >> $file
echo "module unload Openmpi" >> $file
echo "shtcut='/export/applications'" >> $file
echo 'export LD_LIBRARY_PATH=/export/apps/Intel/composer_xe_2013.4.183/compiler/lib/intel64:/export/zimmerman/khyungju/OpenMPI/2.0.2/lib:$LD_LIBRARY_PATH' >> $file
echo 'export PATH=//export/zimmerman/khyungju/OpenMPI/2.0.2/bin:$PATH' >> $file
echo " " >> $file
echo 'cp $name $PBSTMPDIR' >> $file
echo 'cd $PBSTMPDIR' >> $file
echo '/export/zimmerman/khyungju/orca_4_0_0_2_linux_x86-64/orca $name > $PBS_O_WORKDIR/$name.out' >> $file
echo " " >> $file
echo "wait" >> $file
echo "time" >> $file


