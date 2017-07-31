#!/bin/bash

echo "setting up go script"
file=qchem_run.sh


echo "#!/bin/bash" >> $file
echo "#SBATCH --array=" >> $file
echo "#SBATCH --nodes=1" >> $file
echo "#SBATCH --ntasks=1" >> $file
echo "#SBATCH --time=3:00:00" >> $file
echo "#SBATCH -e `pwd`" >> $file
echo "#SBATCH -o `pwd`" >> $file
echo "#SBATCH -p guest --job-name=Qchem_slurm.qsh" >> $file
echo " " >> $file
echo "cd `pwd`" >> $file
echo " " >> $file 
echo 'item=$SLURM_ARRAY_TASK_ID' >> $file
echo 'ID=`printf "%0*d\n" 3 ${item}`' >> $file
echo " " >> $file
echo 'echo "running $ID"' >> $file
echo "hostname" >> $file
echo 'source /export/zimmerman/paulzim/qchem/qchemjan42013cg/paul.set.local0' >> $file
echo " " >> $file
echo "export OMP_NUM_THREADS=2" >> $file
echo 'name=`ls q$ID*.inp`' >> $file
echo 'qchem -np 1 $name $name.out' >> $file
echo " " >> $file
echo 'rm $QCSCRATCH' >> $file
