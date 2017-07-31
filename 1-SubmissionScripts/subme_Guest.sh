#!/bin/bash
#SBATCH --array=25-26
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=96:00:00
#SBATCH -e /export/zimmerman/ipendlet/4-CobaltChemistry/5_Guo-pKaAssessement/4-Screens-Frequency/error.e
#SBATCH -o /export/zimmerman/ipendlet/4-CobaltChemistry/5_Guo-pKaAssessement/4-Screens-Frequency/error.o
#SBATCH -p guest --job-name=4-Screens-Frequency.qsh
 
cd /export/zimmerman/ipendlet/4-CobaltChemistry/5_Guo-pKaAssessement/4-Screens-Frequency
 
item=$SLURM_ARRAY_TASK_ID
ID=`printf "%0*d\n" 3 ${item}`
 
echo "running $ID"
hostname
source /export/zimmerman/paulzim/qchem/qchemjan42013cg/paul.set.local0
 
export OMP_NUM_THREADS=1
name=`ls q$ID*.inp`
qchem -np 1 $name $name.out
 
rm $QCSCRATCH
