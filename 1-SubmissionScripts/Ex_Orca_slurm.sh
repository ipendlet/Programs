#!/bin/bash
#SBATCH --array=1-712
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --time=160:00:00
#SBATCH -e /export/zimmerman/ipendlet/4-CobaltChemistry/3_TestSet_PreliminaryOutFiles/wB97X-D_SolventCalc_ORCA/error.e
#SBATCH -o /export/zimmerman/ipendlet/4-CobaltChemistry/3_TestSet_PreliminaryOutFiles/wB97X-D_SolventCalc_ORCA/error.o
#SBATCH -p guest --job-name=H2Wb.qsh

time

. /etc/profile.d/slurm.sh

ID=`printf "%0*d\n" 4 ${SLURM_ARRAY_TASK_ID}`
cd $SLURM_SUBMIT_DIR
name=`ls q$ID.*.inp`

module unload Openmpi
shtcut="/export/applications"

export LD_LIBRARY_PATH=/export/apps/Intel/composer_xe_2013.4.183/compiler/lib/intel64:/export/zimmerman/khyungju/OpenMPI/2.0.2/lib:$LD_LIBRARY_PATH
export PATH=//export/zimmerman/khyungju/OpenMPI/2.0.2/bin:$PATH

cp $name $SLURM_LOCAL_SCRATCH
cd $SLURM_LOCAL_SCRATCH
/export/zimmerman/khyungju/orca_4_0_0_2_linux_x86-64/orca $name > $SLURM_SUBMIT_DIR/$name.out

wait
