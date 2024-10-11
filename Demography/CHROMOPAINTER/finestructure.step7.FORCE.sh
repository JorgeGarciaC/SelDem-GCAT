#!/bin/bash
#SBATCH -J FMCMC
#SBATCH -e log_%x_%j.err
#SBATCH -o log_%x_%j.out
#SBATCH -t 07-12:00          
#SBATCH --mem-per-cpu=64G 
#SBATCH --array=2-3

i=$SLURM_ARRAY_TASK_ID

Fs_combine=Fs_combine
outdir=FORCE_Fs
FFILE=forcefile.original.txt

module load fineSTRUCTURE/2.1.0;
fs fs -m oMCMC -x 1000000 -y 2000000 -z 10000 -X -Y -F $FFILE -c 0.377082 $Fs_combine/Fs_combine_RES.chunkcounts.out $Fs_combine/Fs_mcmc_RES.chunkcounts.out.$i
fs fs -m oMCMC -x 1000000 -y 2000000 -z 10000 -X -Y -F $FFILE -c 0.377082 $Fs_combine/Fs_combine_RES.chunklengths.out $Fs_combine/Fs_mcmc_RES.chunklengths.out.$i
