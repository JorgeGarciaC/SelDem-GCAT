#!/bin/bash
#SBATCH -J Combine.DvsR
#SBATCH -e log_%x_%j.err
#SBATCH -o log_%x_%j.out
#SBATCH -t 0-08:00 #2-12:00     

#SBATCH --mem-per-cpu=8G 
###SBATCH -a 1-22%1

#i=$SLURM_ARRAY_TASK_ID

scripts=SCRIPTS
dir=ChPruns

Fs_combine=Fs_combine

module load fineSTRUCTURE/2.1.0

fs combine -v -f forcefile.original.txt -F FORCE_Fs/forcefile.txt  -d ${dir}/ -o $Fs_combine/Fs_combine_RES
