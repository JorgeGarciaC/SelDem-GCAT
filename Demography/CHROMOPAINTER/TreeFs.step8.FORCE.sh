#!/bin/bash
#SBATCH -J F_tree
#SBATCH -e log_%x_%j.err
#SBATCH -o log_%x_%j.out
#####SBATCH -t 0-10:00          
#SBATCH --mem-per-cpu=16G 
#SBATCH --array=1-3

i=$SLURM_ARRAY_TASK_ID

Fs_combine=Fs_combine
FFILE=forcefile.original.txt

#mkdir $Fs_combine

module load fineSTRUCTURE/2.1.0;
fs fs -X -Y -F $FFILE -m T -c 0.377082 $Fs_combine/Fs_combine_RES.chunkcounts.out $Fs_combine/Fs_mcmc_RES.chunkcounts.out.$i $Fs_combine/Fs_tree_RES_chunkcounts.$i
fs fs -X -Y -F $FFILE -m T -c 0.377082 $Fs_combine/Fs_combine_RES.chunklengths.out $Fs_combine/Fs_mcmc_RES.chunklengths.out.$i $Fs_combine/Fs_tree_RES_chunklengths.$i
