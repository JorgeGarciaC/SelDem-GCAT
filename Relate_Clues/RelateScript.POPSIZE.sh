#!/bin/bash

#SBATCH -J GCAT_rel_popsize
#SBATCH -e log_%x_%j_%a.err
#SBATCH -o log_%x_%j_%a.out
#SBATCH -t UNLIMITED                                                                                                                                                                                                                                                                                 
#SBATCH --mem-per-cpu=20G 
#SBATCH -n 1 # Number of cores (task per node)
#SBATCH -N 1 # All cores on one machine (nodes)
#SBATCH --cpus-per-task=16
##SBATCH --array=22-22%5 

module load R/3.6.0-foss-2018b
module load GCC/12.3.0

relate_path=~/scratch/relate_v1.2.1_x86_64_dynamic

name=merged.bi.annotate.rename.phased.raw
pop=$1
dir=$2
poplabels=$3

# Previously we have to rename the anc mut to reassemble this nomenclature

${relate_path}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh \
              -i ${dir}/${name} \
              -m 1.25e-8 \
              --poplabels $poplabels \
              --pop_of_interest $pop \
              --threads 16 \
              --num_iter 5 \
              --first_chr 1 \
              --last_chr  22 \
              --noanc 1 \
              -o ${pop}.popsize
