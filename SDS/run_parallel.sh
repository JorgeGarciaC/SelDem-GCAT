#!/bin/bash
#SBATCH -J parallel_SDS
#SBATCH -e log_%x%j.err
#SBATCH -o out_%x%j.out
#SBATCH -t 0-24:00                                                                                                                                                                                                                                                                                
#SBATCH --mem-per-cpu=2G 

module load R

i=$SLURM_ARRAY_TASK_ID
chr=$1

dir=/homes/users/jgarciac/scratch/split/chr${chr}
temp=~/scratch/temp_SDS/chr${chr}
out=~/scratch/SDSs/chr${chr}
b_file=~/SDS/boundaries/centromeres_only/b_${chr}.txt # Needed to take out boundaries
mkdir $out

time Rscript ~/SDS/compute_SDS.R \
                ${temp}/s_file.txt \
                ${dir}/F${i} \
                ~/scratch/o_file.txt \
                $b_file \
                ~/SDS/gcat_thigh.gamma_shapes \
                1e-6 > ${out}/${i}

