#!/bin/bash

#SBATCH -J GCAT_reestimate
#SBATCH -e log_%x_%j_%a.err
#SBATCH -o log_%x_%j_%a.out
#SBATCH -t UNLIMITED                                                                                                                                                                                                                                                                                 
#SBATCH --mem-per-cpu=20G 
#SBATCH -n 1 # Number of cores (task per node)
#SBATCH -N 1 # All cores on one machine (nodes)
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22%5 

#module load GCC/12.3.0

i=$SLURM_ARRAY_TASK_ID
#relate_path=~/scratch/relate_v1.2.1_x86_64_dynamic bug with this version on cluster
relate_path=~/scratch/relate1.1.x

name=merged.bi.annotate.rename.phased.raw_chr
popint=$1
indir=$2
poplabels=$3
haps_dir=$4
outdir=$5

dist=${haps_dir}/chr${i}.merged.bi.annotate.rename.phased.prepared.dist.gz

echo "SUBSET TREES" &&
${relate_path}/bin/RelateExtract \
                 --mode SubTreesForSubpopulation \
                 --anc ${indir}/${name}${i}.anc.gz \
                 --mut ${indir}/${name}${i}.mut.gz \
                 --poplabels ${poplabels} \
                 --pop_of_interest $popint \
                    -o ${outdir}/${popint}_chr${i} &&

echo "Reestimating Branch lengths" 
${relate_path}/scripts/SampleBranchLengths/ReEstimateBranchLengths.sh \
                 -i ${outdir}/${popint}_chr${i} \
                 -m 1.25e-8 \
                 --coal ${popint}.popsize.coal \
                 --threads 8 \
                 --dist $dist \
                 -o ${outdir}/${popint}_upd_chr${i} &&
echo "ENDING"
