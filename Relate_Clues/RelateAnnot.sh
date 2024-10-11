#!/bin/bash

#SBATCH -J GCAT_ann
#SBATCH -e log_%x_%j_%a.err
#SBATCH -o log_%x_%j_%a.out
#SBATCH -t UNLIMITED                                                                                                                                                                                                                                                                                                                                                                                           
#SBATCH --mem-per-cpu=20G 
#SBATCH -n 1 # Number of cores (task per node)
#SBATCH -N 1 # All cores on one machine (nodes)
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22%5 

module load GCC/12.3.0

i=$SLURM_ARRAY_TASK_ID
vcf=/homes/users/jgarciac/scratch/GCAT.1kG38/NORM/PHASED/chr${i}.merged.bi.annotate.rename.phased
haps_dir=~/scratch/RelateGRCh38/HAPS.SAMPLE.3pops
relate_path=~/scratch/relate_v1.2.1_x86_64_dynamic
name=$(basename $vcf .vcf)
ancestral=~/projects/shared_data/1000G_2504_high_coverage/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${i}.fa
mask=~/projects/shared_data/1000G_2504_high_coverage/MASKS.38/20160622.chr${i}.mask.fasta.gz
poplabels=~/scratch/RelateGRCh38/HAPS.SAMPLE.3pops/chr1.merged.bi.annotate.rename.phased.prepared.poplabels
map=~/projects/shared_data/1000G_2504_high_coverage/REF/b38_maps/chr${i}.rename.gmap

coal_leo=~/projects/shared_data/relate/ALL_RESULTS/1000GP_Phase3_mask_prene.coal
dist=${haps_dir}/chr${i}.merged.bi.annotate.rename.phased.prepared.dist.gz
anno=${haps_dir}/chr${i}.merged.bi.annotate.rename.phased.prepared.annot

outdir=ANC.MUTS.RAW.3pops 

cd $outdir

${relate_path}/bin/RelateFileFormats \
                 --mode GenerateSNPAnnotations \
                 --haps ${haps_dir}/${name}.prepared.haps.gz \
                 --sample ${haps_dir}/${name}.prepared.sample.gz \
                 --poplabels ${poplabels} \
                  --ancestor ${ancestral} \
                  --mut merged.bi.annotate.rename.phased.raw_chr${i}.mut.gz \
                 -o merged.bi.annotate.rename.phased.raw_chr${i}
