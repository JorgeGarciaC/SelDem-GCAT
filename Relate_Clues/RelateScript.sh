#!/bin/bash

#SBATCH -J GCAT_rel_par
#SBATCH -e log_%x_%j_%a.err
#SBATCH -o log_%x_%j_%a.out
#SBATCH -t UNLIMITED                                                                                                                                                                                                                                                                                  
#SBATCH --mem-per-cpu=20G 
#SBATCH -n 1 # Number of cores (task per node)
#SBATCH -N 1 # All cores on one machine (nodes)
#SBATCH --cpus-per-task=8
#SBATCH --array=3-11%5 

module load GCC/12.3.0

i=$SLURM_ARRAY_TASK_ID
vcf=/homes/users/jgarciac/scratch/GCAT.1kG38/NORM/PHASED/chr${i}.merged.bi.annotate.rename.phased
haps_dir=~/scratch/RelateGRCh38/HAPS.SAMPLE.3pops
relate_path=~/scratch/relate_v1.2.1_x86_64_dynamic
name=$(basename $vcf .vcf)
ancestral=~/projects/shared_data/1000G_2504_high_coverage/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor_${i}.fa
mask=~/projects/shared_data/1000G_2504_high_coverage/MASKS.38/20160622.chr${i}.mask.fasta.gz
poplabels=~/scratch/RelateGRCh38/relate.GRCh38.poplabels
map=~/projects/shared_data/1000G_2504_high_coverage/REF/b38_maps/chr${i}.rename.gmap

coal_leo=~/projects/shared_data/relate/ALL_RESULTS/1000GP_Phase3_mask_prene.coal
dist=${haps_dir}/chr${i}.merged.bi.annotate.rename.phased.prepared.dist.gz
anno=${haps_dir}/chr${i}.merged.bi.annotate.rename.phased.prepared.annot

outdir=ANC.MUTS.RAW.3pops 

echo "CREATING HAPS-SAMPLE FILES" &&
${relate_path}/bin/RelateFileFormats \
                 --mode ConvertFromVcf \
                 --haps ${haps_dir}/${name}.haps \
                 --sample ${haps_dir}/${name}.sample \
                 -i $vcf &&

echo "Finishing HAPS-SAMPLE FILES. PREPARING INPUT" &&
${relate_path}/scripts/PrepareInputFiles/PrepareInputFiles.sh \
                --haps ${haps_dir}/${name}.haps \
                --sample ${haps_dir}/${name}.sample \
                --ancestor ${ancestral} \
                --mask ${mask} \
                --poplabels ${poplabels} \
                --remove_ids ~/scratch/RelateGRCh38/remove.samples.relate.txt \
                -o ${haps_dir}/${name}.prepared &&

cd $outdir # Neccesary as relate works on current directory

# same N as original 
echo "BEGINNING INFERRING THE ARG" &&
${relate_path}/scripts/RelateParallel/RelateParallel.sh \
      -m 1.25e-8 \
      -N 30000 \
      --haps ${haps_dir}/${name}.prepared.haps.gz  \
      --sample ${haps_dir}/${name}.prepared.sample.gz \
      --map $map \
      --max_memory 19 \
      -o ${name}.raw \
      --coal $coal_leo \
      --dist $dist \
      --anot $anno \
      --threads 8 &&

echo "FINISH INFERRING THE ARG. EXITING"
