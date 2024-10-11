#!/bin/bash
#SBATCH -J input.chromo
#SBATCH -e log_%x%j.err
#SBATCH -o log_%x%j.out
#SBATCH -t 0-24:00     
#SBATCH --mem-per-cpu=12G 
#SBATCH -a 1-22%5

i=$SLURM_ARRAY_TASK_ID
out=07.03.2023

haps_dir=/homes/users/jgarciac/scratch/CHROMOPAINTER/PHASED.DEMO.DATA/
relate_path=~/scratch/relate

${relate_path}/bin/RelateFileFormats \
                 --mode RemoveSamples \
                 --haps $haps_dir/DEMO.chr$i.phased.haps \
                 --sample $haps_dir/DEMO.chr$i.phased.sample \
                 -i id.remove.01.09.2023.txt \
                 -o $out/HAPS_SAMPLE/DEMO.$i

scripts=SCRIPTS
genetic_map=~/projects/shared_data/1000G_2504_high_coverage/REF/b38_maps/chr${i}.b38.gmap.gz
outputSHAPEIT=$out/HAPS_SAMPLE/DEMO.$i.haps

module load Perl/5.18.2;
perl $scripts/impute2chromopainter.pl -r $genetic_map $outputSHAPEIT ${out}/DEMO.chr$i.phased
