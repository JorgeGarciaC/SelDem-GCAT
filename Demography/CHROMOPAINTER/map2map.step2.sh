#!/bin/bash
#SBATCH -J map2map
#SBATCH -e log_%x_%j.err
#SBATCH -o log_%x_%j.out
#SBATCH -t 0-24:00          
#SBATCH --mem-per-cpu=8G 
#SBATCH -a 1-22%10

i=$SLURM_ARRAY_TASK_ID

scripts=~/scratch/CHROMOPAINTER/SCRIPTS

# Needs to add the CHR as chr${i} in order to be recognized by the program.
genetic_map_ChP=~/scratch/CHROMOPAINTER/G.MAP/chr${i}.b38.hapmap.txt
phasefile=~/scratch/CHROMOPAINTER/07.03.2023/DEMO.chr${i}.phased
recfile=$(basename $genetic_map_ChP .txt)
out=~/scratch/CHROMOPAINTER/07.03.2023/recom.rates

module load Perl/5.18.2;
perl $scripts/convertrecfile.pl -M hapmap $phasefile $genetic_map_ChP ${out}/${recfile}.recomrates # - M hapmap when map 4 columns y la 4 es cM o plain when 2 y la 2 es rate
