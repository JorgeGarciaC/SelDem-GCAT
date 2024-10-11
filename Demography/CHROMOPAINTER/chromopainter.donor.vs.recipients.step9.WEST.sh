#!/bin/bash
#SBATCH -J ChP.DvsR
#SBATCH -e log_%x_%j.err
#SBATCH -o log_%x_%j.out
#SBATCH -t 0-08:00 #2-12:00           
#SBATCH --mem-per-cpu=8G 
#SBATCH -a 1-22

i=$SLURM_ARRAY_TASK_ID

scripts=SCRIPTS
dir=07.03.2023
s_dir=WESTE.11.12.23

# Needs to add the CHR as chr${i} in order to be recognized by the program.
phasefile=${dir}/DEMO.chr${i}.phased.phase
recfile=${dir}/rec.files/chr${i}.b38.hapmap.recomrates
labels_clusters=${dir}/${s_dir}/idfile.globetrotter.txt
clust_file=${dir}/${s_dir}/popfile.1.txt
out_chunks=${dir}/${s_dir}/CHUNKS

Ne=136.0208
Mu=0.0007192314

# CHUNKS
sbatch -J ChP.DR.${i}.A -o log_WESTE_%x_%j_${i}.A.out -e log_EASTE_%x_%j_${i}.A.err -t 2-00:00 --wrap="$scripts/ChromoPainterv2 -g $phasefile -r $recfile -t $labels_clusters -f $clust_file 0 0 -o ${out_chunks}/$(basename $phasefile .phase).A -s 0 -n $Ne -M $Mu"

out_samples=${dir}/${s_dir}/SAMPLES
clust_file=${dir}/${s_dir}/popfile.2.txt

# SAMPLES
sbatch -J ChP.DR.${i}.A -o log_WESTE_%x_%j_${i}.A.out -e log_EASTE_%x_%j_${i}.A.err -t 0-16:00 --wrap="$scripts/ChromoPainterv2 -g $phasefile -r $recfile -t $labels_clusters -f $clust_file 0 0 -o ${out_samples}/$(basename $phasefile .phase).A -s 10 -n $Ne -M $Mu"

