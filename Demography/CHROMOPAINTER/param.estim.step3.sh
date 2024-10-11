#!/bin/bash
#SBATCH -J param.estim
#SBATCH -e log_%x_%j.err
#SBATCH -o log_%x_%j.out
#SBATCH -t 0-04:00         
#SBATCH --mem-per-cpu=16G 
##SBATCH -a 20-20%15

i=$SLURM_ARRAY_TASK_ID

for i in 1 4 17 20;
do

scripts=SCRIPTS

#RUN FOR CHR 1,4, 17, 20

dir=07.03.2023

# Needs to add the CHR as chr${i} in order to be recognized by the program.
phasefile=${dir}/DEMO.chr${i}.phased.phase
recfile=${dir}/rec.files/chr${i}.b38.hapmap.recomrates
labels=$dir/idfile.txt # 1181
out=$dir/param.estim

mkdir $out

m=1
n=10

for count in {1..117}; do # num haplotipos / 10 - el úlitmo job si no es múltiplo de diez
sbatch \
    -J estim.${count} \
    -o log_%x_%j_${count}.out \
    -e log_%x_%j_${count}.err \
    -t 2-00:00 \
    --wrap="$scripts/ChromoPainterv2 -g $phasefile -r $recfile -t $labels -o ${out}/$(basename $phasefile .phase).n.${count} -a $m $n -i 10 -in -iM"

m=$(( $m + 10 ))
n=$(( $n + 10 )); done


sbatch \
   -J estim.118 \
   -o log_%x_%j_118.out \
   -e log_%x_%j_118.err \
   -t 2-00:00 \
   --wrap="$scripts/ChromoPainterv2 -g $phasefile -r $recfile -t $labels -o ${out}/$(basename $phasefile .phase).n.118 -a 1171 1181 -i 10 -in -iM"


#only needed if last job can't be separate, i.e. if is not 10 multiple

done
