#!/bin/bash
#SBATCH -J ChP
#SBATCH -e log_%x_%j.err
#SBATCH -o log_%x_%j.out
#SBATCH -t 0-08:00 #2-12:00           
#SBATCH --mem-per-cpu=8G 
#SBATCH -a 1-22%1

i=$SLURM_ARRAY_TASK_ID

scripts=SCRIPTS
dir=07.03.2023
# Needs to add the CHR as chr${i} in order to be recognized by the program.
phasefile=${dir}/DEMO.chr${i}.phased.phase
recfile=${dir}/rec.files/chr${i}.b38.hapmap.recomrates
labels=$dir/idfile.txt
out=$dir/ChPruns
#mkdir $out


# We use the parameters estimate in step3
SWrate=136.0208 # 
Mu=0.0007192314 #
m=1
n=10

#Hacemos cuenta de 1/10 de las muestras que tenemos

for count in {1..117}; do
sbatch -J ChP.Run.${i}.${count} \
    -o log_%x_%j_${count}.out \
    -e log_%x_%j_${count}.err \
    -t 0-06:00 \
    --wrap="$scripts/ChromoPainterv2 -g $phasefile -r $recfile -t $labels -o ${out}/$(basename $phasefile .phase).${count} -a $m $n -i 0 -n $SWrate -M $Mu"
m=$(( $m + 10 ))
n=$(( $n + 10 )); done # Aquí añadimos diez a cada número para que el programa nos cuente de 10 en 10 

# En este caso el ultimo job va con 11 muestras pues sino, nos sobraría una muestra

sbatch -J ChP.Run.${i}.118 -o log_%x_%j_118.out -e log_%x_%j_118.err -t 0-06:00 --wrap="$scripts/ChromoPainterv2 -g $phasefile -r $recfile -t $labels -o ${out}/$(basename $phasefile .phase).118 -a 1171 1181 -i 0 -n $SWrate -M $Mu"
