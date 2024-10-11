#!/bin/bash
#SBATCH -J aCLUES_hg38
#SBATCH -e log_%x_%j_%a.err
#SBATCH -o log_%x_%j_%a.out
#SBATCH -t 0-18:00                                                                                                                  
#SBATCH --mem-per-cpu=8G 
#SBATCH -n 1 # Number of cores (task per node)
#SBATCH -N 1 # All cores on one machine (nodes)

####################################################
# Input: anc, haps, sample, poplabel and mut files #
# Output: timeb file for CLUES                     #
####################################################

module load Python/3.6.6-foss-2018b
module load bcftools 

relate_path=~/scratch/relate1.1.x
clues_path=/homes/users/jgarciac/scratch/CLUES

chr=$1
rsid=$2
bp=$3
out=${rsid}
out_dir=$4
pop=$5
prevpos=$6
peak=$7
symbol=$8

i=$(echo $chr | cut -f 1 -d "," | sed 's/chr//')

##### Input

coal=/homes/users/jgarciac/scratch/RelateGRCh38/${pop}.popsize.coal
anc=/homes/users/jgarciac/scratch/RelateGRCh38/ANC.MUTS.${pop}/${pop}_upd_chr${i}.anc.gz

aDNA_samps=~/scratch/POP_gen_West_EUR/all.ancient.EUROPE.txt
aDNA_times=~/scratch/POP_gen_West_EUR/all.ancient.EUROPE.samplingtimes.txt
aDNA_vcf=~/scratch/POP_gen_West_EUR/1KG_lenient/${i}.filter.vcf.gz

anc_names=$(basename $anc '.anc.gz')
dir=$(dirname $anc)

# Sample Branch Length for the desired bp position
${relate_path}/scripts/SampleBranchLengths/SampleBranchLengths.sh \
    -i ${dir}/${anc_names} \
    -o ${out_dir}/${anc_names}_${out} \
    -m 1.25e-8 \
    --coal $coal \
    --format b \
    --num_samples 100 \
   --first_bp $bp \
    --last_bp $bp &&

# We extract the Genotype Posterior likelihood in a format needed for clues, in whic
# first is the time in generations, next the probablity of 0|0, then 0|1 and 1|1

is_flipped=$( zcat ${dir}/${anc_names}.mut.gz | grep ";$rsid;" | cut -f 8 -d ";" )

bcftools query -S $aDNA_samps -r ${i}:${prevpos} -f '[ %GP\n]' $aDNA_vcf | sed 's/,/ /g' | paste $aDNA_times - -d "" | sort -n > ${out_dir}/${rsid}.likelihoods.txt
echo 0 > ${out_dir}/timeBins.${rsid}.txt
bin_end=$(head -n 1 $aDNA_times)
head -n 1 $aDNA_times >> ${out_dir}/timeBins.${rsid}.txt

# next, if the snp is not covered in the ancient dataset, we will skip
numb_cols=$(awk '{print NF;exit}' ${out_dir}/${rsid}.likelihoods.txt)

if [[ $numb_cols -eq 4 ]]; 
    then  
    echo "Calculating Clues using ancient genotypes"
    echo "$is_flipped"
    if [[ $is_flipped -eq 1 ]]
        then 
        echo "Flipping ancient genotypes"
        awk '{print $1,$4,$3,$2}' ${out_dir}/${rsid}.likelihoods.txt > ${out_dir}/aux.${rsid}.lkl
        mv ${out_dir}/aux.${rsid}.lkl ${out_dir}/${rsid}.likelihoods.txt
    fi
        echo "Beggining inference"
        python ${clues_path}/inference.py \
            --times ${out_dir}/${anc_names}_${out} \
            --coal $coal \
            --out ${out_dir}/anc.${anc_names}_${out}_c \
            --ancientSamps ${out_dir}/${rsid}.likelihoods.txt \
            --timeBins ${out_dir}/timeBins.${rsid}.txt \
            --tCutoff 1000 | tee > ${out_dir}/anc.${rsid}.txt &&
    
        LR=$(grep "logLR" ${out_dir}/anc.${rsid}.txt | cut -f 2 -d "")
        s=$(grep "0-" ${out_dir}/anc.${rsid}.txt | cut -f 2 -d "")
        echo -e "$chr\t$rsid\t$pos\t$LR\t$s" >> ${out_dir}/anc.report.selcoeff.txt
        # We extract the logLR and the selection coefficient to both generate a report and make the plot
        python ${clues_path}/plot_custom_traj.py ${out_dir}/anc.${anc_names}_${out}_c ${out_dir}/anc.${anc_names}_${out}_plot \
            --ext png \
            --title "${rsid} - ${symbol}" \
            --subtitle "$LR ; s: $(echo $s | awk '{print $2}')"   \
            --subtitleR  "Chromosome ${i}, peak ${peak}"
        echo "Finished! Exiting"
    else
    echo "Variant no present in ancient samples. Exiting without infering"
    fi
