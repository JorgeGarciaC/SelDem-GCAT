#!/bin/bash
#SBATCH -J prepare_SDS
#SBATCH -e log_%x%j.err
#SBATCH -o log_%x%j.out
#SBATCH -t 0-04:00
#SBATCH --mem-per-cpu=6G
#SBATCH -a 1-22%5

##########################################
#                                        #
# Jorge García IBE/UPF                   #
# Input: A VCF.                         #
# Output: s_files and t_files for SDS   #
#                                        #
##########################################

module load bcftools 
module load VCFtools/0.1.16-foss-2018b-Perl-5.28.0
module load Python/3.8.6-GCCcore-10.2.0
module load PLINK

i=$SLURM_ARRAY_TASK_ID

################
# Prepare s_file
out=~/scratch/temp_SDS/chr${i}
mkdir $out
vcf=~/scratch/GCAT_NOPHASED_ANN/chr${i}_variants_GCAT_panel_genotypes.ibs.no_miss.switch.annotate.vcf.gz
name=$(basename $vcf '.vcf.gz')
tabix $vcf

echo -e "\nParsing $name\n"
time bcftools view $vcf -m 2 -M 2  | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' |  vcftools --vcf - --singletons --out ${out}/singletons

# delete doubletons and retain true singletons.
time awk 'BEGIN {OFS="\t"} {if ($3=="S") {print $2,$5}}' ${out}/singletons.singletons > ${out}/sing_pos.txt

# Join singletons to each sample with a python dictionary. Returns a tab archive 
echo -e "\nPaso awk terminado\n"
time ~/SDS/Singleton_Mat.py ${out}/sing_pos.txt ${out}/output.csv
echo -e "\nPaso matriz singletons completado\n" 

# Delete headers
time sort ${out}/output.csv | cut -d , -f 2- | sed 's/,/\t/g' > ${out}/s_file.txt
echo -e "\ns_file completado\n"

################################
# Prepare t_file

echo -e "\nTransform vcf $name to 012 genotype matrix\n"
# Here we pre proccess test snps

plink2 --vcf $vcf  \
    --maf 0.05 --max-maf 0.95 \
    --hwe 1e-50 \
    --geno 0.00001 \
    --double-id --allow-extra-chr --make-bed \
    --export vcf --out ${out}/${name} # take out columns count y alt y cm

python3.8 ~/SDS/MakeT_file.py ${out}/${name}.vcf ${out}/t_file.txt
variants_f=$(wc -l ${out}/t_file.txt)

dir=~/scratch/split/chr${i}
mkdir $dir 

# We break the t_file to run SDS in parallel.

awk -v dir=$dir 'NR%5000==1{x="F"++i;}{print > dir"/"x}' ${out}/t_file.txt
echo -e "\nFinal variants = $variants_f \n"

###########
# Launch SDS

sbatch --array=1-${n_files}%5 ~/run_parallel.sh $i
