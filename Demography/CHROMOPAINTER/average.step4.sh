#!/bin/bash
#SBATCH -J param.estim
#SBATCH -e log_%x_%j.err
#SBATCH -o log_%x_%j.out
#SBATCH -t 0-02:00      
#SBATCH --mem-per-cpu=3G 

scripts=SCRIPTS

for i in 1 4 17 20; do cat param.estim/DEMO.chr${i}.phased.n.*.EMprobs.out | awk 'BEGIN {n=0;m=0;count=0;FNR=1} {if (FNR%17==0) {n+=$4;m+=$5;count+=1;FNR=0}} END {print n/count,m/count,count} ' > average.chr${i}; done

### Check that you have made the division correctly

# Next, cat average.chr and add the snp numbers. 
# Make the weight pound as €(n or m * snp number / snp number)
