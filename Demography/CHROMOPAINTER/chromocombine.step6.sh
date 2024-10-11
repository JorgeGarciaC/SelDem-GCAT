#!/bin/bash
#SBATCH -J combine
#SBATCH -e log_%x_%j.err
#SBATCH -o log_%x_%j.out
#SBATCH -t 0-02:00             
#SBATCH --mem-per-cpu=32G 

out=07.03.2023/ChPruns

Fs_combine=07.03.2023/Fs_combine

mkdir $Fs_combine

module load fineSTRUCTURE/2.1.0
fs combine -v -i 07.03.2023/idfile.txt -d ${out}/ -o $Fs_combine/Fs_combine_RES
