#!/bin/bash
#SBATCH --job-name=bwa_index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o ../log_files/%x_%A.out
#SBATCH -e ../log_files/%x_%A.err

hg19=/home/FCAM/nperera/Tutorial/variant_detection_GATK/Illumina/Analysis_2/hg19/hg19.fa

module load bwa/0.7.17

d1=../raw_data

##################################################################
## Create BWA Index
##################################################################
echo "=========== hg19 BWA index  ================="
cd ../hg19/
bwa index ${hg19}

echo "=========== hg19 BWA index Done ================="
