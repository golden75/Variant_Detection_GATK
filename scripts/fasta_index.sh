#!/bin/bash
#SBATCH --job-name=fasta_index
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

##################################################################
## CreateFASTA Index
##################################################################
echo "=========== hg19 fasta file index  ================="
module load samtools/1.7
samtools faidx ${hg19}
echo "=========== fasta file index Done  ================="
