#!/bin/bash
#SBATCH --job-name=singletons
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[1-8]%8
##SBATCH --mail-type=ALL
##SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH -o ../log_files/%x_%A_%a.out
#SBATCH -e ../log_files/%x_%A_%a.err

hg19=/home/FCAM/nperera/Tutorial/variant_detection_GATK/Illumina/hg19/hg19.fa
R1="_1.fastq"
R2="_2.fastq"

d1="raw_data"
d2=align


INPUT_FILES=(SRR1517848 SRR1517878 SRR1517884 SRR1517906 SRR1517991 SRR1518011 SRR1518158 SRR1518253)
INPUT_FILE_NAME="${INPUT_FILES[$SLURM_ARRAY_TASK_ID - 1]}"

echo "host name : " `hostname`
echo "input file name : " $INPUT_FILE_NAME 
echo "SLURM_ARRAY_TASK_ID : " $SLURM_ARRAY_TASK_ID


if [ ! -d ../${d2} ]; then
        mkdir -p ../${d2}
fi

cd ../${d2}

##################################################################
## SAM to BAM
##################################################################
echo "=== SAM to BAM `date` ===
module load samtools/1.7
samtools view -@ 8 -bS ${INPUT_FILE_NAME}.sam > ${INPUT_FILE_NAME}.bam
echo "bam file created for ${INPUT_FILE_NAME} `date`"

##################################################################
## Remove sigletons
##################################################################
echo "=== Remove singletons `date` ==="
samtools view -@ 8 -F 0x04 -b ${INPUT_FILE_NAME}.bam > ${INPUT_FILE_NAME}_filtered.bam
echo "sigletons removed for ${INPUT_FILE_NAME} `date`"
module unload samtools/1.7


echo "=== SAM to BAM and Singletons Reomved: ${INPUT_FILE_NAME}  `date`==="
