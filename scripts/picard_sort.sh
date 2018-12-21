#!/bin/bash
#SBATCH --job-name=picard_sort
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

hg19=/home/FCAM/nperera/Tutorial/variant_detection_GATK/Illumina/Analysis_2/hg19/hg19.fa
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
## Sort 
##################################################################
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

java -jar $PICARD SortSam \
        INPUT=${INPUT_FILE_NAME}_filtered.bam \
        OUTPUT=${INPUT_FILE_NAME}_filtered_sort.bam \
        SORT_ORDER=coordinate \
        CREATE_INDEX=True

echo "Sorted : " ${INPUT_FILE_NAME} `date`

module unload picard/2.9.2

echo "=== BAM file sort : ${INPUT_FILE_NAME}  `date`==="
