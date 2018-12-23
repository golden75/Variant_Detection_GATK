#!/bin/bash
#SBATCH --job-name=remove_duplicates
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
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
d3=noduplicates


INPUT_FILES=(SRR1517848 SRR1517878 SRR1517884 SRR1517906 SRR1517991 SRR1518011 SRR1518158 SRR1518253)
INPUT_FILE_NAME="${INPUT_FILES[$SLURM_ARRAY_TASK_ID - 1]}"

echo "host name : " `hostname`
echo "input file name : " $INPUT_FILE_NAME 
echo "SLURM_ARRAY_TASK_ID : " $SLURM_ARRAY_TASK_ID

##################################################################
## Mark Duplicates    
##################################################################
echo "=== Mark Duplicates starts `date` ==="
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

if [ ! -d ../${d3} ]; then
       mkdir -p ../${d3}
fi

cd ../${d3}/

java -jar $PICARD MarkDuplicates \
        INPUT=../${d2}/${INPUT_FILE_NAME}_filtered_sort.bam \
        OUTPUT=${INPUT_FILE_NAME}_nodup.bam \
        REMOVE_DUPLICATES=Ture \
        METRICS_FILE=${INPUT_FILE_NAME}_metrics.txt \
        CREATE_INDEX=True

echo "=== Mark Duplicates : ${INPUT_FILE_NAME} `date`==="

