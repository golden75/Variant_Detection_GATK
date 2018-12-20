#!/bin/bash
#SBATCH --job-name=data_download
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[1-8]%8
##SBATCH --mail-type=ALL
##SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o ../log_files/%x_%A.out
#SBATCH -e ../log_files/%x_%A.err

d1="raw_data"

INPUT_FILES=(SRR1517848 SRR1517878 SRR1517884 SRR1517906 SRR1517991 SRR1518011 SRR1518158 SRR1518253)
INPUT_FILE_NAME="${INPUT_FILES[$SLURM_ARRAY_TASK_ID - 1]}"

echo "host name : " `hostname`
echo "input file name : " $INPUT_FILE_NAME 
echo "SLURM_ARRAY_TASK_ID : " $SLURM_ARRAY_TASK_ID

##################################################################
## Download data
##################################################################
echo "=========== Download of reads :: ${INPUT_FILE_NAME}  ================="
module load sratoolkit/2.8.2

if [ ! -d ../${d1} ]; then
        mkdir -p ../${d1}
fi

cd ../${d1}

fastq-dump --split-files ${INPUT_FILE_NAME}

echo "Download complete ${INPUT_FILE_NAME}  `date`"
