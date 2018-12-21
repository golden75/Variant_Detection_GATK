#!/bin/bash
#SBATCH --job-name=dictionary
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=neranjan.perera@uconn.edu
#SBATCH -o ../log_files/%x_%A.out
#SBATCH -e ../log_files/%x_%A.err

hg19=/home/FCAM/nperera/Tutorial/variant_detection_GATK/Illumina/Analysis_2/hg19/hg19.fa
ref=hg19


##################################################################
## CreateSequenceDictionary
##################################################################
echo "start: `date`"
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

if [ ! -d ../"$ref" ]; then
       mkdir -p ../"$ref"
fi

cd ../"$ref"/

java -jar $PICARD CreateSequenceDictionary \
        REFERENCE=hg19.fasta \
        OUTPUT=hg19.dict \
        CREATE_INDEX=True

echo "=========== hg19 sequence dictionary done  ================="
        
