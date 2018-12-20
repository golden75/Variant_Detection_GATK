# Variant Detection using GATK 

### Introduction to Variant Detection using Whole Exome Sequencing 
Genome-wide sequencing methods have provided a method of deep understanding in the  sequence variations between genotype and phenotype. The 1000 genome project (http://www.1000genomes.org) have indeed made inroads into the study of population genetics, including the investigation of causal variants of genes for various human syndromes. Next Generation Sequencing (NGS) technology is evolving at a rapid rate and new sequencing platforms are been developed. Whole Genome Sequencing (WGS) have been used comprehensively to detecting genomic variations such as single nucleotide variants (SNVs), Copy number variants (CNVs), insertions and deletions (InDels) and chromosomal rearrangements.  However at low cost, whole exome sequencing (WES) platforms, have performed well in high sequencing coverage and readily interpreting protein coding exons associated compared to WGS platforms. This technique has also led to more samples to be analyzed and can generate targeted DNA sequences and identify substantially more genetic variations.

The exome is the part of the genome composed of exons, the sequences which, when transcribed, remain within the mature RNA after introns are removed by RNA splicing and contribute to the final protein product encoded by that gene. It consists of all DNA that is transcribed into mature RNA in cells of any type, as distinct from the transcriptome, which is the RNA that has been transcribed only in a specific cell population. The exome of the human genome consists of roughly 180,000 exons constituting about 1% of the total genome.


### Sample Data Download 
For this pipe line we are using **NA12878** test sample where it contain SRX655430: Illumina random exon sequencing of genomic DNA paired-end library 'Pond-314378' containing sample 'NA12878'[[link](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1517906)]

<pre>
INPUT_FILES=(SRR1517848 SRR1517878 SRR1517884 SRR1517906 SRR1517991 SRR1518011 SRR1518158 SRR1518253)

d1="raw_data"

module load sratoolkit/2.8.2

if [ ! -d ../${d1} ]; then
        mkdir -p ../${d1}
fi

cd ../${d1}

fastq-dump --split-files ${INPUT_FILE_NAME}
</pre>
The full script for slurm shedular can be found in the scripts folder by the name data_download.sh. <a href="/scripts/data_download.sh">here</a>

In this step we will download the paired-end interleveled fastq files from the NCBI SRA database. While downloading we will be splitting the fastq files into two fastq files; forward and reverse (R1 and R2)strand by using `--split-files` command, which gives us:

<pre>
raw_data/
├── SRR1517848_1.fastq
├── SRR1517848_2.fastq
├── SRR1517878_1.fastq
├── SRR1517878_2.fastq
├── SRR1517884_1.fastq
├── SRR1517884_2.fastq
├── SRR1517906_1.fastq
├── SRR1517906_2.fastq
├── SRR1517991_1.fastq
├── SRR1517991_2.fastq
├── SRR1518011_1.fastq
├── SRR1518011_2.fastq
├── SRR1518158_1.fastq
├── SRR1518158_2.fastq
├── SRR1518253_1.fastq
└── SRR1518253_2.fastq
</pre>


<pre style="color: silver; background: black;">
d1="raw_data"

cd ../${d1}

module load sratoolkit/2.8.2

list="SRR796868 SRR796869 SRR796870 SRR796871 SRR796872 
SRR796873 SRR796874 SRR796875 SRR796876 SRR796877 
SRR796878 SRR796879 SRR796880 SRR796881 SRX265476"

for file in $list; do
        echo $file
        fastq-dump --split-files $file
done
</pre>
