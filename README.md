# Variant Detection using GATK 

This repository is a usable, publicly available tutorial for analyzing differential expression data and creating topological gene networks. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use <a href="https://bioinformatics.uconn.edu/unix-basics">this</a> handy guide for the operating system commands.  In this guide, you will be working with common bio Informatic file formats, such as <a href="https://en.wikipedia.org/wiki/FASTA_format">FASTA</a>, <a href="https://en.wikipedia.org/wiki/FASTQ_format">FASTQ</a>, <a href="https://en.wikipedia.org/wiki/SAM_(file_format)">SAM/BAM</a>, and <a href="https://en.wikipedia.org/wiki/General_feature_format">GFF3/GTF</a>. You can learn even more about each file format <a href="https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/">here</a>. If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one <a href="https://bioinformatics.uconn.edu/contact-us/">here</a>.

<div id="toc_container">
<p class="toc_title">Contents</p>
<ul class="toc_list">
    <li><a href="#Header_1"> 1. Introduction to Variant Detection using Whole Exome Sequencing</>
    <li><a href="#Header_2"> 2. Sample Data Download </>
    <li><a href="#Header_3"> 3. Preparing the Reference Sequence</>
    <li><a href="#Header_4"> 4. Aligning of Reads</>
    <li><a href="#Header_5"> 5. SAM to BAM Conversion and Remove Singletons</>
    <li><a href="#Header_6"> 6. Sort BAM files using PICARD</>
    <li><a href="#Header_7"> 7. Remove PCR Duplicates using PICARD</>
    <li><a href="#Header_8"> 8. Add Read Group Information</>
    <li><a href="#Header_9"> 9. Reorder BAM file</>
    <li><a href="#Header_10"> 10. Variant Calling</>
</ul>
</div>

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
The full script for slurm shedular can be found in the scripts folder by the name <a href="/scripts/data_download.sh">data_download.sh</a>.

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


### Preparing the Reference Sequence
The GATK needs two files when accessing the reference file:  
* A dictionary of the contig names and sizes
* An index file to access the reference fasta file bases

In here we are preparing these files upfront so the GATK will be able to use the FASTA file as a reference.

### Generate the BWA index  

First run the following bwa command to create the index, given that you have reference hg19.fa file already downloaded in to the folder called `hg19`  

<pre style="color: silver; background: black;">
hg19= < full_path_to >hg19.fa

module load bwa/0.7.17

cd ../hg19/
bwa index ${hg19}
</pre>
The full slurm script for creating the index can be found at scripts folder by the name, <a href="/scripts/bwa_index.sh">bwa_index.sh</a> .

This will create the following files:
<pre>
hg19/
├── hg19.fa
├── hg19.fa.amb
├── hg19.fa.ann
├── hg19.fa.bwt
├── hg19.fa.pac
└── hg19.fa.sa
</pre>


### Generate Fasta File Index  
Using `samtools` we will create a index of the reference fasta file.  
<pre>
module load samtools/1.7
samtools faidx ${hg19}
</pre>
The full slurm script for creating the index can be found at scripts folder by the name, <a href="/scripts/fasta_index.sh">fasta_index.sh</a>.
This will create:
<pre>
hg19/
└── hg19.fa.fai
</pre>

It will consist of one record per line for each of the contigs in the fasta file. Where each record is composed of 
* contig name
* size
* location
* bases per lane
* bytes per lane


### Create Sequence Dictionary
Use the picard tools to create the dictionary by the following command:
<pre>
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

java -jar $PICARD CreateSequenceDictionary \
        REFERENCE=${hg19} \
        OUTPUT=hg19.dict \
        CREATE_INDEX=True
</pre>
The full slurm script for creating the dictionary can be found at scripts folder by the name, <a href="/scripts/dictionary.sh">dictionary.sh</a>.

This will create:
<pre>
hg19/
└── hg19.dict
</pre>

This is formated like a SAM file header and when running GATK it automatically looks for these files.

### Aligning of Reads

Using BWA aligner we are going to align the reads to the reference fasta file.  

Since we have paired-end reads command will look like:
<pre>
d1=raw_data
d2=align

if [ ! -d ../${d2} ]; then
        mkdir -p ../${d2}
fi

cd ../${d2}

module load bwa/0.7.17
bwa mem -t 8 ${hg19} ../${d1}/${INPUT_FILE_NAME}${R1} ../${d1}/${INPUT_FILE_NAME}${R2} -o ${INPUT_FILE_NAME}.sam
</pre>
The full scrip is called <a href="/scripts/align.sh">align.sh</a> and can be found in the scripts folder.

<pre>
Usage: bwa mem [options] reference.fasta read1.fa read2.fa -o outputname.sam
mem    The BWA-MEM algorithm performs local alignment  
-t     Number of threads  
</pre>

Alignment will create SAM files, once the alignment is done for all the samples we will end up with:
<pre>
<strong>align</strong>/
├── SRR1517848.sam
├── SRR1517878.sam
├── SRR1517884.sam
├── SRR1517906.sam
├── SRR1517991.sam
├── SRR1518011.sam
├── SRR1518158.sam
└── SRR1518253.sam
</pre>

### SAM to BAM Conversion and Remove Singletons

### SAM to BAM Conversion
The BWA aligner will create the aligned files which is in the SAM format (Sequence Alignment/Map format). SAM files are human readable and can be large files. For the easy processing through the programs we will convert these files to binary format which is the BAM format using SAMtools.  

So the code will look like:
<pre>module load samtools/1.7
samtools view -@ 8 -bS ${INPUT_FILE_NAME}.sam > ${INPUT_FILE_NAME}.bam </pre> 

<pre>Useage: samtools view [options] 
-@    number of treads
-b    output in BAM format
-S    input format auto detected
</pre>

This will create BAM format files:
<pre>
<strong>align/</strong>
├── SRR1517848.bam
├── SRR1517878.bam
├── SRR1517884.bam
├── SRR1517906.bam
├── SRR1517991.bam
├── SRR1518011.bam
├── SRR1518158.bam
└── SRR1518253.bam
</pre>


### Remove Singletons
The unmapped reads are called singletons. The samtools flags can be used to remove these tagged reads.  
Using the following command we will remove the singletons:
<pre>
module load samtools/1.7
samtools view -@ 8 -F 0x04 -b ${INPUT_FILE_NAME}.bam > ${INPUT_FILE_NAME}_filtered.bam
</pre>

This will produce BAM files:
<pre>
<strong>align/</strong>
├── SRR1517848_filtered.bam
├── SRR1517878_filtered.bam
├── SRR1517884_filtered.bam
├── SRR1517906_filtered.bam
├── SRR1517991_filtered.bam
├── SRR1518011_filtered.bam
├── SRR1518158_filtered.bam
└── SRR1518253_filtered.bam
</pre>

The full slurm script for SAM to BAM conversion and singletons removal can be found at scripts folder by the name <a href="/scripts/singletons.sh">singletons.sh</a>


### Sort BAM files using PICARD

Once the singletons removed the BAM files are sorted using PICARD tools. For GATK analysis the BAM files need to be correctly formatted as well. The correct formatting includes:
* It must be aligned
* It must be sorted in coordinate order
* It must list the read groups with sample names in the header
* Every read must belong to a read group.
* The BAM file must pass Picard ValidateSamFile validation


In order to pass the files into Picard for processing in this step we will sort the aligned reads in the coordinate order.

The command is as follows:
<pre>
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

java -jar $PICARD SortSam \
        INPUT=${INPUT_FILE_NAME}_filtered.bam \
        OUTPUT=${INPUT_FILE_NAME}_filtered_sort.bam \
        SORT_ORDER=coordinate \
        CREATE_INDEX=True
</pre>

The full script is called <a href="/scripts/picard_sort.sh">picard_sort.sh</a> and can be found at scripts folder.

This will create sorted BAM format files:
<pre>
<strong>align/</strong>
├── SRR1517848_filtered_sort.bam
├── SRR1517878_filtered_sort.bam
├── SRR1517884_filtered_sort.bam
├── SRR1517906_filtered_sort.bam
├── SRR1517991_filtered_sort.bam
├── SRR1518011_filtered_sort.bam
├── SRR1518158_filtered_sort.bam
└── SRR1518253_filtered_sort.bam
</pre>


#### Tip
After sorting your reads you can check whether the reads are sorted accordingly and does it have the `SO: coordinate` flag to satisfy the GATK requirements by using `samtools` command to check the Header:  
`samtools view -H SRR1517848_filtered_sort.bam`   

which will produce: 
<pre>
@HD	VN:1.5	SO:coordinate
@SQ	SN:chr1	LN:249250621
@SQ	SN:chr2	LN:243199373
@SQ	SN:chr3	LN:198022430
@SQ	SN:chr4	LN:191154276
@SQ	SN:chr5	LN:180915260
@SQ	SN:chr6	LN:171115067
@SQ	SN:chr7	LN:159138663
@SQ	SN:chrX	LN:155270560
@SQ	SN:chr8	LN:146364022
@SQ	SN:chr9	LN:141213431
@SQ	SN:chr10	LN:135534747
@SQ	SN:chr11	LN:135006516
@SQ	SN:chr12	LN:133851895
@SQ	SN:chr13	LN:115169878
@SQ	SN:chr14	LN:107349540
@SQ	SN:chr15	LN:102531392
@SQ	SN:chr16	LN:90354753
@SQ	SN:chr17	LN:81195210
@SQ	SN:chr18	LN:78077248
@SQ	SN:chr20	LN:63025520
@SQ	SN:chrY	LN:59373566
@SQ	SN:chr19	LN:59128983
@SQ	SN:chr22	LN:51304566
@SQ	SN:chr21	LN:48129895
@SQ	SN:chr6_ssto_hap7	LN:4928567
@SQ	SN:chr6_mcf_hap5	LN:4833398
. 
.
</pre>


### Remove PCR Duplicates using PICARD

During the sequencing the same DNA molecules can be sequenced multiple times resulting in duplicates. These reads should not be counted as information in variant detection. In this step we will mark the duplicate reads and will remove them.  

Following command will remove the duplicate reads from each sample file.
<pre>
d3=noduplicates
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
</pre>
The full script is called <a href="/scripts/markduplicates.sh">markduplicates.sh</a> and can be found at scripts folder.

This will result in duplicates removed BAM files which will be:
<pre>
<strong>noduplicates/</strong>
├── SRR1517848_metrics.txt
├── SRR1517848_nodup.bai
├── SRR1517848_nodup.bam
├── SRR1517878_metrics.txt
├── SRR1517878_nodup.bai
├── SRR1517878_nodup.bam
├── SRR1517884_metrics.txt
├── SRR1517884_nodup.bai
├── SRR1517884_nodup.bam
├── SRR1517906_metrics.txt
├── SRR1517906_nodup.bai
├── SRR1517906_nodup.bam
├── SRR1517991_metrics.txt
├── SRR1517991_nodup.bai
├── SRR1517991_nodup.bam
├── SRR1518011_metrics.txt
├── SRR1518011_nodup.bai
├── SRR1518011_nodup.bam
├── SRR1518158_metrics.txt
├── SRR1518158_nodup.bai
├── SRR1518158_nodup.bam
├── SRR1518253_metrics.txt
├── SRR1518253_nodup.bai
└── SRR1518253_nodup.bam
</pre>


### Add Read Group Information

In this section we will be adding meta data about the sample. Adding meta data is very important is downstream analysis of your data, and these information is visible to GATK analysis tools. In here we use the minimal read group information for the samples and some are important tags.

In the SAM/BAM file the read group information is indicated in @RG tag which signify the "read group".

* ID : globally unique string which identify this run. Usually this linked to the lane where the data was run.

* SM : associated name in the DNA sample. This will be the sample identifier and it is the most important tag. In GATK all the analysis is done by sample, and this will selects which sample group it will belong to.

* PL : platform used. eg: "Illumina", "Pacbio", "iontorrent"

* LB : an identifier of the library from this DNA was sequenced. This field is important for future reference and quality control. In the case of errors associated with DNA preparation, this will link the data to the laboratory preparation step.

* PU : platform unit identifier for the run. The generic identifier will allow to go back to the machine, time and where it was run. Usually this is a flowcell-barcode-lane unique identifier.

To learn more about SAM tools tags please refer the [SAM tools format](http://samtools.github.io/hts-specs/SAMv1.pdf).

The following Picard tools command will add the read group information to each sample.
<pre>
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

cd ../${d4}

java -jar $PICARD AddOrReplaceReadGroups \
        INPUT=../${d3}/${INPUT_FILE_NAME}_nodup.bam \
        OUTPUT=${INPUT_FILE_NAME}_rg.bam \
        RGID=group1 \
        RGSM=${INPUT_FILE_NAME} \
        RGPL=illumina \
        RGLB=1 \
        RGPU=barcode \
        CREATE_INDEX=True
</pre>

The full script is called <a href="/scripts/add_readgroups.sh">add_readgroups.sh</a> and can be found in scripts folder.
The above command will add reads groups to each sample and will created BAM files:
<pre>
readgroup/
├── SRR1517848_rg.bai
├── SRR1517848_rg.bam
├── SRR1517878_rg.bai
├── SRR1517878_rg.bam
├── SRR1517884_rg.bai
├── SRR1517884_rg.bam
├── SRR1517906_rg.bai
├── SRR1517906_rg.bam
├── SRR1517991_rg.bai
├── SRR1517991_rg.bam
├── SRR1518011_rg.bai
├── SRR1518011_rg.bam
├── SRR1518158_rg.bai
├── SRR1518158_rg.bam
├── SRR1518253_rg.bai
└── SRR1518253_rg.bam
</pre>

#### How to check the reads have read group information ?
You can do this by quick samtools and unix commands using:  
`samtools view -H SRR1517848_rg.bam | grep '^@RG'`  
which will give you:
<pre>@RG	ID:group1	LB:1	PL:illumina	SM:SRR1517848	PU:barcode</pre>

The presence of the `@RG` tags indicate the presence of read groups. Each read group has a `SM` tag, indicating the sample from which the reads belonging to that read group originate.

In addition to the presence of a read group in the header, each read must belong to one and only one read group. Given the following example reads.


### Reorder BAM file

In this step we will reorder the SAM/BAM file to match the contig ordering in the reference fasta file, as to determine the exact name matching of contigs. Reads which are mapped to contigs which are absent in the new reference file are rejected or dropped. This step can run faster if we provide the indexed BAM file.

The following command will reorder the BAM file using PICARD tools:
<pre>
module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

cd ../${d5}/

java -jar $PICARD ReorderSam \
        INPUT=../${d4}/${INPUT_FILE_NAME}_rg.bam \
        OUTPUT=${INPUT_FILE_NAME}_karyotype.bam \
        REFERENCE=${hg19} \
        CREATE_INDEX=True
</pre>

The full script is called <a href="/scripts/reorder.sh">reorder.sh</a> and can be found in scripts folder.
This will create karyotype BAM files:
<pre>
<strong>reorder/</strong>
├── SRR1517848_karyotype.bai
├── SRR1517848_karyotype.bam
├── SRR1517878_karyotype.bai
├── SRR1517878_karyotype.bam
├── SRR1517884_karyotype.bai
├── SRR1517884_karyotype.bam
├── SRR1517906_karyotype.bai
├── SRR1517906_karyotype.bam
├── SRR1517991_karyotype.bai
├── SRR1517991_karyotype.bam
├── SRR1518011_karyotype.bai
├── SRR1518011_karyotype.bam
├── SRR1518158_karyotype.bai
├── SRR1518158_karyotype.bam
├── SRR1518253_karyotype.bai
└── SRR1518253_karyotype.bam
</pre>


### Variant Calling 

In this step we will call the variants using HaplotypeCaller in GATK software. 

<pre>
module load GATK/4.0
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

cd ../${d6}/ 

gatk HaplotypeCaller \
        --reference ${hg19} \
        --input ../${d5}/${INPUT_FILE_NAME}_karyotype.bam \
        --output ${INPUT_FILE_NAME}_haplotype.vcf
</pre>

The full script is called <a href="/scripts/haplotypeCaller.sh">haplotypeCaller.sh</a> and can be found in scripts folder.
This creates a VCF file called ${INPUT_FILE_NAME}_haplotype.vcf, containing all the variant sites the HaplotypeCaller evaluate including both SNPs and Indels.
<pre>
<strong>haplotypes/</strong>
├── SRR1517848_haplotype.vcf
├── SRR1517848_haplotype.vcf.idx
├── SRR1517878_haplotype.vcf
├── SRR1517878_haplotype.vcf.idx
├── SRR1517884_haplotype.vcf
├── SRR1517884_haplotype.vcf.idx
├── SRR1517906_haplotype.vcf
├── SRR1517906_haplotype.vcf.idx
├── SRR1517991_haplotype.vcf
├── SRR1517991_haplotype.vcf.idx
├── SRR1518011_haplotype.vcf
├── SRR1518011_haplotype.vcf.idx
├── SRR1518158_haplotype.vcf
├── SRR1518158_haplotype.vcf.idx
├── SRR1518253_haplotype.vcf
└── SRR1518253_haplotype.vcf.idx
</pre>

