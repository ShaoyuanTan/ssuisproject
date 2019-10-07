# Comprehensive assessment of Oxford Nanopore MinION sequencing for bacterial characterization
This repository contains bioinformatics methods and codes for our paper: Comprehensive assessment of Oxford Nanopore MinION sequencing for bacterial characterization. The main content of this paper is generation of S.suis whole genomes using oxford nanopore sequencing, and its comparison to Illumina sequencing.
 
**Raw reads data** was deposited in NCBI under accession ID:

Below are the main bioinformatic analysis methods.

## MinION reads assembly
Software version: Canu-1.6

#### ONT reads processing

#### canu assembly

## Illumina reads assembly

## Hybrid assembly using MinION reads and Illumina reads
cat illumina_reads.fastq > illumina.fastq #cat illumina reads into one fastq file
bwa index minion.contigs.fasta # index contigs file from MinION assembly
bwa mem minion.contigs.fasta illumina.fastq > hybrid.sam # mapping illumina reads to minion contigs
samtools view -b hybrid.sam > hybrid.bam # change sam file to bam file
samtools sort hybrid.bam -o hybrid_sorted.bam # sort bam file
samtools index hybrid_sorted.bam # index bam file
java -Xmx16G -jar pilon-1.22.jar --genome minion.contigs.fastq --bam hybrid_sorted.bam --outdir folder for hybrid assembly # generate a hybrid assembly
