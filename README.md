# Comprehensive assessment of Oxford Nanopore MinION sequencing for bacterial characterization
This repository contains bioinformatics methods and codes for our paper: Comprehensive assessment of Oxford Nanopore MinION sequencing for bacterial characterization. The main content of this paper is generation of S.suis whole genomes using oxford nanopore sequencing, and its comparison to Illumina sequencing.
 
**Raw reads data** was deposited in NCBI under accession ID:

Below are the main bioinformatic analysis methods.

## MinION reads assembly


#### ONT reads processing

#### Canu assembly

Software version: Canu-1.6

```
canu -p ssuisX -d ssuis__assembly genomeSize=2m -nanopore-raw ssuisX.fastq useGrid=0 # using MinION raw reads ssuisX.fastq to generate assembly ssuisX.contigs.fasta.
```
#### MinION assembly optimization

## Illumina reads processing and assembly

Software version: Trimmomatic-0.36;  FastQC-0.11.6; SPAdes-3.11.1

```
java -jar trimmomatic-0.36.jar PE -phred33 read_R1_001.fastq read_R2_001.fastq output_forward_paired.fastq output_forward_unpaired.fastq output_reverse_paired.fastq output_reverse_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 #trim Illumina reads using phred quality 33 as cutoff.

fastqc -o folder_for_quality output_forward_paired.fastq output_reverse_paired.fastq # check the trimmed reads quality

python spades.py -1 output_forward_paired.fastq -2 output_reverse_paired.fastq -o folder_for_assembly
```

## Hybrid assembly using MinION reads and Illumina reads

Software version: Samtools-1.8 ; pilon-1.22

```
cat illumina_reads.fastq > illumina.fastq #cat illumina reads into one fastq file

bwa index minion.contigs.fasta # index contigs file from MinION assembly

bwa mem minion.contigs.fasta illumina.fastq > hybrid.sam # mapping illumina reads to minion contigs

samtools view -b hybrid.sam > hybrid.bam # change sam file to bam file

samtools sort hybrid.bam -o hybrid_sorted.bam # sort bam file

samtools index hybrid_sorted.bam # index bam file

java -Xmx16G -jar pilon-1.22.jar --genome minion.contigs.fastq --bam hybrid_sorted.bam --outdir folder for hybrid assembly # generate a hybrid assembly
```

## Assembly quality accessment
Software version: Quast-4.5
```
quast-4.5/quast.py assembly.fasta -o assembly_qualityaccessment # assembly contigs from MinION, Illumina and hybrid were assessed by quast to generat N50, number of contigs, longest contigs etc. 
```
