# Comprehensive assessment of Oxford Nanopore MinION sequencing for bacterial characterization

This repository contains bioinformatics methods and codes for our paper: Comprehensive assessment of Oxford Nanopore MinION sequencing for bacterial characterization. The main content of this paper is generation of S.suis whole genomes using oxford nanopore sequencing, and its comparison to Illumina sequencing.
 
**Raw reads data** was deposited in NCBI under accession ID:

Below are the main bioinformatic analysis methods.

## MinION reads assembly

#### ONT raw reads processing and evaluation

MinION raw sequencing reads in fast5 format were basecalled using Albacore to generat fastq files. 

MinIONQC (https://github.com/roblanf/minion_qc) was used to determine yields, read length, and read quality.

Alignment of the MinION raw fastq file from basecall to the draft genome from Illumina sequencing was performed using Graphmap (https://github.com/isovic/graphmap).

AlignQC (https://github.com/jason-weirather/AlignQC) was used to evaluate the alignment and to determine the single read error rates.

Qualimap (http://qualimap.bioinfo.cipf.es/) was used to characterize the coverage distribution of the alignment across the whole genome.

#### Canu assembly

Software and version: Canu-1.6 (https://github.com/marbl/canu)

```
canu -p ssuisX -d ssuis__assembly genomeSize=2m -nanopore-raw ssuisX.fastq useGrid=0 # using MinION raw reads ssuisX.fastq to generate assembly ssuisX.contigs.fasta.
```

#### MinION assembly optimization

Optimization of de novo assembly was performed by examining different fold coverage levels, adjusting read quality filters, and applying different error correction methods. 

Different fold coverage subsets were obtained by randomly selecting reads from a single sequencing dataset using Seqtk (github.com/lh3/seqtk) and fastq-tools (github.com/dcjones/fastq-tools). 

Different quality filter cut-off subsets were generated using the "sequencing_summary.txt" from basecall and R (version 3.4.0).

Error correction methods included: racon-1.3.1 (https://github.com/isovic/racon), nanopolish-0.10.1 (https://github.com/jts/nanopolish) and Pilon-1.22 (https://github.com/broadinstitute/pilon).

**Racon:**

```
graphmap/bin/Linux-x64/graphmap align -r MinION.contigs.fasta -d MinION.rawreads.fastq -o graphmap.sam 

racon --sam MinION.rawreads.fastq graphmap.sam MinION.contigs.fasta racon_consensus.fasta
```

**Nanopolish:**

```
nanopolish index -d fast5_files/ reads.fasta

graphmap/bin/Linux-x64/graphmap align -r MinION.contigs.fasta -d MinION.rawreads.fastq -o graphmap.sam 

samtools view -b graphmap.sam > graphmap.bam

samtools sort graphmap.bam -o graphmap_sorted.bam

samtools index graphmap_sorted.bam

nanopolish variants --consensus -o polished.vcf -w "tig00000001:0-20000" -r reads.fasta -b reads.sorted.bam -g contigs.fasta

nanopolish/nanopolish vcf2fasta --skip-checks -g contigs.fasta dnanopolished.vcf > polished_genome.fasta
```

**Pilon:** Refer to "hybrid assembly" section

## Illumina reads processing and assembly

Software and version: Trimmomatic-0.36 (http://www.usadellab.org/cms/index.php?page=trimmomatic);  FastQC-0.11.6 (https://github.com/s-andrews/FastQC); SPAdes-3.11.1 (https://github.com/ablab/spades)

```
java -jar trimmomatic-0.36.jar PE -phred33 read_R1_001.fastq read_R2_001.fastq output_forward_paired.fastq output_forward_unpaired.fastq output_reverse_paired.fastq output_reverse_unpaired.fastq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 #trim Illumina reads using phred quality 33 as cutoff.

fastqc -o folder_for_quality output_forward_paired.fastq output_reverse_paired.fastq # check the trimmed reads quality

python spades.py -1 output_forward_paired.fastq -2 output_reverse_paired.fastq -o folder_for_assembly
```

## Hybrid assembly using MinION reads and Illumina reads

Software and version: Samtools-1.8 (https://github.com/samtools/samtools); pilon-1.22 (https://github.com/broadinstitute/pilon)

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

Assembly contigs from MinION, Illumina and hybrid were assessed by Quast to generat N50, number of contigs, longest contigs etc. And Consensus accuracy and assembly quality of each were evaluated using Mummer. 

Software and version: Quast-4.5 (https://github.com/ablab/quast); Mummer-4.0.0 (https://github.com/mummer4/mummer)

```
quast-4.5/quast.py assembly.fasta -o assembly_qualityaccessment 

Mummer-4.00beta2/nucmer -c 100 contigs.fasta reference.fasta -p comparison

Mummer-4.00beta2/dnadiff -d comparison.delta -p comparison_results 
```

## Multilocus sequence typing (MLST) determination

A total of 7 housekeeping genes (aroA, cpn60, dpr, gki, mutS, recA, thrA) for S. suis MLST analysis were downloaded from the PubMLST website (pubmlst.org/). 

Software and version: ncbi_blast

```
makeblastdb -in housekeeping_genes.fasta -title mlst -dbtype nucl -out reference

blastn -db dreference -max_target_seqs 2000 -query contigs.fasta -evalue 1e-3 word_size 11 -outfmt 7 > mlst_results
```

According to the blast results, top match alleles were recorded, and MLST profiles were predicted using PubMLST (pubmlst.org/).
