#! /bin/bash

#Prep directories and obtained data
mkdir gen-711-final-project
cd gen-711-final-project
mkdir raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R1_001.fastq.gz raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R2_001.fastq.gz raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/15_S2_L001_R1_001.fastq.gz raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/15_S2_L001_R2_001.fastq.gz raw_reads
mkdir fastqc_output

#Run FastQC on raw reads
fastqc raw_reads/69_S8_L001_R1_001.fastq.gz -o ./fastqc_output
fastqc raw_reads/69_S8_L001_R2_001.fastq.gz -o ./fastqc_output
fastqc raw_reads/15_S2_L001_R1_001.fastq.gz -o ./fastqc_output
fastqc raw_reads/15_S2_L001_R2_001.fastq.gz -o ./fastqc_output

#Run trimmomatic
conda activate genomics 
trim_scriptV2.sh ./raw_reads/69_S8_L001_R1_001.fastq.gz .raw_reads/69_S8_L001_R2_001.fastq.gz
trim_scriptV2.sh ./raw_reads/15_S2_L001_R1_001.fastq.gz .raw_reads/15_S2_L001_R2_001.fastq.gz
mkdir ./fastqc_trimmed

#Run FastQC on trimmed reads
fastqc trimmed_reads/69_S8_L001_R1_001.fastq.gz -o ./fastqc_trimmed
fastqc trimmed_reads/69_S8_L001_R2_001.fastq.gz -o ./fastqc_trimmed
fastqc trimmed_reads/15_S2_L001_R1_001.fastq.gz -o ./fastqc_trimmed
fastqc trimmed_reads/15_S2_L001_R2_001.fastq.gz -o ./fastqc_trimmed

#Genome assembly 
cd gen-711-final-project/trimmed_reads
nohup spades.py -1 69_S8_L001_R1_001.fastq.gz -2 69_S8_L001_R2_001.fastq.gz -s unpaired-69_S8_L001_R1_001.fastq.gz -s unpaired-69_S8_L001_R2_001.fastq.gz -o 15-spades-assembly-default -t 24 &
nohup spades.py -1 15_S2_L001_R1_001.fastq.gz -2 15_S2_L001_R2_001.fastq.gz -s unpaired-15_S2_L001_R1_001.fastq.gz -s unpaired-15_S2_L001_R2_001.fastq.gz -o 15-spades-assembly-default -t 24 &
