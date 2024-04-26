#! /bin/bash

#Prepped directories and obtained data#
mkdir gen-711-final-project
cd gen-711-final-project
mkdir raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R1_001.fastq.gz raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R2_001.fastq.gz raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/15_S2_L001_R1_001.fastq.gz raw_reads
cp /tmp/gen711_project_data/genome-assembly-fqs/15_S2_L001_R2_001.fastq.gz raw_reads
mkdir fastqc_output

#Run FastQC on raw reads#
cd gen-711-final
