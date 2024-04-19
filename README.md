# FINAL PROJECT  
## Study Background  
other people did the initial study  

## Methods - I used this code   
PREP DIRECTORIES + OBTAIN DATA  
mkdir gen-711-final-project  
cd gen-711-final-project  
mkdir raw_reads  
cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R1_001.fastq.gz raw_reads  
cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R2_001.fastq.gz raw_reads  
cp /tmp/gen711_project_data/genome-assembly-fqs/777_S1_L001_R1_001.fastq.gz raw_reads  
cp /tmp/gen711_project_data/genome-assembly-fqs/777_S1_L001_R2_001.fastq.gz raw_reads  
mkdir fastqc_output  

FASTQC - ON POWERSHELL  
ssh ntm1021@ron.sr.unh.edu  
cd gen-711-final-project  
fastqc 69_S8_L001_R1_001.fastq.gz -o ../fastqc_output  
fastqc 69_S8_L001_R2_001.fastq.gz -o ../fastqc_output  
fastqc 777_S1_L001_R1_001.fastq.gz -o ../fastqc_output  
fastqc 777_S1_L001_R2_001.fastq.gz -o ../fastqc_output  

## Conclusion  
chat GPT can be very helpful   
