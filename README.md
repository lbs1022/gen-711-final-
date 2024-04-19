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

RUN TRIMMOMATIC  
tmux
conda activate genomics   
cd ~/gen-711-final-project/  
trim_scriptV2.sh ./raw_reads/15_S2_L001_R1_001.fastq.gz ./raw_reads/15_S2_L001_R2_001.fastq.gz  
trim_scriptV2.sh ./raw_reads69_S8_L001_R1_001.fastq.gz ./raw_reads/69_S8_L001_R2_001.fastq.gz
ctrl+b  
d  

## Conclusion  
chat GPT can be very helpful   
