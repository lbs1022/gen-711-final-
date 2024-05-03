# FINAL PROJECT  
## Study Background  
All of the raw data was obtained through a study conducted by the Hubbard Genome Center at UNH 

[Here is an overview of the initial study](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6952671/)

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
mkdir ./fastqc_trimmed

FASTQC TRIMMED READS - ON POWERSHELL  
ssh ntm1021@ron.sr.unh.edu  
cd gen-711-final-project  
fastqc ./trimmed_reads/15_S2_L001_R1_001.fastq.gz -o ./fastqc_trimmed  
fastqc ./trimmed_reads/15_S2_L001_R2_001.fastq.gz -o ./fastqc_trimmed  
fastqc ./trimmed_reads/69_S8_L001_R1_001.fastq.gz -o ./fastqc_trimmed  
fastqc ./trimmed_reads/69_S8_L001_R2_001.fastq.gz -o ./fastqc_trimmed  

GENOMA ASSEMBLY  
tmux new -s assembly  
conda activate genomics  
cd ~/gen-711-final-project  
nohup spades.py -1 15_S2_L001_R1_001.fastq.gz -2 15_S2_L001_R2_001.fastq.gz -s unpaired-15_S2_L001_R1_001.fastq.gz -s unpaired-15_S2_L001_R2_001.fastq.gz -o 15-spades-assembly-default -t 24 &  
nohup spades.py -1 69_S8_L001_R1_001.fastq.gz -2 69_S8_L001_R2_001.fastq.gz -s unpaired-69_S8_L001_R1_001.fastq.gz -s unpaired-69_S8_L001_R2_001.fastq.gz -o 69-spades-assembly-default -t 24 &  


## Conclusion  
chat GPT can be very helpful   
