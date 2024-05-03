# FINAL PROJECT  
## Study Background 
The UNH Tech Camp was the group that collected the data and did an initial study on them.
This camp was made of kids from grade 5 to 12.
John Caparso, Alex Goin, Amino Hussein, Tonya Kirichok, Ada Milhauser, Aakansh Mysore, Nana Suarez, Olivia Tatro, Logan Quiter & Yasmin Yan were all the students that worked on this study. 
This specific project was led by Joseph Sevigny, Steve Simpson, Kelley Thomas and Andrea de Assis.
All of the data was collected in the summer of 2022 from Acadia National Park in Maine. 
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

GENOME ASSEMBLY  
tmux new -s assembly  
conda activate genomics  
cd ~/gen-711-final-project  
nohup spades.py -1 15_S2_L001_R1_001.fastq.gz -2 15_S2_L001_R2_001.fastq.gz -s unpaired-15_S2_L001_R1_001.fastq.gz -s unpaired-15_S2_L001_R2_001.fastq.gz -o 15-spades-assembly-default -t 24 &  
nohup spades.py -1 69_S8_L001_R1_001.fastq.gz -2 69_S8_L001_R2_001.fastq.gz -s unpaired-69_S8_L001_R1_001.fastq.gz -s unpaired-69_S8_L001_R2_001.fastq.gz -o 69-spades-assembly-default -t 24 &  

GENOME ASSESMENT  
quast.py 15contigs.fasta -o quast_results  
conda activate busco   
busco -i 15contigs.fasta -m genome -o busco-results -l bacteria

GENOME ANNOTATION  
conda activate genomics  
nohup prokka 15contigs.fasta --outdir prokka_output --cpus 24 --mincontiglen 200 &  
grep -o "product=.*" prokka_output/PROKKA_*.gff | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt  

EXTRACT 16S rRNA SEQUENCES   
extract_sequences "16S ribosomal RNA" prokka_output/PROKKA_05032024.ffn > 16S_sequence.fasta

BLAST  
makeblastdb -in 15contigs.fasta -dbtype nucl -out 15contigs_db  
blastn -query 16S_sequence.fasta -db contigs_db -out 16S_vs_contigs_6.tsv -outfmt 6  
blob_blast.sh contigs.fasta  

READ MAPPING  
bwa index contigs.fasta  




## Conclusion  
chat GPT can be very helpful   
