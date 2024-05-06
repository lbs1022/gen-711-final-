# FINAL PROJECT  
## Study Background 
The UNH Tech Camp was the group that collected the data and did an initial study on them.
This camp was made of kids from grade 5 to 12.
John Caparso, Alex Goin, Amino Hussein, Tonya Kirichok, Ada Milhauser, Aakansh Mysore, Nana Suarez, Olivia Tatro, Logan Quiter & Yasmin Yan were all the students that worked on this study. 
This specific project was led by Joseph Sevigny, Steve Simpson, Kelley Thomas and Andrea de Assis.
All of the data was collected in the summer of 2022 from Acadia National Park in Maine. 
## Methods - I used this code 
<details>
  <summary>Prep Directories and Obtain Data</summary>
   - We made a final project directory called "gen-711-final-project"        
   - We then pulled the forwards and backwards reads for samples 69 and 15  
        <details> 
        <summary>code</summary>
  print(mkdir gen-711-final-project)      
  cd gen-711-final-project        
  mkdir raw_reads        
  cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R1_001.fastq.gz raw_reads        
  cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R2_001.fastq.gz raw_reads        
  cp /tmp/gen711_project_data/genome-assembly-fqs/777_S1_L001_R1_001.fastq.gz raw_reads        
  cp /tmp/gen711_project_data/genome-assembly-fqs/777_S1_L001_R2_001.fastq.gz raw_reads        
  mkdir fastqc_output      
        </details>
  
  
  
  
  mkdir gen-711-final-project    
  
  cd gen-711-final-project    
  
  mkdir raw_reads    
  
  cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R1_001.fastq.gz raw_reads    
  
  cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R2_001.fastq.gz raw_reads    
  
  cp /tmp/gen711_project_data/genome-assembly-fqs/777_S1_L001_R1_001.fastq.gz raw_reads    
  
  cp /tmp/gen711_project_data/genome-assembly-fqs/777_S1_L001_R2_001.fastq.gz raw_reads    
  
  mkdir fastqc_output  
  
</details>


PREP DIRECTORIES + OBTAIN DATA  
mkdir gen-711-final-project  
cd gen-711-final-project  
mkdir raw_reads  
cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R1_001.fastq.gz raw_reads  
cp /tmp/gen711_project_data/genome-assembly-fqs/69_S8_L001_R2_001.fastq.gz raw_reads  
cp /tmp/gen711_project_data/genome-assembly-fqs/777_S1_L001_R1_001.fastq.gz raw_reads  
cp /tmp/gen711_project_data/genome-assembly-fqs/777_S1_L001_R2_001.fastq.gz raw_reads  
mkdir fastqc_output  

FASTQC   
cd raw_reads  
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
mkdir ./fastqc_trimmed

FASTQC TRIMMED READS  
cd gen-711-final-project  
fastqc ./trimmed_reads/15_S2_L001_R1_001.fastq.gz -o ./fastqc_trimmed  
fastqc ./trimmed_reads/15_S2_L001_R2_001.fastq.gz -o ./fastqc_trimmed  
fastqc ./trimmed_reads/69_S8_L001_R1_001.fastq.gz -o ./fastqc_trimmed  
fastqc ./trimmed_reads/69_S8_L001_R2_001.fastq.gz -o ./fastqc_trimmed  

GENOME ASSEMBLY    
conda activate genomics  
cd ~/gen-711-final-project  
nohup spades.py -1 15_S2_L001_R1_001.fastq.gz -2 15_S2_L001_R2_001.fastq.gz -s unpaired-15_S2_L001_R1_001.fastq.gz -s unpaired-15_S2_L001_R2_001.fastq.gz -o 15-spades-assembly-default -t 24 &  
nohup spades.py -1 69_S8_L001_R1_001.fastq.gz -2 69_S8_L001_R2_001.fastq.gz -s unpaired-69_S8_L001_R1_001.fastq.gz -s unpaired-69_S8_L001_R2_001.fastq.gz -o 69-spades-assembly-default -t 24 &  

GENOME ASSESMENT    
quast.py 15contigs.fasta -o quast_results_15    
quast.py 69contigs.fasta -o quast_results_69  
conda activate busco     
busco -i 15contigs.fasta -m genome -o busco-results-15 -l bacteria  
busco -i 69contigs.fasta -m genome -o busco-results-69 -l bacteria  

GENOME ANNOTATION    
conda activate genomics    
nohup prokka 15contigs.fasta --outdir prokka_output_15 --cpus 24 --mincontiglen 200 &     
nohup prokka 69contigs.fasta --outdir prokka_output_69 --cpus 24 --mincontiglen 200 &      
grep -o "product=.*" prokka_output_15/PROKKA_* | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances_15.txt     
grep -o "product=.*" prokka_output_69/PROKKA_* | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances_69.txt    

EXTRACT 16S rRNA SEQUENCES   
extract_sequences "16S ribosomal RNA" prokka_output_15/PROKKA_05032024.ffn > 16S_sequence_15.fasta    
extract_sequences "16S ribosomal RNA" prokka_output_69/PROKKA_05032024.ffn > 16S_sequence_69.fasta  

BLAST  
makeblastdb -in 15contigs.fasta -dbtype nucl -out 15contigs_db        
makeblastdb -in 69contigs.fasta -dbtype nucl -out 69contigs_db  
blastn -query 16S_sequence_15.fasta -db 15contigs_db -out 16S_vs_15contigs_6.tsv -outfmt 6        
blastn -query 16S_sequence_69.fasta -db 69contigs_db -out 16S_vs_69contigs_6.tsv -outfmt 6  
blob_blast.sh contigs.fasta    
blob_blast.sh contigs.fasta  

READ MAPPING  
bwa index 15contigs.fasta         
bwa index 69contigs.fasta     
bwa mem -t 24 15contigs.fasta ./trimmed_reads/15_S2_L001_R1_001.fastq.gz ./trimmed_reads/15_S2_L001_R2_001.fastq.gz > 15raw_mapped.sam      
bwa mem -t 24 69contigs.fasta ./trimmed_reads/69_S8_L001_R1_001.fastq.gz ./trimmed_reads/69_S8_L001_R2_001.fastq.gz > 69raw_mapped.sam    
samtools view -@ 24 -Sb  15raw_mapped.sam  | samtools sort -@ 24 -o 15sorted_mapped.bam  
samtools view -@ 24 -Sb  69raw_mapped.sam  | samtools sort -@ 24 -o 69sorted_mapped.bam  
samtools flagstat 15sorted_mapped.bam    
samtools flagstat 69sorted_mapped.bam    
samtools index 15sorted_mapped.bam    
samtools index 69orted_mapped.bam      
bedtools genomecov -ibam 15sorted_mapped.bam > 15coverage.out    
bedtools genomecov -ibam 69sorted_mapped.bam > 69coverage.out     
gen_input_table.py  --isbedfiles 15contigs.fasta 15coverage.out >  15coverage_table.tsv  
gen_input_table.py  --isbedfiles 69contigs.fasta 69coverage.out >  69coverage_table.tsv 

NONTARGET CONTIG REMOVAL  
blobtools create -i 15contigs.fasta -b 15sorted_mapped.bam -t 15contigs.fasta.vs.nt.cul5.1e5.megablast.out -o 15blob_out    
blobtools create -i 69contigs.fasta -b 69sorted_mapped.bam -t 69contigs.fasta.vs.nt.cul5.1e5.megablast.out -o 69blob_out  
blobtools view -i 15blob_out.blobDB.json -r all -o 15blob_taxonomy    
blobtools view -i 69blob_out.blobDB.json -r all -o 69blob_taxonomy  
blobtools plot -i 15blob_out.blobDB.json -r genus    
blobtools plot -i 69blob_out.blobDB.json -r genus

FILTERING GENOME  
#in order to decide the filtering criteria, we first tested the outcomes of filtering with different lengths and coverage numbers using the below code:  
grep -v '#' 15blob_taxonomy.15blob_out.blobDB.table.txt | awk -F'\t' '$2 > <insert-legnth>' | awk -F'\t' '$5 < <insert-coverage>' | awk '{print $18}'     
grep -v '#' 69blob_taxonomy.69blob_out.blobDB.table.txt | awk -F'\t' '$2 > <insert-legnth>' | awk -F'\t' '$5 < <insert-coverage>' | awk '{print $18}'  
#once we decided on the filtering criteria, we ran the following code:  
grep -v '##' 15blob_taxonomy.15blob_out.blobDB.table.txt | awk -F'\t' '$2 > 500' | awk -F'\t' '$5 > 20' | awk -F'\t' '{print $1}' > 15list_of_contigs_to_keep_len500_cov20.txt    
grep -v '##' 69blob_taxonomy.69blob_out.blobDB.table.txt | awk -F'\t' '$2 > 600' | awk -F'\t' '$5 > 10' | awk -F'\t' '{print $1}' > 69list_of_contigs_to_keep_len600_cov10.txt  





## Conclusion  
chat GPT can be very helpful   
