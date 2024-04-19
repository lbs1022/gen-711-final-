#!/bin/bash

FORWARD=$1
REVERSE=$2
THREADS=$3
ADAPTERS='/home/genome/joseph7e/scripts/adapters.fa'

source activate

mkdir trimmed-reads/

/opt/anaconda/anaconda3/bin/trimmomatic PE -threads $THREADS $FORWARD $REVERSE trimmed-reads/$FORWARD trimmed-reads/unpaired-$FORWARD trimmed-reads/$REVERSE trimmed-reads/unpaired-$REVERSE ILLUMINACLIP:$ADAPTERS:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:36
