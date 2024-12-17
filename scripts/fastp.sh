#!/bin/bash

#make a directoy for your trimmed reads
mkdir ../data/trimmed

while read j; do
        echo $j
        fastp --in1 "../data/rawreads"$j"_R1_001.fastq" --in2 "../data/rawreads"$j"_R2_001.fastq" --out1 "../data/trimmed/"$j"_1.trimmed.fastq" --out2 "../data/trimmed/"$j"_2.trimmed.fastq"
done <../data/SRA_IDS


