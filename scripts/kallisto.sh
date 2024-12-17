#!/bin/bash

#make an index for your reference transcriptome
REF=../data/ref/GCF_000181275.1_SorAra2.0_rna.fna
IDX=../data/ref/GCF_000181275.1_SorAra2.0_rna.idx
kallisto index -i $IDX $REF

#mkdir ../data/kallisto_quant

while read j; do
        echo $j
        kallisto quant -i ../data/ref/GCF_000181275.1_SorAra2.0_rna.idx -o ../data/kallisto_quant/$j "../data/trimmed/"$j"_1.trimmed.fastq" "../data/trimmed/"$j"_2.trimmed.fastq"
done <../data/SRA_IDS
