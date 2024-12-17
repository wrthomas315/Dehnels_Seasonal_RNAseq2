#!/bin/bash

#make a new directory for these in the data folder
mkdir ../data/rawreads

#run a script that loops through id chart and downloads each file
while read j; do
        echo $j
        prefetch -v ../data/rawreads/$j
        fastq-dump -I --split-files --outdir ../data/rawreads ../data/rawreads/$j".sra"
done <../data/SRA_IDS
