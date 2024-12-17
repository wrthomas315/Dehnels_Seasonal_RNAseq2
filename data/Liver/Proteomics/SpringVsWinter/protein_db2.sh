#!/bin/bash

#Preprocessing for GSEA to rank genes and get gene symbolbs
#Subsequent GSEA run code on r; 

#first copy and paste three columns from Johns file  to txt file (1st column:-logpvalue,2nd:fold change, 3rd:protein accession names
#for loop that cylces through each entry and prints out in rank form (direction*-log(pvalue), with logfold change, gene symbol
rm RANKFORM_$1
while read j; do
	#clips off all protein accession numbers for three columns
	THREE=$(echo $j | awk 'BEGIN { FS = ";" } ; { print $1}')
	#clips off to single protein accession number so use as search patters
	PROT_A=$(echo $j | awk 'BEGIN { FS = ";" } ; { print $1}' | awk '{print $3}')
	#takes protein accession number and greps it to the protein gtf with following 80 lines, searches for gene name, clips off excess parts of gene name
	PROT_GS=$(cat ../../GCF_027595985.1_mSorAra2.pri_protein.gpff | grep $PROT_A -A80 | grep "gene=" | sort -u | awk 'BEGIN { FS = "=" } ; { print $2}' | sed 's/^"//' | sed 's/"$//')
	#print it to screen
	echo $PROT_A
	#figure out sign of change
	SIGN=$(echo $j | awk 'BEGIN { FS = ";" } ; { print $1}' | awk '{print $2}' | bc)
	#if statements in bash don't handle decimals well, so see if negative sign is there in line, if true then add it to pvalue 
	if [[ $(echo $SIGN | grep "-") -eq 0 ]]
	then
		echo $THREE $PROT_GS "+" >> RANKFORM_$1
	else
		echo "-"$THREE $PROT_GS "-" >> RANKFORM_$1
	fi
done <$1 
#
#Now you have a rank file, but for gsea, we do not want loci, so remove them
cat RANKFORM_$1 | grep -v " LOC" > NOLOCRANKFORM_$1
#
#And need to remove duplicate genes that are a byproduct of the original bioinformatics
#Get list of genes and how many times found
awk '{print $4}' ./NOLOCRANKFORM_$1 | sort | uniq -c | sort -nr > geneDUPlist.txt

#for loop to get rid of any greater than two
while read j; do
	#save how many times found
	VALUE=$(echo $j | awk '{print $1}')
	#save gene as value
	GENE=$(echo $j | awk '{print $2}')
	#if statement if greater than 1, if so remove from document
	if [[ $VALUE -gt 1 ]]
	then
		GENE2=$(echo " " $GENE " ")
		echo $GENE2
		cat NOLOCRANKFORM_$1 | grep -v " $GENE " > tmp
		cp tmp NOLOCRANKFORM_$1
	fi
done <geneDUPlist.txt
#can run for other seasons
