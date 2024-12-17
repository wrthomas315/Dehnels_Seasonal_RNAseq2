# Dehnels_Seasonal_RNAseq
Scripts and code to reproduce RNAseq analysis for looking at changes in expression through Dehnel's phenomenon; specifically in the cortex, hippocampus, and liver.

###Goals and strategy

The objective of this project is to evaluate how expression changes 
throughout the size plasticity of Dehnel’s phenomenon in *Sorex araneus*. 
We will be analyzing three different tissue types; the liver, hippocampus, 
and cortex. These regions will help to understand both size and metabolic 
changes underpinning Dehnel’s phenomenon. First, we will have to access 
the quality of our RNA-seq data, filter low quality reads and trim 
adapters, map to the transcriptome and quantify abundance, followed by 
normalization. Then we will 1) analyze differential expression between 
stages of Dehnel’s phenomenon using DESeq2, 2) characterize temporal 
patterns in metabolite concentrations, 3) build gene correlation networks 
and identify correlation between network structure and traits, 4) analyze differential conentrations of proteins and compare with transcriptomics. 
Throughout the analysis, we will look at resultant genes and test whether they enrich 
KEGG pathways using fgsea.

### Data

RNA-seq analyses require alignment to a reference and quantification of 
reads. The genome and original unfiltered reads can be downloaded as 
described below. However, these steps could be skipped when reproducing, 
as count data has been saved in ./data/Liver/GeneCounts.

The reference (mSorAra2; GCF_027595985.1) can be download from straight 
from NCBI, or using the code below.

```
mkdir ./data/ref/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/595/985/GCF_027595985.1_mSorAra2.pri/GCF_027595985.1_mSorAra2.pri_genomic.gff.gz ./data/ref/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/595/985/GCF_027595985.1_mSorAra2.pri/GCF_027595985.1_mSorAra2.pri_rna.fna.gz   ./data/ref/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/027/595/985/GCF_027595985.1_mSorAra2.pri/GCF_027595985.1_mSorAra2.pri_protein.gpff.gz   ./data/ref/
gunzip ../ref/GCF_027595985.1_mSorAra2.pri_rna.fna.gz  
```

RNA-seq data from this project can also be found on NCBI Sequencing Read Archive. The list of samples and associated accession numbers can be found in the data folder. These can be downloaded manually, or using the getter.sh script with the help of sratoolkit (https://github.com/ncbi/sra-tools). Note, all scripts are meant to be ran from the scripts folder with indirect paths contained in this git.

```
bash get_rawseq.sh
```

### Quality control, filtering, trimming
Again, these scripts can be skipped if reproducing from counts. If not proceed! Here we will trim adapters from our reads and remove low quality reads using default settings and fastp. Will need to download fastp to your local environment (https://github.com/OpenGene/fastp).

```
bash fastp.sh
```

### Mapping and quantification
Reads that have went through quality control are then mapped to the reference transcriptome and quantified using pseudoalignment. This method does not directly map reads to the genome, but can infer counts despite similarities between different coding regions (https://pachterlab.github.io/kallisto/about).

```
bash kallisto.sh
```
Note: This will create new transcript abundances separate ffrom the ones used in this analysis. Further scripts will use the ones I generated ./data/TISSUE/TranscriptAbundances and naming convention, but feel free to update the paths in the scripts with the ones you generated.

### Prepare proteomics results
Need to turn protein IDs to coding gene names to run fgsea and compare with transcriptomics results. Run for both comparisons and alter for different folders outputs.

```
bash protein_db2.sh ~/data/Liver/Proteomics/SpringVsWinter/JN_prot_liver_WiSp.txt
```

### Analyses
Each analysis was conducted using the R code below for each tissue type. For best results, run in RStudio, as each matrix and figure is not set to print out in a best attempt to not overwrite results. If this is your desired outcome, edit code to include saving.
 
```
R Dehnel_Liver.R
```

### MetaboAnalyst5.0
Both the above programs were done online at the below links. In a perfect world these should be scripted, however, due to conflicts in packages and Rversions they were not. https://www.metaboanalyst.ca https://david.ncifcrf.gov/summary.jsp
