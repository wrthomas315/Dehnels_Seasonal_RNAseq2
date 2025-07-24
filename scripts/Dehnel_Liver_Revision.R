## Here is the R code for processing the temporal transcriptomics of Dehnel's phenomenon liver
#Reads were preprocessed on Noctillio server, where they were trimmed using fastp and aligned/quantified using Kalisto
#First I want to get the quantification files onto R and convert transcript abundance to gene abundance

###First set up all your libraries
library(stringr)
library(readr)
library(assertr)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(regionReport)
library(rhdf5)
library(edgeR)
library(tibble)
library( "genefilter" )
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(TCseq)
library(cluster)
library(EnhancedVolcano)
library(fgsea)



## STEP1: Import quant files
#Create a transcript to gene key using the gff file for the shrew
sortxdb <- makeTxDbFromGFF("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/refs/GCF_027595985.1_mSorAra2.pri_genomic.gtf.gz")
k<- keys(sortxdb, keytype="TXNAME")
sortx2gene <- AnnotationDbi::select(sortxdb, k, "GENEID", "TXNAME")
sortx2gene <- sortx2gene[!duplicated(sortx2gene[,1]),]
sortx2gene <- na.omit(sortx2gene)
#gets rid of the XRs, which are some misc_rnas, and do not apppear to be associated with genes
#load tidyverse here as it does not work if loaded above
library(tidyverse)

#now import quant files previously created using kallisto
cs_liv_samples <- read.table("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/SampleList", header = T)
cs_liv_files <-file.path("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/TranscriptAbundances/", cs_liv_samples$Sample_name, "abundance.tsv")
names(cs_liv_files) <- paste0("sample_", cs_liv_samples$Sample_name)
all(file.exists(cs_liv_files))

#write to tables and can save with hashed out code below, but note will overwrite
cs_liv.count.tsv <- tximport(cs_liv_files, type = "kallisto", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
cs_liv.tpm.tsv <- tximport(cs_liv_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
write.table(cs_liv.tpm.tsv$abundance, "/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/GeneCounts/GeneCounts.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)



## STEP2: Graph Liver mass for panel A
SampleData_Rreadable <- read_csv("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/SampleData_Rreadable.csv")
livermass <- SampleData_Rreadable$Liver_Mass
bodymass <- SampleData_Rreadable$Body_Mass
stagemass <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
livermean <- rep(mean(livermass,na.rm=TRUE))
livermass_matrix <- as.data.frame(cbind(livermass,bodymass))
livermass_matrix$norm <- (livermass_matrix$livermass/livermass_matrix$bodymass)
livermass_matrix$rel_mass <- c(rep(median(livermass_matrix$norm[1:5],na.rm=TRUE),24))
indiv_livmat <- as.data.frame(cbind(livermass_matrix,stagemass))
mean(indiv_livmat$livermass[10:14])-mean(indiv_livmat$livermass[15:19])
#Stage3vs4 shows significant size change
t.test(indiv_livmat$livermass[10:14],indiv_livmat$livermass[15:19],)
#plot
livmassplot <-ggplot(indiv_livmat, aes(x = stagemass, y =  livermass, fill = stagemass)) +
  geom_boxplot(size=.4,show.legend = FALSE)+
  xlab("Stage")+
  ylab("Liver Mass (g)")+
  scale_fill_manual(values =  alpha(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2"),.80))+
  theme_classic()
ggsave("~/data/Liver/BodyMass.png", livmassplot,width = 4, height = 3.14, dpi =300,)



## STEP3: Normalization
#Begin by setting up design matrix
colnames(cs_liv.count.tsv$counts) <- c("Stg1_1","Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
liv_stages <- factor(c(cs_liv_samples$Run))
liv_organs <- factor(c(cs_liv_samples$Condition))
liv_full1 <- factor(c(cs_liv_samples$Run))
liv_sex <- factor(c(cs_liv_samples$Sex))
liv_stages_organ_frame <-cbind(as.data.frame(liv_stages),as.data.frame(liv_organs),as.data.frame(liv_full1),as.data.frame(liv_sex))
#Now to normalize our reads
dds_liv_all <- DESeqDataSetFromMatrix(round(cs_liv.count.tsv$counts), DataFrame(liv_stages_organ_frame), ~liv_sex + liv_full1)
mcols(dds_liv_all) <- cbind(mcols(dds_liv_all), row.names(cs_liv.count.tsv$counts))
rownames(dds_liv_all) <- row.names(cs_liv.count.tsv$counts)
dds_liv_all <- DESeq(dds_liv_all)
#And look at pcas and heatmaps of counts
vst_dds_liv_all <- vst(dds_liv_all)
pcaData_liv_all<- plotPCA(vst_dds_liv_all,intgroup=c("liv_stages","liv_organs"), ntop=3336, returnData=TRUE)
percentVar <- attr(pcaData_liv_all, "percentVar")
pca <- ggplot(pcaData_liv_all, aes(x = PC1, y = PC2, color = factor(liv_stages))) +
  geom_point(size = 6) +
  scale_color_manual(
    name = "Life Stage",  # This is the new legend title
    values = alpha(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"), 0.80),
    labels = c("Summer Juvenile", "Autumn Juvenile", "Winter Juvenile", "Spring Adult", "Summer Adult")) +
  theme_bw()
ggsave("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/Supplements/Figures/Supplementary_Figure_Transcriptome_PCA.pdf", plot = pca, width = 10, height = 6)
ggsave("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/Supplements/Figures/Supplementary_Figure_Transcriptome_PCA.png", pca,width = 6, height = 5, dpi =300,)
###LiverHeatmap###
livsampleDists <- dist(t(assay(vst_dds_liv_all)))
livsampleDistMatrix <- as.matrix(livsampleDists)
colnames(livsampleDistMatrix) <- NULL
#make the heatmap
feet<-pheatmap(livsampleDistMatrix, clustering_distance_rows=livsampleDists,
         clustering_distance_cols = livsampleDists, color = colorRampPalette(rev(brewer.pal(n = 9, name ="Reds")))(255))
ggsave("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/Supplements/Figures/Supplementary_Figure_Transcriptome_pheatmap.png", feet,width = 6, height = 5, dpi =300,)
#make a heatmap of genes now
livtopVarGenes <- head( order( rowVars( assay(vst_dds_liv_all) ), decreasing=TRUE ), 3336 )
heatmap.2( assay(vst_dds_liv_all)[ livtopVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
)



## STEP 4 : Differential Expression
#LivStage3vs1
liv13res <- results(dds_liv_all, contrast = c("liv_full1","Stage3","Stage1"))
liv13resSig <- subset(liv13res,liv13res$padj<.05)
liv13up <- subset(liv13resSig,(liv13resSig$log2FoldChange)>=0)
liv13down <- subset(liv13resSig,(liv13resSig$log2FoldChange)<=0)
write.table(liv13res, "~/data/Liver/DifferentialExp/Stage3vs1/Stage3vsStage1.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(liv13resSig, "~/data/Liver/DifferentialExp/Stage3vs1/Stage3vsStage1_Sig.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(liv13up, "~/data/Liver/DifferentialExp/Stage3vs1/Stage3vsStage1Upsiglog.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(liv13down, "~/data/Liver/DifferentialExp/Stage3vs1/Stage3vsStage1Downsiglog.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
#LivStage4vs3
liv34res <- results(dds_liv_all, contrast = c("liv_full1","Stage4","Stage3"))
liv34resSig <- subset(liv34res,liv34res$padj<.05)
liv34up <- subset(liv34resSig,(liv34resSig$log2FoldChange)>=0)
liv34down <- subset(liv34resSig,(liv34resSig$log2FoldChange)<=0)
write.table(liv34res, "~/data/Liver/DifferentialExp/Stage4vs3/Stage4vsStage3_DESeq.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(liv34resSig, "~/data/Liver/DifferentialExp/Stage4vs3/Stage4vsStage3_DESeq_Sig.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(liv34up, "~/data/Liver/DifferentialExp/Stage4vs3/Stage4vsStage3Upsiglog.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(liv34down, "~/data/Liver/DifferentialExp/Stage4vs3/Stage4vsStage3Downsiglog.txt", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)



###STEP 5: FSGEA for transcriptomic data sets
#Stage 3 vs 1
Stage3vsStage1_DESeq <- read_delim("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Stage3vs1/Stage3vsStage1_DESeq.tsv", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
hres <- Stage3vsStage1_DESeq %>% 
  dplyr::select(Gene, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene) %>% 
  summarize(stat=mean(stat))
hranks <- deframe(hres)
tibble(gene = "NDUFB10", rank = hranks["NDUFB10"])
fgsea_Stage3vsStage1_DESeq <- fgsea(pathways=gmtPathways("/Users/bill/ShrewProjects/github/Sorex_Genome2/data/000_miscFiles/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt"), hranks) %>% 
  as_tibble() %>% 
  arrange(padj)
paths <- gmtPathways("/Users/bill/ShrewProjects/github/Sorex_Genome2/data/000_miscFiles/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt")
fgsea_Stage3vsStage1_DESeqTidy <- fgsea_Stage3vsStage1_DESeq %>%
  as_tibble() %>%
  arrange(desc(NES))
supp_fgsea_a1<-plotEnrichment(paths[["KEGG_OXIDATIVE_PHOSPHORYLATION"]],
               hranks) + labs(title="KEGG_OXIDATIVE_PHOSPHORYLATION")
supp_fgsea_a2<-plotEnrichment(paths[["KEGG_GLYCOLYSIS_GLUCONEOGENESIS"]],
               hranks) + labs(title="KEGG_GLYCOLYSIS_GLUCONEOGENESIS")
supp_fgsea_a3<-plotEnrichment(paths[["KEGG_FATTY_ACID_METABOLISM"]],
               hranks) + labs(title="KEGG_FATTY_ACID_METABOLISM")
# Show in a nice table:
fgsea_Stage3vsStage1_DESeqTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()
livexp_Stage3vs1_fsgea <-ggplot(subset(fgsea_Stage3vsStage1_DESeqTidy,padj<0.05 ), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES>0)) +
  coord_flip() +
  labs(x="KEGGPathwayLiver", y="Normalized Enrichment Score") + 
  scale_fill_manual(values =  c("red", "blue"))+
  theme_bw()
livexp_Stage3vs1_fsgea
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/fsgea_SumvsWinter.png", livexp_Stage3vs1_fsgea,width = 7.1, height = 6, dpi =300,)
fgsea_Stage3vsStage1_DESeqTidy2 <- apply(fgsea_Stage3vsStage1_DESeqTidy,2,as.character)
write.table(fgsea_Stage3vsStage1_DESeqTidy2, file='~/data/Liver/DifferentialExp/Stage3vs1/fsgea_liverStage3vs1Tidy.tsv', quote=FALSE, sep='\t')
fsgea1 <- read_delim("ShrewProjects/github/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Stage3vs1/fsgea_liverStage3vs1Tidy.tsv",
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
fsgea1 <- fsgea1 %>%
  mutate(
    # Remove any prefix like "83\t" and keep only the R-style c(...) string
    leadingEdge_clean = str_extract(leadingEdge, "c\\(.*\\)"),
    
    # Parse the cleaned R-style vector string into an actual R list
    leadingEdge_list = map(leadingEdge_clean, ~ tryCatch(eval(parse(text = .x)), error = function(e) NA))
  )
topPathwaysUp <- fsgea1 %>%
  filter(ES > 0) %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  pull(pathway)

topPathwaysDown <- fsgea1 %>%
  filter(ES < 0) %>%
  filter(padj < 0.05) %>%
  arrange(pval) %>%
  pull(pathway)

topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(paths[topPathways], hranks, fsgea1, 
              gseaParam=0.5)
# Open PDF device
pdf("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/Supplements/Figures/Supplementary_Figure_GSEA.pdf", width = 8.5, height = 11)

# First page: 3 enrichment plots
grid.arrange(supp_fgsea_a1, supp_fgsea_a2, supp_fgsea_a3, ncol = 1)

# New page for the table
grid.newpage()
grid.table(plotGseaTable(paths[topPathways], hranks, fsgea1, 
                         gseaParam=0.5))

# Close the device
dev.off()

#Stage 4 vs 3
Stage4vsStage3_DESeq <- read_delim("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Stage4vs3/Stage4vsStage3_DESeq.tsv", 
                                   delim = "\t", escape_double = FALSE, 
                                   trim_ws = TRUE)
hres2 <- Stage4vsStage3_DESeq %>% 
  dplyr::select(Gene, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(Gene) %>% 
  summarize(stat=mean(stat))
hranks2 <- deframe(hres2)
fgsea_Stage4vsStage3_DESeq <- fgsea(pathways=gmtPathways("~/data/refs/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt"), hranks2) %>% 
  as_tibble() %>% 
  arrange(padj)
fgsea_Stage4vsStage3_DESeqTidy <- fgsea_Stage4vsStage3_DESeq %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgsea_Stage4vsStage3_DESeqTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()
livexp_Stage4vs3_fsgea <-ggplot(subset(fgsea_Stage4vsStage3_DESeqTidy,padj<0.05), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES>0)) +
  coord_flip() +
  labs(x="KEGGPathwayLiver", y="Normalized Enrichment Score") + 
  scale_fill_manual(values =  c("red", "blue"))+
  theme_bw()
livexp_Stage4vs3_fsgea
ggsave("~/data/Liver/DifferentialExp/Stage4vs3/fsgea_SumvsWinter.png", livexp_Stage4vs3_fsgea,width = 6.5, height = 6, dpi =300,)
fgsea_Stage4vsStage3_DESeqTidy2 <- apply(fgsea_Stage4vsStage3_DESeqTidy,2,as.character)
write.table(fgsea_Stage4vsStage3_DESeqTidy2, file='~/data/Liver/DifferentialExp/Stage4vs3/fsgea_liverStage4vs3Tidy.tsv', quote=FALSE, sep='\t')



###STEP 6: FSGEA for proteomics data sets
#Stage 3 vs 1
Stage3vsStage1_PROT <- read_delim("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/data/Liver/Proteomics/WintervsSummerJ/NOLOCRANKFORM_NJ", 
                                  col_names = FALSE)
hres3 <- Stage3vsStage1_PROT %>% 
  dplyr::select(X4, X1) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(X4) %>% 
  summarize(stat=mean(X1))
hranks3 <- deframe(hres3)
fgsea_Stage3vsStage1_PROT <- fgsea(pathways=gmtPathways("~/data/refs/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt"), hranks3) %>% 
  as_tibble() %>% 
  arrange(padj)
fgsea_Stage3vsStage1_PROTTidy <- fgsea_Stage3vsStage1_PROT %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgsea_Stage3vsStage1_PROTTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()
livexp_Stage3vs1_fsgeaPROT <-ggplot(subset(fgsea_Stage3vsStage1_PROTTidy,padj<0.05 & (NES > 0 | NES < -1.55)), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES>0)) +
  coord_flip() +
  labs(x="KEGGPathwayLiver", y="Normalized Enrichment Score") + 
  scale_fill_manual(values =  c("red", "blue"))+
  theme_bw()
livexp_Stage3vs1_fsgeaPROT
ggsave("~/data/Liver/Proteomics/WintervsSummerJ/fsgea_WinterVsSummerJ_prot.png", livexp_Stage3vs1_fsgeaPROT,width = 6.5, height = 6, dpi =300,)
fgsea_Stage3vsStage1_PROTTidy2 <- apply(fgsea_Stage3vsStage1_PROTTidy,2,as.character)
write.table(fgsea_Stage3vsStage1_PROTTidy2, file='~/data/Liver/Proteomics/WintervsSummerJ/fsgea_liverStage3vs1Tidy_prot.tsv', quote=FALSE, sep='\t')

#Stage 4 vs 3
Stage4vsStage3_PROT <- read_delim("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/data/Liver/Proteomics/SpringVsWinter/NOLOCRANKFORM_JN_prot_liver_WiSp.txt", 
                                  col_names = FALSE)
hres4 <- Stage4vsStage3_PROT %>% 
  dplyr::select(X4, X1) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(X4) %>% 
  summarize(stat=mean(X1))
hranks4 <- deframe(hres4)
fgsea_Stage4vsStage3_PROT <- fgsea(pathways=gmtPathways("~/data/refs/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt"), hranks4) %>% 
  as_tibble() %>% 
  arrange(padj)
fgsea_Stage4vsStage3_PROTTidy <- fgsea_Stage4vsStage3_PROT %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgsea_Stage4vsStage3_PROTTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()
livexp_Stage4vs3_fsgeaPROT <-ggplot(subset(fgsea_Stage4vsStage3_PROTTidy,padj<0.05 & (NES > 0 | NES < -1.55)), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES>0)) +
  coord_flip() +
  labs(x="KEGGPathwayLiver", y="Normalized Enrichment Score") + 
  scale_fill_manual(values =  c("red", "blue"))+
  theme_bw()
livexp_Stage4vs3_fsgeaPROT
ggsave("~/data/Liver/Proteomics/SpringVsWinter/fsgea_SpringvsWinter_prot.png", livexp_Stage4vs3_fsgeaPROT,width = 6.5, height = 6, dpi =300,)
fgsea_Stage4vsStage3_PROTTidy2 <- apply(fgsea_Stage4vsStage3_PROTTidy,2,as.character)
write.table(fgsea_Stage4vsStage3_PROTTidy2, file='/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq3/data/Liver/Proteomics/SpringVsWinter/fsgea_liverStage4vs3Tidy_prot.tsv', quote=FALSE, sep='\t')





###STEP 7: Comparisons between proteomics and transcriptomics
#Stage 3 vs 1
#overlapping genes
overlap_gene_3vs1 <- list(Proteomics = Stage3vsStage1_PROT$X4, Transcriptomics=rownames(liv13res))
upset_gene_3vs1 <- fromList(overlap_gene_3vs1)
plot2_3vs1 <- upset(upset_gene_3vs1,nsets = 2,order.by = "freq")
plot2_3vs1
#overlapping direction and sigificance
#first subset by overlapping gene names
reduced_3vs1 <-Reduce(intersect, list(Stage3vsStage1_PROT$X4,rownames(liv13res)))
OLgene_3vs1 <- c()
dir_prot_3vs1 <- c()
dir_trans_3vs1 <- c()
sig_prot_3vs1 <- c()
sig_trans_3vs1 <- c()
for(i in 1:length(reduced_3vs1)) {
  OLgene_3vs1[i] <- reduced_3vs1[i]
  dir_prot_3vs1[i] <- subset(Stage3vsStage1_PROT,Stage3vsStage1_PROT$X4 == OLgene_3vs1[i])[1,2]
  dir_trans_3vs1[i] <- subset(liv13res,rownames(liv13res) == OLgene_3vs1[i])[1,2]
  sig_prot_3vs1[i] <- abs(subset(Stage3vsStage1_PROT,Stage3vsStage1_PROT$X4 == OLgene_3vs1[i])[1,1])
  sig_trans_3vs1[i] <- subset(liv13res,rownames(liv13res) == OLgene_3vs1[i])[1,6]
}
#overlapping significant
PT_combo_3vs1  <- as.data.frame(cbind(as.data.frame(OLgene_3vs1),as.data.frame(unlist(dir_prot_3vs1)),as.data.frame(dir_trans_3vs1),as.data.frame(unlist(sig_prot_3vs1)),as.data.frame(sig_trans_3vs1)))
colnames(PT_combo_3vs1) <- c("OLgene","dir_prot","dir_trans","sig_prot","sig_trans")
write.table(PT_combo_3vs1, file='~/data/Liver/Proteomics/WintervsSummerJ/DESEQplusPROTEO.tsv', quote=FALSE, sep='\t')
plot3_3vs1 <- ggplot(PT_combo_3vs1, aes(dir_prot,dir_trans))+
  geom_point(size =1)+
  theme_bw()+
  geom_smooth(method = "lm") +
  labs(x="Difference Proteomics", y="LFC Transcriptomics")+geom_hline(yintercept=0)+geom_vline(xintercept=0)
summary(lm(dir_prot ~ dir_trans, data = PT_combo_3vs1))
PT_combo_3vs1$plot_color <- with(PT_combo_3vs1, ifelse(
  sig_trans < 0.05 & sig_prot > -log10(0.05), "#CA9961",
  ifelse(sig_trans < 0.05, "#FD8C3E",
         ifelse(sig_prot >-log10(0.05), "#009284", "black"))
))
PT_combo_3vs1$plot_color <- factor(
  PT_combo_3vs1$plot_color,
  levels = c("black", "#FD8C3E", "#009284", "#CA9961")
)
PT_combo_3vs1 <- PT_combo_3vs1[order(PT_combo_3vs1$plot_color), ]
plot3_3vs1 <- ggplot(PT_combo_3vs1, aes(dir_prot, dir_trans)) +
  geom_point(size = 4.5,aes(color = plot_color)) +
  scale_color_identity()+
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(x = NULL, y = NULL) + 
  theme_bw() +
  theme(
    axis.text = element_blank(),
    panel.grid.major = element_line(linewidth = 2.5),
    panel.grid.minor = element_blank()
  )
ggsave("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Stage3vs1corr.png", plot3_3vs1,width = 12.5, height = 10, dpi =300,)

#Stage 4 vs 3
#overlapping genes
overlap_gene_4vs3 <- list(Proteomics = Stage4vsStage3_PROT$X4, Transcriptomics=rownames(liv34res))
upset_gene_4vs3 <- fromList(overlap_gene_4vs3)
plot2_4vs3 <- upset(upset_gene_4vs3,nsets = 2,order.by = "freq")
plot2_4vs3
#overlapping direction and sigificance
#first subset by overlapping gene names
reduced_4vs3 <-Reduce(intersect, list(Stage4vsStage3_PROT$X4,rownames(liv34res)))
OLgene_4vs3 <- c()
dir_prot_4vs3 <- c()
dir_trans_4vs3 <- c()
sig_prot_4vs3 <- c()
sig_trans_4vs3 <- c()
for(i in 1:length(reduced_4vs3)) {
  OLgene_4vs3[i] <- reduced_4vs3[i]
  dir_prot_4vs3[i] <- subset(Stage4vsStage3_PROT,Stage4vsStage3_PROT$X4 == OLgene_4vs3[i])[1,2]
  dir_trans_4vs3[i] <- subset(liv34res,rownames(liv34res) == OLgene_4vs3[i])[1,2]
  sig_prot_4vs3[i] <- abs(subset(Stage4vsStage3_PROT,Stage4vsStage3_PROT$X4 == OLgene_4vs3[i])[1,1])
  sig_trans_4vs3[i] <- subset(liv34res,rownames(liv34res) == OLgene_4vs3[i])[1,6]
}
#overlapping significant
PT_combo_4vs3  <- as.data.frame(cbind(as.data.frame(OLgene_4vs3),as.data.frame(unlist(dir_prot_4vs3)),as.data.frame(dir_trans_4vs3),as.data.frame(unlist(sig_prot_4vs3)),as.data.frame(sig_trans_4vs3)))
colnames(PT_combo_4vs3) <- c("OLgene","dir_prot","dir_trans","sig_prot","sig_trans")
write.table(PT_combo_4vs3, file='~/data/Liver/Proteomics/SpringVsWinter/DESEQplusPROTEO.tsv', quote=FALSE, sep='\t')
plot3_4vs3 <- ggplot(PT_combo_4vs3, aes(dir_prot,dir_trans))+
  geom_point()+
  theme_bw()+
  geom_smooth(method = "lm") +
  labs(x="Difference Proteomics", y="LFC Transcriptomics")+geom_hline(yintercept=0)+geom_vline(xintercept=0)
summary(lm(dir_prot ~ dir_trans, data = PT_combo_4vs3))
PT_combo_4vs3$plot_color <- with(PT_combo_4vs3, ifelse(
  sig_trans < 0.05 & sig_prot > -log10(0.05), "#CA9961",
  ifelse(sig_trans < 0.05, "#FD8C3E",
         ifelse(sig_prot >-log10(0.05), "#009284", "black"))
))
PT_combo_4vs3$plot_color <- factor(
  PT_combo_4vs3$plot_color,
  levels = c("black", "#FD8C3E", "#009284", "#CA9961")
)
PT_combo_4vs3 <- PT_combo_4vs3[order(PT_combo_4vs3$plot_color), ]
plot3_4vs3 <- ggplot(PT_combo_4vs3, aes(dir_prot, dir_trans)) +
  geom_point(size = 4.5,aes(color = plot_color)) +
  scale_color_identity()+
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  labs(x = NULL, y = NULL) + 
  theme_bw() +
  theme(
    axis.text = element_blank(),
    panel.grid.major = element_line(linewidth = 2.5),
    panel.grid.minor = element_blank()
  )
plot3_4vs3
ggsave("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Stage4vs3corr.png", plot3_4vs3,width = 12.5, height = 10, dpi =300,)

#Step 8 : Idenfity interesting overlapped genes
#Stage 3 vs 1 NDUFS
#NDUF Trans fsgea
which(fgsea_Stage3vsStage1_DESeqTidy$pathway == "KEGG_OXIDATIVE_PHOSPHORYLATION")
oxphos_fsgea <- unlist(fgsea_Stage3vsStage1_DESeqTidy$leadingEdge[1])
oxphos_fsgea_deseq <- liv13res[rownames(liv13res) %in% oxphos_fsgea, ]
oxphos_fsgea_deseqSig <- subset(oxphos_fsgea_deseq,padj<0.05)
#NDUF Proteomics
which(fgsea_Stage3vsStage1_PROTTidy$pathway == "KEGG_OXIDATIVE_PHOSPHORYLATION")
oxphos_fsgea_protlist <- unlist(fgsea_Stage3vsStage1_PROTTidy$leadingEdge[1])
oxphos_fsgea_prot <- Stage3vsStage1_PROT[Stage3vsStage1_PROT$X4 %in% oxphos_fsgea_protlist, ]
oxphos_fsgea_protSig <-subset(oxphos_fsgea_prot,X1>-log10(.05))
#overlapping genes in each set
Reduce(intersect, list(oxphos_fsgea_protlist,oxphos_fsgea))
#NDUF overlap
NDUF_sig_overlap <-Reduce(intersect, list(oxphos_fsgea_protSig$X4,rownames(oxphos_fsgea_deseqSig)))
#NDUF same direction + significant in one
NDUF_vector <- c(oxphos_fsgea_protSig$X4,rownames(oxphos_fsgea_deseqSig))
NDUF_frame <- PT_combo_3vs1[PT_combo_3vs1$OLgene %in% NDUF_vector, ]
#11/11 in same direction, no overlap in signficance

#Stage 3 vs 1 Gluconeogenesis
#Gluconeogenesis Trans fsgea
which(fgsea_Stage3vsStage1_DESeqTidy$pathway == "KEGG_GLYCOLYSIS_GLUCONEOGENESIS")
glyc_fsgea <- unlist(fgsea_Stage3vsStage1_DESeqTidy$leadingEdge[3])
glyc_fsgea_deseq <- liv13res[rownames(liv13res) %in% glyc_fsgea, ]
glyc_fsgea_deseqSig <- subset(glyc_fsgea_deseq,padj<0.05)
#Gluconeogenesis overlap, not generated from proteomics but from Sig list
glyc_fsgea_deseqSig_frame <- PT_combo_3vs1[PT_combo_3vs1$OLgene %in% rownames(glyc_fsgea_deseqSig), ]
#Gluconeogenesis just direction
glyc_frame <- PT_combo_3vs1[PT_combo_3vs1$OLgene %in% rownames(glyc_fsgea_deseq), ]
# 13/14 same

#Stage 3 vs 1 Fatty Acid metabolism
#FaM Trans fsgea
which(fgsea_Stage3vsStage1_DESeqTidy$pathway == "KEGG_FATTY_ACID_METABOLISM")
fattyacid_fsgea <- unlist(fgsea_Stage3vsStage1_DESeqTidy$leadingEdge[15])
fattyacid_fsgea_deseq <- liv13res[rownames(liv13res) %in% fattyacid_fsgea, ]
fattyacid_fsgea_deseqSig <- subset(fattyacid_fsgea_deseq,padj<0.05)
#FaM overlap
fattyacid_fsgea_deseqSig_frame <- PT_combo_3vs1[PT_combo_3vs1$OLgene %in% rownames(fattyacid_fsgea_deseqSig), ]
#FaM same direction
fam_frame <- PT_combo_3vs1[PT_combo_3vs1$OLgene %in% rownames(fattyacid_fsgea_deseq), ]
#11 of 18, more difficult to tell, But look at CPT1A to help be more conclusive!

#Complete overlap
#Significant (trans = 172, prot = 124 , overlap = 34, one in opposite directio DAG1 )
length(subset(PT_combo_3vs1,sig_trans<0.05)[,1])
length(subset(PT_combo_3vs1,sig_prot>-log10(.05) & dir_prot<0 )[,1])
length(subset(PT_combo_3vs1,sig_prot>-log10(.05) & sig_trans<0.05)[,1])
subset(PT_combo_3vs1,sig_prot>-log10(.05) & sig_trans<0.05)
#notable ACOX1, SLC25A20
#Same direction and one significant
length(subset(PT_combo_3vs1,(sig_prot>-log10(.05) | sig_trans<0.05) &((dir_prot>0 & dir_trans>0)|(dir_prot<0 & dir_trans<0)))[,1])
subset(PT_combo_3vs1,(sig_prot>-log10(.05) | sig_trans<0.05) &((dir_prot>0 & dir_trans>0)|(dir_prot<0 & dir_trans<0)))
#notable ACOX1, NDUFS6, ALDH3A2, NDUFB4

#Stage 4 vs 3
#Significant (trans = 274, prot = 183 , overlap = 41, but 10 in opposite direction)
length(subset(PT_combo_4vs3,sig_trans<0.05 )[,1])
length(subset(PT_combo_4vs3,sig_prot>-log10(.05) & dir_prot<0)[,1])
length(subset(PT_combo_4vs3,sig_prot>-log10(.05) & sig_trans<0.05)[,1])
subset(PT_combo_4vs3,sig_prot>-log10(.05) & sig_trans<0.05)
#notable None
#Same direction and one significant
length(subset(PT_combo_4vs3,(sig_prot>-log10(.05) | sig_trans<0.05) &((dir_prot>0 & dir_trans>0)|(dir_prot<0 & dir_trans<0)))[,1])
subset(PT_combo_4vs3,(sig_prot>-log10(.05) | sig_trans<0.05) &((dir_prot>0 & dir_trans>0)|(dir_prot<0 & dir_trans<0)))
subset(PT_combo_4vs3, OLgene == "ACOX1")
#notable NDUFB4, PGK1, ALDHA1A1,ALDH3A2, ACOX1


###Step 9: Plotting Heatmaps
#heatmap with overlap in proteomics and transcriptomics
#stage 1 vs 3
liv13resX <- results(dds_liv_all, contrast = c("liv_full1","Stage3","Stage1"))
#stage 3 vs 4
#make new copy of data
liv34resX <- results(dds_liv_all, contrast = c("liv_full1","Stage4","Stage3"))
#And plot with a volcano plot
#First change any read to a triangle if off the plot (>abs(6))
for (i in 1:length(liv13resX$padj)) {
  if  (liv13resX$padj[i]<1e-10 & !is.na (liv13resX$padj[i])) {
    liv13resX$padj[i] <- 1e-10
  }
  if (liv13resX$log2FoldChange[i]>4 & !is.na (liv13resX$log2FoldChange[i])) {
    liv13resX$log2FoldChange[i] <- 4
  }
  if (liv13resX$log2FoldChange[i]< -4 & !is.na (liv13resX$log2FoldChange[i])) {
    liv13resX$log2FoldChange[i] <- -4
  }
}
#create custom key-value pairs for different cell-types
#this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(liv13resX$log2FoldChange) == 4, 17,
  ifelse(liv13resX$padj==1e-10, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- 'PBMC'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
#And also change color for high and low effect, and not significant
keyvals <- ifelse(
  liv13resX$padj > 0.05, 'grey',
  ifelse(liv13resX$log2FoldChange <= -1.58, 'red',
         ifelse(liv13resX$log2FoldChange >= 1.58, 'blue',
                ifelse(liv13resX$log2FoldChange >= 0, 'blue',
                       'red'))))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'red'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
#Plot with genes identified in analysis of interet
threeone_overlap <-subset(PT_combo_3vs1,sig_prot>-log10(.05) & sig_trans<0.05 &((dir_prot>0 & dir_trans>0)|(dir_prot<0 & dir_trans<0)))
rownames(liv13resX[rownames(liv13resX) %in% threeone_overlap$OLgene, ])
top_up <- threeone_overlap %>%
  filter(dir_trans > 0) %>%
  arrange(sig_trans) %>%
  slice_head(n = 10)
top_down <- threeone_overlap %>%
  filter(dir_trans < 0) %>%
  arrange(sig_trans) %>%
  slice_head(n = 10)
top_combined <- bind_rows(
  top_up %>% mutate(direction = "up"),
  top_down %>% mutate(direction = "down")
)
Eh1<-EnhancedVolcano(liv13resX,
                     lab = rownames(liv13resX),
                     xlim=c(-4 ,4),
                     ylim=c(0,10),
                     x = 'log2FoldChange',
                     y = 'padj',
                     shapeCustom = keyvals.shape,
                     selectLab = top_combined$OLgene,
                     pCutoff = .05,
                     FCcutoff = 15,
                     colCustom = keyvals,
                     colAlpha = .5,
                     boxedLabels = TRUE,
                     pointSize = 4.0,
                     legendPosition = 'none',
                     drawConnectors = TRUE,
                     gridlines.major  = FALSE,
                     widthConnectors = 0.5)
Eh1<-EnhancedVolcano(
  liv13resX,
  lab = rownames(liv13resX),
  x = 'log2FoldChange',
  y = 'padj',
  xlab = "",
  ylab = "",
  xlim = c(-4, 4),
  ylim = c(0, 10),
  shapeCustom = keyvals.shape,
  selectLab = top_combined$OLgene,
  pCutoff = 0.00000000000001,
  FCcutoff = 15,
  colCustom = keyvals,
  colAlpha = 0.5,
  boxedLabels = TRUE,
  pointSize = 12.0,
  labSize = 12,
  drawConnectors = TRUE,
  widthConnectors = 0.8,
  legendPosition = 'none',
  gridlines.minor = FALSE,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  hline = 0.05,
  hlineWidth = 3
)+theme(
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  panel.grid.major = element_line(linewidth = 2.5)
)
Eh1
ggsave("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Stage3vs1VOLC_trans.png", Eh1,width = 12.5, height = 10, dpi =300,)
#Stage 3 vs 1 proteomics
###volcano
PT_combo_3vs1X <- PT_combo_3vs1
PT_combo_3vs1X$sig_prot <- 10^(-PT_combo_3vs1X$sig_prot)
PT_combo_3vs1X$sig_prot
for (i in 1:length(PT_combo_3vs1X$sig_prot)) {
  if  (PT_combo_3vs1X$sig_prot[i]<1e-3 & !is.na (PT_combo_3vs1X$sig_prot[i])) {
    PT_combo_3vs1X$sig_prot[i] <-1e-3
  }
  if (PT_combo_3vs1X$dir_prot[i]>3 & !is.na (PT_combo_3vs1X$dir_prot[i])) {
    PT_combo_3vs1X$dir_prot[i] <- 3
  }
  if (PT_combo_3vs1X$dir_prot[i]< -3 & !is.na (PT_combo_3vs1X$dir_prot[i])) {
    PT_combo_3vs1X$dir_prot[i] <- -3
  }
}
#create custom key-value pairs for different cell-types
#this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(PT_combo_3vs1X$dir_prot) == 3, 17,
  ifelse(PT_combo_3vs1X$sig_prot==1e-3, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- 'PBMC'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
#And also change color for high and low effect, and not significant
keyvals <- ifelse(
  PT_combo_3vs1X$sig_prot > 0.05, 'grey',
  ifelse(PT_combo_3vs1X$dir_prot <= -1.58, 'red',
         ifelse(PT_combo_3vs1X$dir_prot >= 1.58, 'blue',
                ifelse(PT_combo_3vs1X$dir_prot >= 0, 'blue',
                       'red'))))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'red'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
top_up <- threeone_overlap %>%
  filter(dir_prot > 0) %>%
  arrange(sig_prot) %>%
  slice_head(n = 10)
top_down <- threeone_overlap %>%
  filter(dir_prot < 0) %>%
  arrange(sig_prot) %>%
  slice_head(n = 10)
top_combined <- bind_rows(
  top_up %>% mutate(direction = "up"),
  top_down %>% mutate(direction = "down")
)
Eh2<-EnhancedVolcano(PT_combo_3vs1X,
                lab = PT_combo_3vs1X$OLgene,
                xlim=c(-3 ,3),
                ylim=c(0,3),
                x = 'dir_prot',
                y = 'sig_prot',
                shapeCustom = keyvals.shape,
                selectLab = c(top_combined$OLgene),
                pCutoff = .05,
                FCcutoff = 15,
                colCustom = keyvals,
                colAlpha = .5,
                boxedLabels = TRUE,
                pointSize = 4.0,
                legendPosition = 'none',
                drawConnectors = TRUE,
                gridlines.major  = FALSE,
                widthConnectors = 0.5,)
Eh2<-EnhancedVolcano(PT_combo_3vs1X,
                     lab = PT_combo_3vs1X$OLgene,
                     xlim=c(-3 ,3),
                     ylim=c(0,3),
                     x = 'dir_prot',
                     y = 'sig_prot',
                     xlab = "",
                     ylab = "",
                     shapeCustom = keyvals.shape,
                     selectLab = c(top_combined$OLgene),
                     FCcutoff = 15,
                     colCustom = keyvals,
                     colAlpha = .5,
                     boxedLabels = TRUE,
                     pointSize = 12.0,
                     labSize = 12,
                     drawConnectors = TRUE,
                     gridlines.minor  = FALSE,
                     widthConnectors = 0.8,
                     legendPosition = 'none',
                     title = NULL,
                     subtitle = NULL,
                     caption = NULL,
                     hline = 0.05,
                     hlineWidth = 3)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_line(linewidth = 2.5)
  )
Eh2
ggsave("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Stage3vs1VOLC_prot.png", Eh2,width = 12.5, height = 10, dpi =300,)
#stage 3 vs 4 transcriptomics
#make new copy of data
liv34resX <- results(dds_liv_all, contrast = c("liv_full1","Stage4","Stage3"))
#And plot with a volcano plot
#First change any read to a triangle if off the plot (>abs(6))
for (i in 1:length(liv34resX$padj)) {
  if  (liv34resX$padj[i]<1e-10 & !is.na (liv34resX$padj[i])) {
    liv34resX$padj[i] <- 1e-10
  }
  if (liv34resX$log2FoldChange[i]>3 & !is.na (liv34resX$log2FoldChange[i])) {
    liv34resX$log2FoldChange[i] <- 3
  }
  if (liv34resX$log2FoldChange[i]< -3 & !is.na (liv34resX$log2FoldChange[i])) {
    liv34resX$log2FoldChange[i] <- -3
  }
}
#create custom key-value pairs for different cell-types
#this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(liv34resX$log2FoldChange) == 3, 17,
  ifelse(liv34resX$padj==1e-10, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- 'PBMC'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
#And also change color for high and low effect, and not significant
keyvals <- ifelse(
  liv34resX$padj > 0.05, 'grey',
  ifelse(liv34resX$log2FoldChange <= -1.58, 'red',
         ifelse(liv34resX$log2FoldChange >= 1.58, 'blue',
                ifelse(liv34resX$log2FoldChange >= 0, 'blue',
                       'red'))))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'red'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
#Plot with genes identified in analysis of interet
fourthree_overlap <-subset(PT_combo_4vs3,sig_prot>-log10(.05) & sig_trans<0.05 &((dir_prot>0 & dir_trans>0)|(dir_prot<0 & dir_trans<0)))
top_up <- fourthree_overlap %>%
  filter(dir_trans > 0) %>%
  arrange(sig_trans) %>%
  slice_head(n = 10)
top_down <- fourthree_overlap %>%
  filter(dir_trans < 0) %>%
  arrange(sig_trans) %>%
  slice_head(n = 10)
top_combined <- bind_rows(
  top_up %>% mutate(direction = "up"),
  top_down %>% mutate(direction = "down")
)
Eh3<-EnhancedVolcano(liv34resX,
                lab = rownames(liv34resX),
                xlim=c(-3 ,3),
                ylim=c(0,10),
                x = 'log2FoldChange',
                y = 'padj',
                xlab = "",
                ylab = "",
                shapeCustom = keyvals.shape,
                selectLab = c(top_combined$OLgene),
                pCutoff = 0.00000000000001,
                FCcutoff = 15,
                colCustom = keyvals,
                colAlpha = .5,
                boxedLabels = TRUE,
                pointSize = 12.0,
                labSize = 12,
                drawConnectors = TRUE,
                gridlines.minor  = FALSE,
                widthConnectors = 0.8,
                legendPosition = 'none',
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                hline = 0.05,
                hlineWidth = 3)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_line(linewidth = 2.5)
  )
Eh3
ggsave("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Stage4vs3VOLC_trans.png", Eh3,width = 12.5, height = 10, dpi =300,)

#stage 3 vs 4 proteomics
PT_combo_4vs3X <- PT_combo_4vs3
PT_combo_4vs3X$sig_prot <- 10^(-PT_combo_4vs3X$sig_prot)
PT_combo_4vs3X$sig_prot
for (i in 1:length(PT_combo_4vs3X$sig_prot)) {
  if  (PT_combo_4vs3X$sig_prot[i]<1e-3 & !is.na (PT_combo_4vs3X$sig_prot[i])) {
    PT_combo_4vs3X$sig_prot[i] <-1e-3
  }
  if (PT_combo_4vs3X$dir_prot[i]>3 & !is.na (PT_combo_4vs3X$dir_prot[i])) {
    PT_combo_4vs3X$dir_prot[i] <- 3
  }
  if (PT_combo_4vs3X$dir_prot[i]< -3 & !is.na (PT_combo_4vs3X$dir_prot[i])) {
    PT_combo_4vs3X$dir_prot[i] <- -3
  }
}
#create custom key-value pairs for different cell-types
#this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(PT_combo_4vs3X$dir_prot) == 3, 17,
  ifelse(PT_combo_4vs3X$sig_prot==1e-3, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- 'PBMC'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
#And also change color for high and low effect, and not significant
keyvals <- ifelse(
  PT_combo_4vs3X$sig_prot > 0.05, 'grey',
  ifelse(PT_combo_4vs3X$dir_prot <= -1.58, 'red',
         ifelse(PT_combo_4vs3X$dir_prot >= 1.58, 'blue',
                ifelse(PT_combo_4vs3X$dir_prot >= 0, 'blue',
                       'red'))))
keyvals[is.na(keyvals)] <- 'grey'
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'red'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
top_up <- fourthree_overlap %>%
  filter(dir_prot > 0) %>%
  arrange(sig_prot) %>%
  slice_head(n = 10)
top_down <- fourthree_overlap %>%
  filter(dir_prot < 0) %>%
  arrange(sig_prot) %>%
  slice_head(n = 10)
top_combined <- bind_rows(
  top_up %>% mutate(direction = "up"),
  top_down %>% mutate(direction = "down")
)
Eh4<-EnhancedVolcano(PT_combo_4vs3X,
                lab = PT_combo_4vs3X$OLgene,
                xlim=c(-3 ,3),
                ylim=c(0,3),
                x = 'dir_prot',
                y = 'sig_prot',
                xlab = "",
                ylab = "",
                shapeCustom = keyvals.shape,
                selectLab = c(top_combined$OLgene),
                FCcutoff = 15,
                colCustom = keyvals,
                colAlpha = .5,
                boxedLabels = TRUE,
                pointSize = 12.0,
                labSize = 12,
                legendPosition = 'none',
                drawConnectors = TRUE,
                widthConnectors = 0.8,
                gridlines.minor  = FALSE,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                hline = 0.05,
                hlineWidth = 3)+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid.major = element_line(linewidth = 2.5)
  )
Eh4
ggsave("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Stage4vs3VOLC_prot.png", Eh4,width = 12.5, height = 10, dpi =300,)



###Step 10: PLotting Interesting Genes
interesting_liv <- DESeq2::counts(dds_liv_all, normalized=TRUE)
#GLUCOSE
#ALDH3A2
plotCounts(dds_liv_all, gene="ALDH3A2", intgroup="liv_full1")
ALDH3A2_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "ALDH3A2", ])
ALDH3A2_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(ALDH3A2_liv) <- c("Expression", "Stage")
#plot
summary_ALDH3A2 <- ALDH3A2_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_ALDH3A2
# Step 3: Plot the data
ALDH3A2_plot<-ggplot(summary_ALDH3A2 , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "lightskyblue1", size =3) +
  geom_point(size = 5.5, color = "lightskyblue1") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
ALDH3A2_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/ALDH3A2_plot.png", ALDH3A2_plot,width = 6, height = 3.7, dpi =300,)
#PGK1
plotCounts(dds_liv_all, gene="PGK1", intgroup="liv_full1")
PGK1_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "PGK1", ])
PGK1_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(PGK1_liv) <- c("Expression", "Stage")
#plot
summary_PGK1 <- PGK1_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_PGK1
# Step 3: Plot the data
PGK1_plot<-ggplot(summary_PGK1 , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "lightskyblue1", size =3) +
  geom_point(size = 5.5, color = "lightskyblue1") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
PGK1_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/PGK1_plot.png", PGK1_plot,width = 6, height = 3.7, dpi =300,)
#ALDH1B1
plotCounts(dds_liv_all, gene="ALDH1B1", intgroup="liv_full1")
ALDH1B1_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "ALDH1B1", ])
ALDH1B1_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(ALDH1B1_liv) <- c("Expression", "Stage")
#plot
summary_ALDH1B1 <- ALDH1B1_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_ALDH1B1
# Step 3: Plot the data
ALDH1B1_plot<-ggplot(summary_ALDH1B1 , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "lightskyblue1", size =3) +
  geom_point(size = 5.5, color = "lightskyblue1") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
ALDH1B1_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/ALDH1B1_plot.png", ALDH1B1_plot,width = 6, height = 3.7, dpi =300,)
#G6PC1
plotCounts(dds_liv_all, gene="G6PC1", intgroup="liv_full1")
G6PC1_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "G6PC1", ])
G6PC1_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(G6PC1_liv) <- c("Expression", "Stage")
#plot
summary_G6PC1 <- G6PC1_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_G6PC1
# Step 3: Plot the data
G6PC1_plot<-ggplot(summary_G6PC1 , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "lightskyblue1", size =3) +
  geom_point(size = 5.5, color = "lightskyblue1") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
G6PC1_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/G6PC1_plot.png", G6PC1_plot,width = 6.1, height = 3.7, dpi =300,)

#FAM
#ACOX1
plotCounts(dds_liv_all, gene="ACOX1", intgroup="liv_full1")
ACOX1_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "ACOX1", ])
ACOX1_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(ACOX1_liv) <- c("Expression", "Stage")
#plot
summary_ACOX1 <- ACOX1_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_ACOX1
# Step 3: Plot the data
ACOX1_plot<-ggplot(summary_ACOX1 , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "cornflowerblue", size =3) +
  geom_point(size = 5.5, color = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
ACOX1_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/ACOX1_plot.png", ACOX1_plot,width = 6.1, height = 3.7, dpi =300,)
#HADHA
plotCounts(dds_liv_all, gene="HADHA", intgroup="liv_full1")
HADHA_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "HADHA", ])
HADHA_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(HADHA_liv) <- c("Expression", "Stage")
#plot
summary_HADHA <- HADHA_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_HADHA
# Step 3: Plot the data
HADHA_plot<-ggplot(summary_HADHA , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "cornflowerblue", size =3) +
  geom_point(size = 5.5, color = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
HADHA_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/HADHA_plot.png", HADHA_plot,width = 6, height = 3.7, dpi =300,)
#ACADVL
plotCounts(dds_liv_all, gene="ACADVL", intgroup="liv_full1")
ACADVL_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "ACADVL", ])
ACADVL_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(ACADVL_liv) <- c("Expression", "Stage")
#plot
summary_ACADVL <- ACADVL_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_ACADVL
# Step 3: Plot the data
ACADVL_plot<-ggplot(summary_ACADVL , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "cornflowerblue", size =3) +
  geom_point(size = 5.5, color = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
ACADVL_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/ACADVL_plot.png", ACADVL_plot,width = 6.1, height = 3.7, dpi =300,)
#ACADM
plotCounts(dds_liv_all, gene="ACADM", intgroup="liv_full1")
ACADM_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "ACADM", ])
ACADM_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(ACADM_liv) <- c("Expression", "Stage")
#plot
summary_ACADM <- ACADM_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_ACADM
# Step 3: Plot the data
ACADM_plot<-ggplot(summary_ACADM , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "cornflowerblue", size =3) +
  geom_point(size = 5.5, color = "cornflowerblue") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
ACADM_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/ACADM_plot.png", ACADM_plot,width = 6, height = 3.7, dpi =300,)

#OXPHOS
#SDHB frame
plotCounts(dds_liv_all, gene="SDHB", intgroup="liv_full1")
SDHB_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "SDHB", ])
SDHB_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
SDHB_liv
colnames(SDHB_liv) <- c("Expression", "Stage")
#plot
summary_SDHB <- SDHB_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_SDHB
# Step 3: Plot the data
SDHB_plot<-ggplot(summary_SDHB , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "navy", size =3) +
  geom_point(size = 5.5, color = "navy") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/SDHB_plot.png", SDHB_plot,width = 6, height = 3.7, dpi =300,)
#NFUDB4
plotCounts(dds_liv_all, gene="NDUFB4", intgroup="liv_full1")
NFUDB4_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "NDUFB4", ])
NFUDB4_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(NFUDB4_liv) <- c("Expression", "Stage")
#plot
summary_NFUDB4 <- NFUDB4_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_NFUDB4
# Step 3: Plot the data
NFUDB4_plot<-ggplot(summary_NFUDB4 , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "navy", size =3) +
  geom_point(size = 5.5, color = "navy") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
NFUDB4_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/NFUDB4_plot.png", NFUDB4_plot,width = 6, height = 3.7, dpi =300,)
#PPA1
plotCounts(dds_liv_all, gene="PPA1", intgroup="liv_full1")
PPA1_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "PPA1", ])
PPA1_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(PPA1_liv) <- c("Expression", "Stage")
#plot
summary_PPA1 <- PPA1_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_PPA1
# Step 3: Plot the data
PPA1_plot<-ggplot(summary_PPA1 , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "navy", size =3) +
  geom_point(size = 5.5, color = "navy") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
PPA1_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/PPA1_plot.png", PPA1_plot,width = 6.1, height = 3.7, dpi =300,)
#ATP5MG
plotCounts(dds_liv_all, gene="ATP5MG", intgroup="liv_full1")
ATP5MG_liv <-as.data.frame(interesting_liv[rownames(interesting_liv) %in% "ATP5MG", ])
ATP5MG_liv$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
colnames(ATP5MG_liv) <- c("Expression", "Stage")
#plot
summary_ATP5MG <- ATP5MG_liv %>%
  group_by(Stage) %>%
  summarise(
    mean_expression = mean(Expression),
    sd_expression = sd(Expression)
  )
summary_ATP5MG
# Step 3: Plot the data
ATP5MG_plot<-ggplot(summary_ATP5MG , aes(x = Stage, y = mean_expression)) +
  geom_line(group = 1, color = "navy", size =3) +
  geom_point(size = 5.5, color = "navy") +
  geom_errorbar(aes(ymin = mean_expression - sd_expression, 
                    ymax = mean_expression + sd_expression), size = 1,
                width = 0.1) +
  theme_classic()
ATP5MG_plot
ggsave("~/data/Liver/DifferentialExp/Stage3vs1/Figure2plots/ATP5MG_plot.png", ATP5MG_plot,width = 6, height = 3.7, dpi =300,)




#STEP11: WGCNA
#OG vignette taken down, here is copy
#https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/RNA-SeqIntro/wgcna.html#gsc.tab=0
library(WGCNA)
options(stringsAsFactors = FALSE)
#get normalized results
liv_wg <- DESeq2::counts(dds_liv_all, normalized=TRUE)
#transpose dataframe
liv_wgt = as.data.frame(t(liv_wg))
#check to see if all genes and samples are okay
liv_gsg = goodSamplesGenes(liv_wgt, verbose = 3)
#if not okay from below run...
liv_gsg$allOK
#...me
if (!liv_gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!liv_gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(liv_wgt)[!liv_gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(liv_wgt)[!liv_gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  liv_wgt = liv_wgt[liv_gsg$goodSamples, liv_gsg$goodGenes]
}
#all good now?
liv_gsg = goodSamplesGenes(liv_wgt, verbose = 3)
liv_gsg$allOK
# check to see if there are any obvious outliers
liv_wgt_sampleTree = hclust(dist(liv_wgt), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(liv_wgt_sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#Remove outliers, each stage trait
#remove outliers
# Plot a line to show the cut
abline(h = 800000, col = "red");
# Determine cluster under the line
liv_wgt_clust = cutreeStatic(liv_wgt_sampleTree, cutHeight = 800000, minSize = 10)
liv_wgt_clust
#all good, but would run ... if not
# clust 1 contains the samples we want to keep.
liv_wgt_keepSamples = (liv_wgt_clust==1)
liv_wgt_datExpr = liv_wgt[liv_wgt_keepSamples, ]
liv_wgt_nGenes = ncol(liv_wgt_datExpr)
liv_wgt_nSamples = nrow(liv_wgt_datExpr) 
liv_wgt_nSamples
#Now put in traits
liv_trait<-as.data.frame(bodymass)
rownames(liv_trait) <- rownames(liv_wgt_datExpr)
#cluster
liv_wgt_sampleTree2 = hclust(dist(liv_wgt_datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(liv_trait, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(liv_wgt_sampleTree2, traitColors,
                    groupLabels = names(liv_trait),
                    main = "Sample dendrogram and trait heatmap")
save(liv_wgt_datExpr, liv_trait, file = "/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/WGCNA/BodyMass/Liv-01-dataInput.RData")

#WGCNA Part2
options(stringsAsFactors = FALSE);
enableWGCNAThreads()
# Load the data saved in the first part, or save as new name if conducting in same script
lnames = load(file = "/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/WGCNA/BodyMass/Liv-01-dataInput.RData");
#Note above changes for each set
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(liv_wgt_datExpr, powerVector = powers, verbose = 5)
#FOR 3 AND 4 ABOVE IT liv_wgt
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
##NOW NEED TO PICK POWER THRESHOLD
#7 for removed outliers as it crosses .8 threshold and according to https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html we are close
#link no longer works above, but see here https://www.biostars.org/p/9550786/
softPower = 5;
adjacency = adjacency(liv_wgt_datExpr, power = softPower);
#Creating TOM  Topological Overlap Matrix,
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
#Clustering using TOM
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
table(dynamicMods)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(liv_wgt_datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(liv_wgt_datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
#Plot the new merged set
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs
save(MEs, moduleLabels, moduleColors, geneTree, file = "/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/WGCNA/BodyMass/liv_1-02-networkConstruction-stepByStep.RData")

#WGCNA Part3
#load names again if running in different window
options(stringsAsFactors = FALSE)
lnames_01 = load(file = "/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/WGCNA/BodyMass//Liv-01-dataInput.RData")
nGenes = ncol(liv_wgt_datExpr)
nSamples = nrow(liv_wgt_datExpr);
MEs0 = moduleEigengenes(liv_wgt_datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, liv_trait, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(liv_trait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#Figure 3A
library(viridis)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(liv_trait),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = rev(blueWhiteRed(50)),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               bg.lab.y = NULL,
               main = paste("Module-trait relationships"))
viridis(50, alpha = 1, begin = 0, end = 1, direction = 1, option = "C")
#

bodymass2 = as.data.frame(liv_trait$bodymass)
names(bodymass2) = "bodymass2"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(liv_wgt_datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(liv_wgt_datExpr, bodymass2, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(bodymass2), sep="");
names(GSPvalue) = paste("p.GS.", names(bodymass2), sep="");
module = "grey60"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for bodymass2",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# Create the starting data frame
geneInfo0 = data.frame(Genes = names(liv_wgt_datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
geneInfo0$Genes
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, bodymass2, use = "p")));
modOrder
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.bodymass2));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/WGCNA/BodyMass//liv_geneInfo.csv")
# Select modules
modules = c("red");
# Select module probes
probes = names(liv_wgt_datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("~/data/Liver/WGCNA/CytoscapeInput-edges-Liv_grey_thresh", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("~/data/Liver/WGCNA/CytoscapeInput-nodes-Liv_grey_thresh", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);
filtered_df <- df[df$nodeName %in% gene_names, ]
liv_geneInfo <- read_csv("~/data/Liver/WGCNA/liv_geneInfo.csv")
sub_livgeneInfoS<-subset(liv_geneInfo,moduleColor== "red" & p.MM.red < 0.05)
sub_livgeneInfo<-subset(liv_geneInfo,moduleColor== "red")
sub_livgeneInfo
write.table(sub_livgeneInfoS, file='~/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/WGCNA/red_sig.tsv', quote=FALSE, sep='\t')
write.table(sub_livgeneInfo, file='~/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/WGCNA/red.tsv', quote=FALSE, sep='\t')





#Revision Likelihood ratio test
dds_lrt <- DESeq(dds_liv_all, test="LRT", reduced = ~ 1)
res_LRT <- results(dds_lrt)
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
write.table(res_LRT_tb, file='/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/Supplements/Data/Supplemental_Data_4.tsv', quote=FALSE, sep='\t')
# Subset to return genes with padj < 0.05
sigLRT_genes <- res_LRT_tb %>% 
  filter(padj < 0.01)
vst_mat <- assay(vst_dds_liv_all)
cluster_vst <- vst_mat[sigLRT_genes$gene, ]
library(DEGreport)
library(lasso2)
meta<-liv_stages_organ_frame <- data.frame(
  liv_full1 = liv_full1,
  row.names = colnames(cs_liv.count.tsv$counts)
)

dim(cluster_vst)
head(colnames(cluster_vst))
clusters <- degPatterns(cluster_vst, metadata = meta, time = "liv_full1",minc = 150)
write.table(clusters$df, file='/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/Supplements/Data/Supplemental_Data_5.tsv', quote=FALSE, sep='\t')
#figure
degplot <- degPlotCluster(clusters$normalized, "liv_full1", color = NULL, min_genes = 150,
               process = FALSE, points = TRUE, boxes = "FALSE", smooth = TRUE,
               lines = TRUE, facet = TRUE, cluster_column = "cluster")
degplot2<-degplot + 
  scale_color_manual(values = "black") +  # dots black
  stat_smooth(method = "loess", se = FALSE, color = "blue") +  # smooth in blue
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_blank(),         # remove axis titles
    axis.text = element_blank(),          # remove axis text
    axis.ticks = element_blank(),         # remove axis ticks
    legend.position = "none",             # remove legend
    strip.text = element_blank(),         # remove facet titles
    panel.grid.major = element_line(),    # show major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    plot.title = element_blank()
  )
ggsave("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/Supplements/Figures/degplot2.png", degplot2,width = 10, height = 5, dpi =300,)
#
cluster_groups <- clusters$df
group6 <- clusters$df %>%
  filter(cluster == 6)
group8 <- clusters$df %>%
  filter(cluster == 8)
group1 <- clusters$df %>%
  filter(cluster == 1)
group4 <- clusters$df %>%
  filter(cluster == 4)
#group overlap with wgcna clusters
#import WGCNA again
liv_geneInfo <- read_csv("ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/WGCNA/BodyMass/liv_geneInfo.csv")
#combine
all_groups <- bind_rows(
  group1,
  group4,
  group6,
  group8
)
all_groups <- group8

# Join with module information
combined_with_module <- all_groups %>%
  left_join(liv_geneInfo, by = c("genes" = "Genes")) %>%
  filter(!is.na(moduleColor))
combined_with_module <- all_groups %>%
  left_join(liv_geneInfo, by = c("genes" = "Genes"))
combined_with_module <- combined_with_module %>%
  filter(!is.na(moduleColor))
prop_df <- combined_with_module %>%
  group_by(cluster, moduleColor) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(proportion = count / sum(count))
# Step 2: Plot
global_order <- prop_df %>%
  arrange(desc(proportion)) %>%
  pull(moduleColor) %>%
  unique()
global_order <- rev(global_order)
#reorder
prop_df <- prop_df %>%
  mutate(moduleColor = factor(moduleColor, levels = global_order))

#plot with stacking
library(forcats)
clus8_stack<-ggplot(prop_df, aes(x = factor(cluster), y = proportion, fill = moduleColor)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_identity(guide = "legend") +
  labs(
    x = "Cluster",
    y = "Proportion of Genes",
    fill = "Module Color"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )
ggsave("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/Supplements/Figures/Supplementary_Figure_Transcriptome_clus8_stack.pdf", plot = clus8_stack, width = 2, height = 5)
#check
prop_df %>%
  filter(moduleColor == "red")


####
liv34resSig_df <- as.data.frame(liv34resSig) %>%
  mutate(GeneName = rownames(.))

# Join with module colors
sig_genes_with_module <- liv34resSig_df %>%
  left_join(liv_geneInfo, by = c("GeneName" = "Genes"))

# Count how many genes fall in each module
deg_per_module <- sig_genes_with_module %>%
  group_by(moduleColor) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

deg_per_module

########

#revision hamster comparison
library(readxl)
hiber_exp <- read_excel("ShrewProjects/Dehnels_Seasonal_RNAseq2/data/hiber_exp.xlsx")
hiber_exp <- hiber_exp %>%
  mutate(GeneName = toupper(GeneName))
hiber_up<-hiber_exp %>%
  filter(FDR < 0.05, logFC < 0) %>%
  pull(GeneName)
hiber_down<-hiber_exp %>%
  filter(FDR < 0.05, logFC > 0) %>%
  pull(GeneName)
#Highlight similar significance
liv13_tbl <- as_tibble(liv13res, rownames = "GeneName")
# Filter, join, and select all in one pipe
hiber_dehnel_df <- liv13_tbl %>%
  filter(padj < 0.05) %>%
  inner_join(hiber_exp %>% filter(FDR < 0.05), by = "GeneName") %>%
  select(GeneName, padj, log2FoldChange, FDR, logFC) %>%
  rename(
    Dehnel_padj = padj,
    Dehnel_LFC = log2FoldChange,
    Hibernation_FDR = FDR,
    Hibernation_LFC = logFC
  ) %>%
  # Multiply Hibernation_LFC by -1
  mutate(Hibernation_LFC = -1 * Hibernation_LFC) %>%
  # Create a group variable based on sign of LFCs
  mutate(
    LFC_sign_group = case_when(
      Dehnel_LFC > 0 & Hibernation_LFC > 0 ~ "++",
      Dehnel_LFC > 0 & Hibernation_LFC < 0 ~ "+-",
      Dehnel_LFC < 0 & Hibernation_LFC > 0 ~ "-+",
      Dehnel_LFC < 0 & Hibernation_LFC < 0 ~ "--",
      TRUE ~ "other"
    )
  ) %>%
  arrange(LFC_sign_group, Dehnel_padj) %>%
  select(-LFC_sign_group)
hiber_dehnel_df
write.table(hiber_dehnel_df, file='/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/Supplements/Tables/Supplemental_Table_6', quote=FALSE, sep='\t')
#
intersect(rownames(liv13up), hiber_up)
intersect(rownames(liv13down), hiber_up)
intersect(rownames(liv13up), hiber_down)
intersect(rownames(liv13down), hiber_down)
overlap_hiber <- list(HibernationUp = hiber_up, HibernationDown=hiber_down, DehnelUp=rownames(liv13up),DehnelDown=rownames(liv13down))
upset_overlap_hiber <- fromList(overlap_hiber)
hiber_upset <- upset(upset_overlap_hiber,nsets = 4,order.by = "freq")
ggsave("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Hibernation/upset_hibernation.png", hiber_upset, width = 10, height = 4, dpi =300,)
#fgsea
hiberres <- hiber_exp %>% 
  dplyr::select(GeneName, logFC) %>% 
  na.omit() %>% 
  distinct() %>% 
  mutate(logFC = -1 * logFC) %>% 
  group_by(GeneName) %>% 
  summarize(stat = mean(logFC))
hiberranks <- deframe(hiberres)
fgsea_hiber <- fgsea(pathways=gmtPathways("/Users/bill/ShrewProjects/github/Sorex_Genome2/data/000_miscFiles/c2.cp.kegg.v2023.1.Hs.symbols.gmt.txt"), hiberranks) %>% 
  as_tibble() %>% 
  arrange(padj)
fgsea_hiberTidy <- fgsea_hiber %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgsea_hiberTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()
plot_fgsea_hiber <-ggplot(subset(fgsea_hiberTidy,padj<0.05), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=NES>0)) +
  coord_flip() +
  labs(x="KEGGPathwayHibernation", y="Normalized Enrichment Score") + 
  scale_fill_manual(values =  c("red", "blue"))+
  theme_bw()
ggsave("/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Hibernation/fsgea_hibernation.png", plot_fgsea_hiber,width = 10, height = 3, dpi =300,)
fgsea_hiberTidy2 <- apply(fgsea_hiberTidy,2,as.character)
write.table(fgsea_hiberTidy2, file='/Users/bill/ShrewProjects/github/Dehnels_Seasonal_RNAseq2/data/Liver/DifferentialExp/Hibernation/fsgea_hibernation.tsv', quote=FALSE, sep='\t')





#Correlation between Gene expression and protein raw
match_trans<-DESeq2::counts(dds_liv_all, normalized=TRUE)
matched_LFQ <- read_table("~/Dehnels_Seasonal_RNAseq2/data/Liver/Proteomics/Compare2RNA/matched_LFQ_results.txt")
str(matched_LFQ)

#Prepare normalized RNA data
rna_df <- as.data.frame(match_trans) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "Sample", values_to = "RNA")

#Prepare protein data
prot_df <- matched_LFQ %>%
  pivot_longer(cols = starts_with("Stg"), names_to = "Sample", values_to = "Protein") %>%
  dplyr::select(Gene, Sample, Protein)

#Join RNA and Protein data by Gene + Sample
merged_df <- inner_join(rna_df, prot_df, by = c("Gene", "Sample"))
merged_df <- merged_df %>%
  mutate(logRNA = log10(RNA + 1),
         logProtein = log10(Protein))
#Calculate Pearson correlation
cor_results <- merged_df %>%
  group_by(Gene) %>%
  summarise(correlation = cor(logRNA, logProtein, use = "complete.obs", method = "pearson"),
            .groups = "drop")
#sort by correlation
cor_results <- cor_results %>%
  arrange(desc(correlation))

#total them
sum(cor_results$correlation > 0, na.rm = TRUE)
sum(cor_results$correlation < 0, na.rm = TRUE)

#plot
#hist
rna_prot_hist<-ggplot(cor_results, aes(x = correlation, fill = correlation > 0)) +
  geom_histogram(binwidth = 0.05, color = "white") +
  scale_fill_manual(values = c("FALSE" = "red", "TRUE" = "blue"), 
                    labels = c("Negative", "Positive"),
                    name = "Correlation Sign") +
  theme_bw() +
  labs( x = "Pearson r", y = "Gene Count")
ggsave("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/Proteomics/Compare2RNA/rna_prot_hist.png", rna_prot_hist,width = 5, height = 4, dpi =300,)

#scatter
summary(lm(logProtein ~ logRNA, data = merged_df))
rna_prot_scat<-ggplot(merged_df, aes(x = logRNA, y = logProtein)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(x = "log10 Normalized RNA Counts", y = "log10 LFQ Protein Abundance"
  ) +
  theme_bw(base_size = 14)
ggsave("/Users/bill/ShrewProjects/Dehnels_Seasonal_RNAseq2/data/Liver/Proteomics/Compare2RNA/rna_prot_scat.png", rna_prot_scat,width = 5, height = 4, dpi =300,)
