####Graphing metabolomics normalized data 08/01/2024
metab <- read_csv("~/data/Metaboanaylst_08_01_24/data_processed.csv")
#transpose the frame
metab <- as.data.frame(metab)
metab_future_cols <- metab$...1
rownames(metab) <- metab$...1
metab <- metab[, -1]
metab_trans <- as.data.frame(t(metab))
metab_trans$Stage <- c(rep("Stage1",5),rep("Stage2",4),rep("Stage3",5),rep("Stage4",5),rep("Stage5",5))
metab_trans$Sample <- rownames(metab_trans)



###Repeat below for each desired metabolite
ggplot(metab_trans,aes(x = Sample, y = `Arachidonic acid`, fill =Stage))+
  geom_bar(stat="identity")+
  theme_classic()
arach <- as.data.frame(cbind(metab_trans$Sample,metab_trans$Stage,metab_trans$`Arachidonic acid`))
arach$V3 <- as.numeric(as.character(arach$V3))
arach_df <- arach %>%
  group_by(V2) %>%
  summarise(
    mean_V3 = mean(V3, na.rm=TRUE),
    sd_V3 = sd(V3)
  )

# Plot the mean with error bars
arach_plot<-ggplot(arach_df, aes(x = V2, y = mean_V3, fill = V2)) +
  geom_bar(stat = "identity",show.legend = FALSE) +
  geom_errorbar(aes(ymin = mean_V3 - sd_V3, ymax = mean_V3 + sd_V3), width = 0.2) +
  scale_fill_manual(values =  alpha(c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2"),.80))+
  theme_classic()+
  scale_y_continuous(expand = expansion(mult = 0))#+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.title.y=element_blank())
ggsave("~/data/Metaboanaylst_08_01_24/Arach.png", arach_plot,width = 2.5, height = 3.09, dpi =300,)
