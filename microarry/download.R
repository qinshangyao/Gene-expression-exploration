## git@github.com:qinshangyao/Reactive_astrocyte.git
### 1. download the expression data from GEO

BiocManager::install("GEOquery")
library(GEOquery)
library(tidyverse)
install.packages("ggprism")
install.packages("ggsci")
install.packages("patchwork")
library(ggsci)
library(ggprism)
library(patchwork)
vignette("GEOquery")
rm(list = ls())

##1.download the data 
options( 'download.file.method.GEOquery' = 'libcurl' ) ## enable the network download
gse <- getGEO("GSE35338",getGPL = T,AnnotGPL = T)
eset <- gse[[1]]
rm(gse)


## 2.annotion data 
anno_data <- eset@featureData@data[,c("ID","Gene symbol")]
head(anno_data) 

#######GPL1261
BiocManager::install("mouse430a2.db")
library(mouse430a2.db)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
vignette("clusterProfiler")

keytypes(mouse430a2.db)
columns(mouse430a2.db)
#############mouse430a2SYMBOL %>% toTable()
bitr(anno_data$ID,fromType = "PROBEID",toType = "SYMBOL",OrgDb = "mouse430a2.db" ) -> anno_data2

### the experiment annotation
experiment_df <- pData(eset) %>% select(c("title","geo_accession")) 
experiment_df %>% separate(col = "title",into = c("cell","treatment","batch"),sep = ",") %>% select(2,4) -> experiment_df
experiment_df %>% filter(str_detect(treatment,"7")) -> experiment_df

##filter
eset <- eset[,sampleNames(eset) %in% experiment_df$geo_accession]
eset

## 3.expression data  
exprs_df <- exprs(eset)
dim(exprs_df)
exprs_df[1:3,1:6]
## 4. check the gene expression

tf_gene <- c("Neurog2","Ascl1","Pou3f2","Foxg1","Neurod1","Neurod2","Dlx2","Lmx1a","Mmgt1","Lhx1","Nr4a2")

stem_gene <- c("Sox2","Klf4","Pax6","Nanog","Utf1","Fbxo15","Cript","Fabp7","Tnc")


dat=tibble(
  prob = anno_data$ID[anno_data$`Gene symbol` %in% stem_gene],
  gene = anno_data$`Gene symbol`[anno_data$ID %in% prob],
  as_tibble(exprs_df[prob,])
)

## 5. plot boxplot

dat <- dat %>% gather(3:8,key = "GSM",value = "Exp") %>% mutate(group = if_else(GSM %in% c("GSM866334","GSM866335","GSM866336"),"reactive","normal"))
dat %>% group_by(gene,GSM,group) %>% summarise(Exp = median(Exp)) -> dat

dat %>%
  ggplot(aes(x = gene,y = log2(Exp),fill = group))  + stat_summary(geom= "col",fun = mean,width = 0.65,position = position_dodge(0.65)) + 
  stat_summary(geom = "errorbar",width=0.35,fun.data  = mean_se,position = position_dodge(0.65)) +
  geom_point(size = 2,color = "grey",position = position_dodge(0.65)) + 
  ggprism::theme_prism(base_size = 16,axis_text_angle = 0) + 
  theme(legend.position = c(0.9,0.2),axis.title.x = element_blank())+
  ggsci::scale_fill_d3() +
  scale_y_continuous(expand = c(0,0),guide = guide_prism_minor()) + ylab("log2 Expresion")   + scale_x_discrete(guide = "prism_bracket")+ 
  facet_wrap(~ gene,scales = "free",nrow = 3) + guides(fill = "none")
  





