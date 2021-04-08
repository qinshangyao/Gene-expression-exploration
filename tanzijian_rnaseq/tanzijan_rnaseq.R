library(tidyverse)
library(ggsci)
library(ggthemes)
library(ggprism)

library(patchwork)

df <- read_delim("P2VSDEMO.All.txt",delim = "\t")
head(df)
names(df) <- c("gene","dpi_17","dpi_18","ac","d22","neuron","id")
##############
g1 <-c("ENSMUSG00000020052","ENSMUSG00000034701","ENSMUSG00000038255","ENSMUSG00000095139","ENSMUSG00000027967","ENSMUSG00000026686","ENSMUSG00000061911","ENSMUSG00000020950") #foxg1:ENSMUSG00000020950:df[23402,]
df[match(g1,df$id),]

g2 <- c("Dcx","Rbfox3","Syn2","Kcnip3","Pvalb","Map2")
df[match(g2,df$gene),]

g3 <- c("Sox2","Pax6","Fabp7")
df[match(g3,df$gene),]

g4 <- c("Gfap","S100b","Aqp4","Aldh1l1","Slc1a3","Slc1a2","Aldoc","Vim","Cd44","Nfia","Nfib","Sox9")

df[match(g4,df$gene),]
##################
##################

df[match(g1,df$id),]%>% dplyr::select( c("gene","ac","dpi_17","dpi_18","d22","neuron")) %>% 
  gather("group","expression",-gene) %>% 
  ggplot(aes(x= group,y= expression,fill = group)) + geom_col(width = 0.7) + facet_wrap(~ gene,scales = "free_y") + 
  scale_fill_nejm() + theme_few(base_size = 16)


df[match(g2,df$gene),]%>% dplyr::select( c("gene","ac","dpi_17","dpi_18","d22","neuron")) %>% 
  gather("group","expression",-gene) %>% 
  ggplot(aes(x= group,y= expression,fill = group)) + geom_col(width = 0.7) + facet_wrap(~ gene,scales = "free_y") + 
  scale_fill_nejm() + theme_few(base_size = 16)

df[match(g3,df$gene),]%>% dplyr::select( c("gene","ac","dpi_17","dpi_18","d22","neuron")) %>% 
  gather("group","expression",-gene) %>% 
  ggplot(aes(x= group,y= expression,fill = group)) + geom_col(width = 0.7) + facet_wrap(~ gene,scales = "free_y") + 
  scale_fill_nejm() + theme_few(base_size = 16)

df[match(g4,df$gene),]%>% dplyr::select( c("gene","ac","dpi_17","dpi_18","d22","neuron")) %>% 
  gather("group","expression",-gene) %>% 
  ggplot(aes(x= group,y= expression,fill = group)) + geom_col(width = 0.7) + facet_wrap(~ gene,scales = "free_y") + 
  scale_fill_nejm() + theme_few(base_size = 16)

###
df[match(g1,df$id),]%>% dplyr::select( c("gene","ac","dpi_17","dpi_18","d22","neuron")) %>% 
  gather("group","expression",-gene) %>% 
  ggplot(aes(x= group,y= expression,fill = group)) + geom_col(width = 0.7)  + 
  scale_fill_nejm() + theme_prism(palette = "floral" ,border = T) + facet_wrap(~ gene,scales = "free",nrow = 2) + 
  scale_y_continuous(guide = "prism_offset_minor") + scale_x_discrete(guide = "prism_bracket")+ 
  coord_cartesian(clip = "off")

df[match(g1,df$id),]%>% dplyr::select( c("gene","ac","dpi_17","dpi_18","d22","neuron")) %>% 
  gather("group","expression",-gene) %>% 
  ggplot(aes(x= group,y= expression,fill = group)) + geom_col(width = 0.7)  + 
  scale_fill_nejm() + theme_prism(palette = "floral" ,border = T) + facet_wrap(~ gene,scales = "free",nrow = 2) + 
  scale_y_continuous(guide = "prism_minor",expand = c(0,0)) + 
  coord_cartesian(clip = "off")

##############################################################
genes <- c("Sox3","Hes5")
df[df$gene_id %in% genes,] %>% select(gene_id,DEMO_1, DEMO_2, DEMO_3,  P2_1 , P2_2,  P2_3) -> df

df[df$gene_id == "Sox3",]$gene_id  <- "Sox2" 

df %>% 
  gather("group","expression",-gene_id) %>% 
  tidyr::separate(group,c("group","num")) %>% 
  ggplot(aes(x= group,y= expression,fill = group)) + stat_summary(geom = "col",fun = "mean",width = 0.7) +stat_summary(fun.data = mean_se, geom = "errorbar",width = 0.5)+ facet_wrap(~ fct_rev(gene_id),scales = "free_y") + 
  scale_fill_nejm() + theme_few(base_size = 16)


 ggplot(test2) + geom_roc(aes(d = D,m = group,color = M)) + ggsci::scale_color_lancet() + ggprism::theme_prism() + scale_x_continuous(expand = c(0,0)) + xlab("1- specificity") + geom_abline(intercept = 0,slope = 1,linetype = "dashed")

 ####################WB

