library(GOfuncR)
library("vctrs")
library(dplyr)

mg <- read.csv("~/daima/python/enhancer_network_May_2020/microglia_edge_df.csv")
neu <- read.csv("~/daima/python/enhancer_network_May_2020/other_tissue/Neuron/neuron_edge_df.csv")
oli <- read.csv("~/daima/python/enhancer_network_May_2020/other_tissue/Oligo/oligo_edge_df.csv")



enrich <- function(df,numb){
  tf<- names(sort(table(df$Motif.Name),decreasing = T)[numb])
  gene_ids<-df[df$Motif.Name==tf,]$Gene.Name
  input_hyper = data.frame(gene_ids, is_candidate=1)
  res_hyper = go_enrich(input_hyper, n_randset=100)
  return(res_hyper)
}

res<-enrich(mg,1)

 
  
plotGO<-function(res){
  library(tidyverse)
  
  top<- res$results[1:10,]
  top$logp<- -log(top$raw_p_overrep)
  
  p<-top %>% 
    ggplot(aes(reorder(node_name,logp ), logp)) + 
    geom_col(aes(fill = logp)) +
    scale_fill_gradient2(low = "white", 
                         high = "blue", 
                         midpoint = median(top$logp)) + coord_cartesian(ylim = c(min(top$logp)-5,max(top$logp)+5))+ 
                          xlab("GO term") +ylab("logp")+ coord_flip() 
  return(p)
}  

plotGO(res)
