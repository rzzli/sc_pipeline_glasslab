library(Seurat)
setwd("/data/archive/20-01-21-Stevens-mouse-microglia/counts/dge_txt")

require("biomaRt")

replace_fux <- function(gene){
  return(as.character(new_v2[gene,]))
}


#retrive all mouse ids https://www.biostars.org/p/147351/
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
annot<-as.character(as.matrix(getBM(c( "mgi_symbol"), mart=ensembl)))

convertMouseGeneList <- function(x){
  
  human_ensembl <<-useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse_ensembl <<- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 <<- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse_ensembl, attributesL = c("hgnc_symbol"), martL = human_ensembl, uniqueRows=T)
  #humanx <<- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  return(genesV2)
}

## set two column table
homolog_V2<-convertMouseGeneList(annot)    # col1: mouse col2:human

x<-homolog_V2[!duplicated(homolog_V2$MGI.symbol),] #unique mouse
x1 <- x[!duplicated(x$HGNC.symbol),] # human and mouse unique

new_v2 <- data.frame(x1[,-1], row.names=x1[,1])
colnames(new_v2) <- "human"


mouse_file <- list.files()
merged_mouse<- NULL

for (sample in mouse_file){
  
  
  loaddf <<- read.table(sample,header = T,row.names = 1)
  sample_name <<- strsplit(sample,split = '[.]')[[1]][1]
  list_sample_name <<- strsplit(sample_name,"_")[[1]]
  name_to_use <<- paste(list_sample_name[-1],collapse = '_')
  
  rnames_mouse <<- rownames(loaddf)
  new_rnames<<-as.data.frame(sapply(rnames_mouse,FUN = replace_fux))
  
  temp_df <- cbind(new_rnames[,1],loaddf)
  
  new_df <- temp_df[!is.na(temp_df[,1]),]
  rownames(new_df) <- NULL
  row.names(new_df) <- make.names(new_df[,1],TRUE)
  new_df <- new_df[,-1]
  temp_obj <- CreateSeuratObject(new_df,project = name_to_use,min.cells = 3, min.features = 200)
  temp_obj<-RenameCells(object = temp_obj, add.cell.id = name_to_use)
  if (is.null(merged_mouse)) {
    merged_mouse<-temp_obj
    rm(temp_obj)
  } else {
    previous_merged<- merged_mouse
    merged_mouse <- merge(previous_merged,temp_obj)
    rm(temp_obj)
    rm(previous_merged)
  }
}
