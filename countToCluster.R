require(Seurat)
require(dplyr)
require(jsonlite)
args <- commandArgs(trailingOnly = TRUE)

cur_sample_path <- args[1] #/*/*...tumor1.sample full path
out_dir <- args[2] #full path to outdir without / at the end
json_path <- args[3] # full path to json file
para<-read_json(json_path,simplifyVector = T)

res_cluster=para$res_cluster
max_pc=para$max_pc
top_feature= para$top_feature
nF_range= para$nF_range
pMit=para$pMit

#cur_sample <-basename(cur_sample_path)

#cur_sample_path<-'/gpfs/data01/glasslab/home/zhl022/testTumor/Tumor1.sample'
#out_dir <-"/gpfs/data01/glasslab/home/zhl022/testTumor"
#need mkdir ./intermediate


cur_sample <-basename(cur_sample_path)#Tumpr1.sample
sID <- gsub("\\.sample","",cur_sample) #Tumor1
strF <- paste(out_dir,"/Results/",sID,"/",sID,".rds",sep="")   #path for very first rds file


# if the first rds file does not exist, then make one
if(!file.exists(strF)){
  if(!dir.exists(dirname(strF))) {dir.create(dirname(strF)) }#create dir if not exist
  D <- read.table(cur_sample_path,sep="\t",header=T,as.is=T,row.names = 1) # read Tumor1.samplesheet, row name
  strD <- paste(out_dir,"/",sID,"/",rownames(D),sep="")     # read lines
  mergedSeurat <- NULL
  id1 <- NULL

  # to loop through every subsamples within sample such as microglia1,2,m3 under Tumor1
  for(j in strD){  # j is the rowname, microglia1, microglia2 etc
    data10X <- Read10X(paste(j,"/outs/filtered_feature_bc_matrix/",sep=""))  #read10X object
    data10X <- CreateSeuratObject(data10X,project=basename(j),min.cells = 3, min.features = 200)
    if(is.null(mergedSeurat)){
      mergedSeurat <- data10X
      id1 <- basename(j)
    }else{
      if(is.null(id1)){
        mergedSeurat <- merge(mergedSeurat,data10X,add.cell.id2=basename(j),do.normalize=F,project=sID)
      }else{
        mergedSeurat <- merge(mergedSeurat,data10X,add.cell.id1=id1,add.cell.id2=basename(j),do.normalize=F,project=sID)
        id1 <- NULL
      }
    }
  }
  saveRDS(mergedSeurat,file=strF)
  X<-mergedSeurat
  remove(mergedSeurat)

} else {
X<-readRDS(strF)    # load first rds if not exist
}


# compute %mt dna for qc
X[["percent.mt"]] <- PercentageFeatureSet(X, pattern = "^MT-")

# plot qc, three plots, nCount,nFeature,percent.mt, saved in intermediate folder
pdf(paste(out_dir,'/intermediate/',sID,'/qc_vinplot.pdf',sep=''),width=9)
VlnPlot(X, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(X, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(X, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()


# subset based on qc,need :
#nF_range, defalut is 0.95, then takes 0.025-0.975 nFeature
#pMit, percentile of mitochonria DNA, default 0.99
X <- subset(X, subset = nFeature_RNA > quantile(X@meta.data$nFeature_RNA,(1-nF_range)/2) & nFeature_RNA < quantile(X@meta.data$nFeature_RNA,1-(1-nF_range)/2) & percent.mt < quantile(X@meta.data$percent.mt,pMit))
X <- NormalizeData(X) # normalize data

# annotate variable feature
X <- FindVariableFeatures(X, selection.method = "vst", nfeatures  = top_feature,verbose = F)

#scale, just highly variable genes since later steps only use these Genes
#mean 0, var 1
allgenes <- rownames(X)
X <- ScaleData(X, features = allgenes)

#perform pca
X <- RunPCA(X, features = VariableFeatures(object = X),npcs=max_pc)
max_pc <- length(X@reductions$pca@stdev)

#plot pc and genes of early principle components

pdf(paste(out_dir,'/intermediate/',sID,'/pc_dim.pdf',sep=''),width=9)
VizDimLoadings(X, dims = 1:4, reduction = "pca")
DimHeatmap(X, dims = 1:9, cells = 500, balanced = TRUE)
dev.off()

#choose number of pc for clustering
#use jackstraw method, resample 100 times
#max_pc requires user input, default 100
X <- JackStraw(X, num.replicate = 100,dims=max_pc)
X <- ScoreJackStraw(X, dims = 1:max_pc)

#plot jackstraw and elbow
pdf(paste(out_dir,'/intermediate/',sID,'/pc_jack.pdf',sep=''),width=9)
JackStrawPlot(X, dims = 1:max_pc)
ElbowPlot(X,ndims = max_pc)
dev.off()

# compute number of pc to use for cluster
# pc must have p less than 0.001 from jackstraw
# require max_pc, if somewhere in max_pc greater than 0.001, n_pc is that stop point
#otherwise n_pc ==max_pc
if(min(which(X@reductions$pca@jackstraw@overall.p.values[,2]>0.001)==Inf)){
  n_pc=max_pc
} else {
  n_pc = min(which(X@reductions$pca@jackstraw@overall.p.values[,2]>0.001))
}

# cluster
# require para res_cluster, default=0.5
X <- FindNeighbors(X, dims = 1:n_pc)
X <- FindClusters(X, resolution = res_cluster)

#tsne
#RunTSNE on n_pc, save t_sne plot
X <- RunTSNE(X, dims = 1:n_pc)
pdf(paste(out_dir,'/Results/',sID,'/',sID,'_tsne.pdf',sep=''),width=9)
DimPlot(X, reduction = "tsne")
#FeaturePlot(X, features = c("APOE")) optional plot genes
dev.off()

# save cluster rds
saveRDS(X, file = paste(out_dir,"/Results/",sID,"/",sID,"_cluster.rds",sep="")


X.markers <- FindAllMarkers(X, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
#X.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

top_cluster_gene<- X.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
write.csv(top_cluster_gene, file = paste(out_dir,"/Results/",sID,"/",sID,"_topGenes.csv",sep=''))
