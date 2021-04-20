# libraries ---------------------------------------------------------------

library(future)
library(magrittr)
library(tidyverse)
library(Seurat)

options(future.globals.maxSize = 15*1024^3)
plan('multisession')


# load --------------------------------------------------------------------

# files <- list.files('data/',full.names=T)
lower <- Read10X_h5("data/GSM4728858_SC174raw_feature_bc_matrix.h5")
upper <- Read10X_h5("data/GSM4728859_SC175raw_feature_bc_matrix.h5")

# lower <- lower[rowSums(lower)>10,]
# upper <- upper[rowSums(upper)>10,]

lower <- CreateSeuratObject(counts=lower,project='lower')
upper <- CreateSeuratObject(counts=upper,project='upper')

lower$propMt <- PercentageFeatureSet(lower,'^MT-')
upper$propMt <- PercentageFeatureSet(upper,'^MT-')

plotBasic <- function(data){
  a <- ggplot(data@meta.data)
  x <- a+geom_density(aes(nFeature_RNA),color='red')+scale_x_log10()+
    geom_vline(xintercept=300)
  y <- a+geom_density(aes(nCount_RNA),color='red')+scale_x_log10()+
    geom_vline(xintercept=500)
  z <- a+geom_density(aes(propMt),color='red')+scale_x_log10()+
    geom_vline(xintercept=20)
  cowplot::plot_grid(x,y,z)
}
plotBasic(lower)
plotBasic(upper)
lower <- subset(lower,subset=(nFeature_RNA>150)&(nCount_RNA>300)&(propMt<20))
upper <- subset(upper,subset=(nFeature_RNA>150)&(nCount_RNA>300)&(propMt<20))
# why are they just all cells with like 1 UMI, this data is whack


# other stuff -------------------------------------------------------------

lower %<>% SCTransform(vars.to.regress=c('nCount_RNA','nFeature_RNA','propMt'))
upper %<>% SCTransform(vars.to.regress=c('nCount_RNA','nFeature_RNA','propMt'))

ul <- c(upper,lower)
rm(upper)
rm(lower)
anchorFeats <- SelectIntegrationFeatures(ul)
ul <- PrepSCTIntegration(ul,anchor.features=anchorFeats)
anchorNames <- FindIntegrationAnchors(ul,anchor.features=anchorFeats,
                                      normalization.method='SCT')
rm(anchorFeats)
ul <- IntegrateData(anchorNames,normalization.method='SCT')
rm(anchorNames)
saveRDS(ul,'data/ul.rds')
ul <- readRDS('data/ul.rds')

ul %<>% RunPCA()
ElbowPlot(ul,ndims=50)
# ul %<>% RunUMAP(dims=1:40)
# UMAPPlot(ul)
ul %<>% FindNeighbors(dims=1:40)
ul %<>% FindClusters(resolution=0.7)
saveRDS(ul,'data/ul.rds')
ul %<>% RunUMAP(dims=1:40)
UMAPPlot(ul,label=T)

ul <- readRDS('data/human_lung.rds')

hpca <- celldex::HumanPrimaryCellAtlasData()
sr <- SingleR::SingleR(ul[['RNA']]@data,hpca,hpca$label.fine)

allMark <- readRDS('data/human_lung_allmark.rds')
allMark <- FindAllMarkers(ul)
topMark <- map(group_split(allMark,cluster),~arrange(.x,desc(avg_log2FC)) %>% .[1:15,]) %>% 
  set_names(group_by(allMark,cluster) %>% group_keys() %>% unlist()) %>% 
  map_dfr(~.x)

source('panglao.r')
panglao2 <- function(clust,n=10){
  allMark %>% filter(cluster==clust) %>% arrange(desc(avg_log2FC)) %>% 
    .$gene %>% .[1:n] %>% 
    panglaoMain()
}

explot <- function(genes){
  genes <- toupper(genes)
  ul$temp <- PercentageFeatureSet(ul,features=genes,assay='RNA')
  try(ul$temp <- PercentageFeatureSet(ul,features=genes,assay='SCT'),
      silent=T)
  FeaturePlot(ul,'temp',max.cutoff='q95',min.cutoff='q5')
}

