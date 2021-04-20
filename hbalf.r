# libraries ---------------------------------------------------------------

library(future)
library(magrittr)
library(tidyverse)
library(Seurat)

options(future.globals.maxSize = 1000*1024^3)
plan('multisession')


# load --------------------------------------------------------------------

dat <- readRDS('data/Allcells.counts.rds')
meta <- data.table::fread('data/Allcells.meta.data.csv')

dat <- dat[,meta$Cell[meta$Disease=='control']] #only control samples
meta %<>% filter(Disease=='control') %>% column_to_rownames('Cell')
dat <- dat[rowSums(dat)>10,] #filter features
dim(dat)

dat <- CreateSeuratObject(counts=dat,project='project',meta.data=meta)

dat$propMt <- PercentageFeatureSet(dat,'^MT-')

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
plotBasic(dat)

dat <- subset(dat,subset=(nFeature_RNA>150)&(nCount_RNA>300)&(propMt<20))


# other stuff -------------------------------------------------------------

dat <- SplitObject(dat,'PatientNumber')
dat %<>% map(~SCTransform(.x,vars.to.regress=c('nCount_RNA',
                                               'nFeature_RNA','propMt')))

anchorFeats <- SelectIntegrationFeatures(dat)
dat <- PrepSCTIntegration(dat,anchor.features=anchorFeats)
anchorNames <- FindIntegrationAnchors(dat,anchor.features=anchorFeats,
                                      normalization.method='SCT')
rm(anchorFeats)
dat <- IntegrateData(anchorNames,normalization.method='SCT')
rm(anchorNames)

dat %<>% RunPCA()
ElbowPlot(dat,ndims=50)
# dat %<>% RunUMAP(dims=1:40)
# UMAPPlot(dat)
dat %<>% FindNeighbors(dims=1:30)
dat %<>% FindClusters(resolution=0.5)
dat %<>% RunUMAP(dims=1:50)
UMAPPlot(dat,label=T)
saveRDS(dat,'data/human_balf.rds')
dat <- readRDS('data/human_balf.rds')

hpca <- celldex::HumanPrimaryCellAtlasData()
sr <- SingleR::SingleR(dat[['RNA']]@data,hpca,hpca$label.fine)

allMark <- FindAllMarkers(dat)
topMark <- map(group_split(allMark,cluster),~arrange(.x,desc(avg_log2FC)) %>% .[1:15,]) %>% 
  set_names(group_by(allMark,cluster) %>% group_keys() %>% unlist()) %>% 
  map_dfr(~.x)

allMark <- readRDS('data/human_balf_allmark.rds')
source('panglao.r')
panglao2 <- function(clust){
  allMark %>% filter(cluster==clust) %>% arrange(desc(avg_log2FC)) %>% 
    .$gene %>% .[1:10] %>% 
    panglaoMain()
}

explot <- function(genes){
  genes <- toupper(genes)
  dat$temp <- PercentageFeatureSet(dat,features=genes,assay='RNA')
  try(dat$temp <- PercentageFeatureSet(dat,features=genes,assay='SCT'),
      silent=T)
  FeaturePlot(dat,'temp',max.cutoff='q95',min.cutoff='q5')
}


# cluster annot -----------------------------------------------------------

DimPlot(dat,group.by='CellType')
DotPlot(dat,assay='RNA',features='CD69',group.by='CellType')