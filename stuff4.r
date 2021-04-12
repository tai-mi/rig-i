# libraries ---------------------------------------------------------------

library(future)
library(magrittr)
library(tidyverse)
library(Seurat)

options(future.globals.maxSize = 15*1024^3)
memory.limit(size=20000)
plan('multisession')


# load --------------------------------------------------------------------

dat <- Read10X('data/GSM4820777/')
dat <- CreateSeuratObject(counts=dat,project='project')

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

dat %<>% SCTransform(vars.to.regress=c('nCount_RNA','nFeature_RNA','propMt'))

dat %<>% RunPCA()
ElbowPlot(dat,ndims=50)
PCAPlot(dat)
# dat %<>% RunUMAP(dims=1:40)
# UMAPPlot(dat)
dat %<>% FindNeighbors(dims=1:30)
dat %<>% FindClusters(resolution=0.5)
dat %<>% RunUMAP(dims=1:10)
UMAPPlot(dat,label=T)
saveRDS(dat,'data/mouse_balf.rds')

hpca <- celldex::HumanPrimaryCellAtlasData()
sr <- SingleR::SingleR(dat[['RNA']]@data,hpca,hpca$label.fine)


allMark <- FindAllMarkers(dat)
topMark <- map(group_split(allMark,cluster),~arrange(.x,desc(avg_log2FC)) %>% .[1:15,]) %>% 
  set_names(group_by(allMark,cluster) %>% group_keys() %>% unlist()) %>% 
  map_dfr(~.x)


# cluster annot -----------------------------------------------------------

clustMap <- c(
  'macro/neutro',
  'endothelial/fibro',
  'endothelial',
  'macro/neutro',
  'fibro',
  't/tmem',
  'macro',
  'b',
  'fibro/ec',
  'ec',
  'alveolar',
  'macro?',
  'plasma',
  'macro?',
  'fibro',
  'nk?',
  'neutro',
  'neuron/ec/alveolar',
  't/fibro',
  'macro',
  'fibro',
  'macro?',
  'neuron?',
  'fibro/dividing',
  'histones?',
  25
)
