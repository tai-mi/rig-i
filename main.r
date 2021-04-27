# libraries ---------------------------------------------------------------

library(future)
library(magrittr)
library(tidyverse)
library(Seurat)

options(future.globals.maxSize = 1000*1024^3)
plan('multisession')


# load --------------------------------------------------------------------

gs <- readxl::read_xlsx('data/rig-i genes.xlsx',sheet='Sheet1') %>% 
  select(-protein) %>% 
  filter(gene!='IFNB3')
ggs <- function(title,w=6,h=4){
  ggsave(paste0('figures/hlung_',title,'.png'),width=w,height=h)
}
explot <- function(genes,x1='q95',x2='q1'){
  genes <- toupper(genes)
  # genes <- tolower(genes) %>% Hmisc::capitalize()
  dat$temp <- PercentageFeatureSet(dat,features=genes,assay='RNA')
  try(dat$temp <- PercentageFeatureSet(dat,features=genes,assay='SCT'),
      silent=T)
  FeaturePlot(dat,'temp',max.cutoff=x1,min.cutoff=x2)
}
sf <- function(genes){
  genes %<>% toupper()
  if(length(genes)==1){
    rownames(dat[['RNA']]@data)[str_detect(rownames(dat[['RNA']]@data),genes)]
  } else genes[genes %in% rownames(dat[['RNA']]@data)]
}


# cycle -------------------------------------------------------------------

dat <- readRDS('data/human_lung.rds')
for(g in gs$gene){
  mxcutoff <- 'q95'
  mncutoff <- 'q1'
  # g %<>% tolower() %>% Hmisc::capitalize()
  # while(T){
  tryCatch({
    FeaturePlot(dat,g,max.cutoff=mxcutoff,min.cutoff=mncutoff) %>% 
      print()
    ggs(g)
  },error=function(x){
             try({
               explot(g,mxcutoff,mncutoff) %>% print()
               ggs(g)
             },silent=T)
           })
    
    # h1 <- readline('')
    # if(h1=='') {ggs(g); break}
    # else{
    #   h1 %<>% str_split('-') %>% unlist()
    #   mxcutoff <- if_else(h1[1]=='',mxcutoff,paste0('q',h1[1]))
    #   mncutoff <- if_else(h1[2]=='',mncutoff,paste0('q',h1[2]))
    # }
  # }
}

# save scaled data for all and for groups
dat <- ScaleData(dat,features=gs$gene,assay='RNA',
                 vars.to.regress=c('nCount_RNA','nFeature_RNA','propMt'))
rn <- rownames(dat[['RNA']]@scale.data) %>% sort()
dat$temp <- colSums(dat[['RNA']]@scale.data)
FeaturePlot(dat,'temp',max.cutoff='q95',min.cutoff='q5')
ggs('scaled_all.png')
dat$temp <- colSums(dat[['RNA']]@scale.data[gs$gene[1:3][gs$gene[1:3] %in% rn],])
FeaturePlot(dat,'temp',max.cutoff='q95',min.cutoff='q5')
ggs('scaled_rna_sensors.png')
dat$temp <- colSums(dat[['RNA']]@scale.data[gs$gene[4:8][gs$gene[4:8] %in% rn],])
FeaturePlot(dat,'temp',max.cutoff='q95',min.cutoff='q5')
ggs('scaled_downstream_sensor_tfs.png')
dat$temp <- colSums(dat[['RNA']]@scale.data[gs$gene[9:24][gs$gene[9:24] %in% rn],])
FeaturePlot(dat,'temp',max.cutoff='q95',min.cutoff='q5')
ggs('scaled_IFN1.png')
dat$temp <- colSums(dat[['RNA']]@scale.data[gs$gene[25:26][gs$gene[25:26] %in% rn],])
FeaturePlot(dat,'temp',max.cutoff='q95',min.cutoff='q5')
ggs('scaled_IFN1R.png')
dat$temp <- colSums(dat[['RNA']]@scale.data[gs$gene[27:29][gs$gene[27:29] %in% rn],])
FeaturePlot(dat,'temp',max.cutoff='q95',min.cutoff='q5')
ggs('scaled_downstream_ifn_tfs.png')

