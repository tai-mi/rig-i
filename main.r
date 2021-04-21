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
  ggsave(paste0('figures/mlung_',title,'.png'),width=w,height=h)
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

dat <- readRDS('data/mouse_lung.rds')
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
