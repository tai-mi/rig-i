panglao <- function(genes,species1='human',tumor=F,nonadult=F){
  require(tidyverse)
  require(rvest)
  
  # assemble URL
  genes <- genes %>% str_remove_all('-') %>% str_remove('\\..*') %>% 
    paste0(collapse=',')
  species1 <- tolower(species1)
  if(species1 %in% c('h','human','hs','hu')) species1 <- 3
  if(species1 %in% c('m','mouse','mm','musmus','mu')) species1 <- 2
  if(species1 %in% c('any','both','hm','humu','all')) species1 <- 1
  tumor <- ifelse(tumor,1,0)
  nonadult <- ifelse(nonadult,1,0)
  url1 <- paste0('https://panglaodb.se/search.html?query="',genes,
                 '"&species=',species1,'&tumor=',tumor,'&nonadult=',nonadult)
  
  # get data
  hrefs <- read_html(url1) %>% html_elements('div') %>% 
    html_elements('div') %>% html_elements('a') %>% html_attr('href')
  url2 <- hrefs[str_detect(hrefs,'\\.csv$')]
  if(length(url2)!=1) break
  url2 <- paste0('https://panglaodb.se/',url2)
  
  # retrieve
  data.table::fread(url2) %>% 
    set_names(c('species','gene','source','index','tissue','type','rank')) %>% 
    select(-source,-index) %>% 
    # filter(type!='Unknown') %>% 
    mutate(type=case_when(type=='Gamma delta T cells'~'GD T cells',
        type=='Erythroid-like and erythroid precursor cells'~'Erythroid',
        type=='Plasmacytoid dendritic cells'~'pDC',
        type=='Dendritic cells'~'DC',
        type=='Luminal epithelial cells'~'Luminal epi cells',
        type=='Endothelial cells'~'EC',
        type=='Pulmonary alveolar type II cells'~'Alveolar type II',
        T ~ type))
}

panglaoMain <- function(...){
  data <- panglao(...)
  table(data$type) %>% sort(decreasing=T) %>% as.data.frame() %>% .[1:10,] %>% 
    ggplot()+geom_col(aes(Var1,Freq,fill=Var1))+
    theme(axis.text.x=element_text(angle=90),
          legend.position='none')
}
