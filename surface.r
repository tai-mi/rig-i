# libraries ---------------------------------------------------------------

library(magrittr)
library(tidyverse)


# stuff -------------------------------------------------------------------

x <- readxl::read_xls(path='data/surface.xls',sheet='Sheet0') %>% 
  rename(well=`...1`)

# x <- x[,str_detect(names(x),'Freq. of Parent',negate=T)]
names(x) %<>% str_remove(' \\| Freq\\. of CD44\\+$') %>% str_remove('^cells/live/CD44\\+')
x <- column_to_rownames(x,'well') %>% select(-`/Ki67-`)
x %<>% mutate(across(.fns=~str_extract(.x,'.*(?= \\%)') %>% as.numeric()))

depth <- tibble(n=str_remove_all(colnames(x),'\\+') %>% 
                  str_remove('^/') %>% str_split('/'),
                lvl=str_count(colnames(x),'/'))
table(depth$lvl) # alright we good

dat1 <- tibble(well=str_extract(rownames(x),'^..'),
              c1=0,c2=0,c3=0,c4=0,c5=0,c6=0)
y <- x

#' for each unique lvl in depth
#' for each column at that lvl
#' if a given column is at a lower lvl and that column contains all of the markers in the higher lvl thing
#' then subtract the higher column from the lower column
for(l in rev(2:6)){
  for(i in which(depth$lvl==l)){
    for(i2 in 1:nrow(depth)){
      if((depth$lvl[i2]<l) & all(unlist(depth$n[i2]) %in% unlist(depth$n[i]))){
        x[i2] <- x[i2] - x[i]
      }
    }
  }
}
#'  well, I'm getting some negatives, but the most is -2.34, not sure why
#'  it isn't double counting i don't think?

for(i in 1:ncol(x)){
  dat1[depth$lvl[i]+1] <- dat1[depth$lvl[i]+1] + x[i]
}

dat1 %<>% mutate(c0=100-c1-c2-c3-c4-c5-c6,.after=well)
dat1 %>% pivot_longer(c(c0,c1,c2,c3,c4,c5,c6),names_to='n',values_to='val') %>% 
  rename(`Exhaustion markers expressed`=n) %>% 
  ggplot()+geom_col(aes(well,val,fill=`Exhaustion markers expressed`))


# stuff2 ------------------------------------------------------------------

x <- readxl::read_xls(path='data/surface.xls',sheet='Sheet0') %>% 
  set_names(c('depth','path','value','cells')) %>% 
  mutate(depth=str_count(depth,'\\>') %>% replace_na(0),
         value=case_when(str_detect(path,'fcs$')~100,
                         is.na(value)~0, 
                         T~as.numeric(value)))
