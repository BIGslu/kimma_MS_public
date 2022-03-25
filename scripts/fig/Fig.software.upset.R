library(tidyverse)
library(UpSetR)

#### Data ####
load("results/condition_fit.RData")
load("results/rstr_fit.RData")

load("results/condition_subset_fit.RData")
load("results/rstr_fit.RData")

#### Condition upset ####
fdr.cutoff <- 0.05
condition.ls <- list()

for(s in unique(condition_result$software)){
  for(p in unique(condition_result$paired)){
    for(k in unique(condition_result$kinship)){
      for(w in unique(condition_result$weights)){
        nme <- paste(k,p,w,s, sep="_")
        temp <- condition_result %>% 
          filter(software == s & paired == p & kinship == k &
                   weights == w & FDR < fdr.cutoff) %>% 
          distinct(gene) %>% 
          pull(gene)
        
        if(length(temp > 0)){ condition.ls[[nme]] <- temp }
      }}}}

#Make metadata for sets
condition.meta <- data.frame(set = names(condition.ls)) %>% 
  separate(set, into=c("kinship","paired","weights","software"),
                       sep="_", remove=FALSE)
#No kinship
noK <- sort(names(condition.ls)[grepl("no kinship", names(condition.ls))])
condition.meta.noK <- condition.meta %>% filter(set %in% noK)

upset(fromList(condition.ls[noK]),
      nsets = 12, keep.order=TRUE, sets=noK, order.by="freq",
      set.metadata = list(data = condition.meta.noK, 
                          plots = list(list(type = "heat", 
                                            column = "kinship", assign = 3, 
                                            colors = c(`no kinship` = "grey", 
                                                       kinship = "black")),
                                       list(type = "heat", 
                                            column = "paired", assign = 3, 
                                            colors = c(unpaired = "grey", 
                                                       paired = "black")),
                                       list(type = "heat", 
                                            column = "weights", assign = 3, 
                                            colors = c(`no weights` = "grey", 
                                                       weights = "black")),
                                       list(type = "heat", 
                                            column = "software", assign = 3, 
                                            colors = c(DESeq2 = "#fc8d62", limma = "#66c2a5",
                                                       kimma = "#e78ac3", dream="#8da0cb")))),
      queries = list(list(query = intersects,
                          params = list(noK[-c(1,8)]), color="red", active=TRUE),
                     list(query = intersects,
                          params = list(noK[-c(1)]), color="red", active=TRUE),
                     list(query = intersects,
                          params = list(noK[-c(8)]), color="red", active=TRUE),
                     list(query = intersects,
                          params = list(noK), color="red", active=TRUE),
                     list(query = intersects,
                          params = list(noK[c(5,7)]), color="green", active=TRUE)))

#With kinship
wK <- sort(names(condition.ls)[!grepl("no kinship", names(condition.ls))])

upset(fromList(condition.ls[wK]),
      nsets = 12, keep.order=TRUE, sets=wK, order.by="freq")

#### RSTR venn ####
fdr.cutoff <- 0.3
rstr.ls <- list()

for(s in unique(rstr_result$software)){
  for(p in unique(rstr_result$paired)){
    for(k in unique(rstr_result$kinship)){
      for(w in unique(rstr_result$weights)){
        nme <- paste(k,p,w,s, sep="_")
        temp <- rstr_result %>% 
          filter(software == s & paired == p & kinship == k &
                   weights == w & FDR < fdr.cutoff) %>% 
          distinct(gene) %>% 
          pull(gene)
        
        if(length(temp > 0)){ rstr.ls[[nme]] <- temp }
      }}}}

noK2 <- sort(names(rstr.ls)[grepl("no kinship", names(rstr.ls))])

upset(fromList(rstr.ls[noK2]),
      keep.order=TRUE, sets=noK2, order.by="freq")

wK2 <- sort(names(rstr.ls)[!grepl("no kinship", names(rstr.ls))])

upset(fromList(rstr.ls[wK2]),
      keep.order=TRUE, sets=wK2, order.by="freq")
