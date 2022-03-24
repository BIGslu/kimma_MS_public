library(tidyverse)
library(edgeR)
library(limma)
library(kimma)

#### Data ####
attach("data/RSTR_data_all.RData")

#Subset related and unrelated individuals
kin0 <- kin
diag(kin0) <- 0

#Look for related indiv
kin.summ <- as.data.frame(kin0) %>% 
  #Move rownames into df
  rownames_to_column() %>% 
  #Move column names into df
  pivot_longer(-rowname, names_to = "ptID") %>% 
  #Calculate max kinship value per person
  group_by(ptID) %>% 
  summarise(maxK = max(value)) %>% 
  #Create groups for individuals with and without related indiv in the data set
  mutate(group = ifelse(maxK >= 0.125, 'related',"unrelated"))

#Get FULLID of indiv to filter by
## Indiv with at least 1 related person in data set
related <- kin.summ %>% 
  filter(group == "related") %>% 
  pull(ptID)
## Indiv with no related persons in data set
unrelated <- kin.summ %>% 
  filter(group == "unrelated") %>% 
  pull(ptID)

#Subset unrelated
#Subset to MEDIA vs TB in LTBI samples
datV.media.unrelated <- datV
datV.media.unrelated$targets <- datV.media.unrelated$targets %>% 
  filter(ptID %in% unrelated & condition == "MEDIA")
datV.media.unrelated$E <- as.data.frame(datV.media.unrelated$E) %>% 
  select(all_of(datV.media.unrelated$targets$libID))
# identical(colnames(datV.media.unrelated$E), datV.media.unrelated$targets$libID)
datV.media.unrelated$weights <- as.data.frame(datV.media.unrelated$weights) %>% 
  select(all_of(datV.media.unrelated$targets$libID)) %>% 
  as.matrix()

kin.unrelated <- kin[unrelated,unrelated]

#Subset related
#Subset to MEDIA vs TB in LTBI samples
datV.media.related <- datV
datV.media.related$targets <- datV.media.related$targets %>% 
  filter(ptID %in% related & condition == "MEDIA")
datV.media.related$E <- as.data.frame(datV.media.related$E) %>% 
  select(all_of(datV.media.related$targets$libID))
# identical(colnames(datV.media.related$E), datV.media.related$targets$libID)
datV.media.related$weights <- as.data.frame(datV.media.related$weights) %>% 
  select(all_of(datV.media.related$targets$libID)) %>% 
  as.matrix()

kin.related <- kin[related,related]

#### unrelated ####
# Paired = N, kinship = N, weights = N
kimma_nnn_rstr_unrelated <- kmFit(dat = datV.media.unrelated, 
                                  model = "~ Sample_Group", run.lm = TRUE, 
                                  use.weights = FALSE)
kimma_nnn_rstr_unrelated$lm <- kimma_nnn_rstr_unrelated$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights",
         subset="unrelated")

# Paired = N, kinship = N, weights = Y
kimma_nny_rstr_unrelated <- kmFit(dat = datV.media.unrelated, 
                                  model = "~ Sample_Group", run.lm = TRUE, 
                                  use.weights = TRUE)
kimma_nny_rstr_unrelated$lm <- kimma_nny_rstr_unrelated$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "weights",
         subset="unrelated")

# Paired = N, kinship = Y, weights = N
kimma_nyn_rstr_unrelated <- kmFit(dat = datV.media.unrelated, kin=kin.unrelated,
                                       model = "~ Sample_Group + (1|ptID)", run.lmekin = TRUE, 
                                       use.weights = FALSE)
kimma_nyn_rstr_unrelated$lmekin <- kimma_nyn_rstr_unrelated$lmekin %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "no weights",
         subset="unrelated")

# Paired = Y, kinship = Y, weights = Y
kimma_nyy_rstr_unrelated <- kmFit(dat = datV.media.unrelated, kin=kin.unrelated,
                                       model = "~ Sample_Group + (1|ptID)", run.lmekin = TRUE, 
                                       use.weights = TRUE)
kimma_nyy_rstr_unrelated$lmekin <- kimma_nyy_rstr_unrelated$lmekin %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "weights",
         subset="unrelated")

#### related ####
# Paired = N, kinship = N, weights = N
kimma_nnn_rstr_related <- kmFit(dat = datV.media.related, 
                                  model = "~ Sample_Group", run.lm = TRUE, 
                                  use.weights = FALSE)
kimma_nnn_rstr_related$lm <- kimma_nnn_rstr_related$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights",
         subset="related")

# Paired = N, kinship = N, weights = Y
kimma_nny_rstr_related <- kmFit(dat = datV.media.related, 
                                  model = "~ Sample_Group", run.lm = TRUE, 
                                  use.weights = TRUE)
kimma_nny_rstr_related$lm <- kimma_nny_rstr_related$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "weights",
         subset="related")

# Paired = N, kinship = Y, weights = N
kimma_nyn_rstr_related <- kmFit(dat = datV.media.related, kin=kin.related,
                                     model = "~ Sample_Group + (1|ptID)", run.lmekin = TRUE, 
                                     use.weights = FALSE)
kimma_nyn_rstr_related$lmekin <- kimma_nyn_rstr_related$lmekin %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "no weights",
         subset="related")

# Paired = N, kinship = Y, weights = Y
kimma_nyy_rstr_related <- kmFit(dat = datV.media.related, kin=kin.related,
                                     model = "~ Sample_Group + (1|ptID)", run.lmekin = TRUE, 
                                     use.weights = TRUE)
kimma_nyy_rstr_related$lmekin <- kimma_nyy_rstr_related$lmekin %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "weights",
         subset="related")

#### Save ####

rstr.subset.ls <- list()

for(d in c("kimma_nnn_rstr_unrelated", "kimma_nny_rstr_unrelated",
           "kimma_nyn_rstr_unrelated", "kimma_nyy_rstr_unrelated",
           "kimma_nnn_rstr_related", "kimma_nny_rstr_related",
           "kimma_nyn_rstr_related", "kimma_nyy_rstr_related")){
  rstr.subset.ls[[d]] <- get(d)
}

save(rstr.subset.ls, file="results/rstr_subset_fit.RData")
