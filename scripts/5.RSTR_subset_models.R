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
kimma_nnn_rstr_unrelated$lm.fit <- kimma_nnn_rstr_unrelated$lm.fit %>% 
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
kimma_nny_rstr_unrelated$lm.fit <- kimma_nny_rstr_unrelated$lm.fit %>% 
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
kimma_nyn_rstr_unrelated$lmekin.fit <- kimma_nyn_rstr_unrelated$lmekin.fit %>% 
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
kimma_nyy_rstr_unrelated$lmekin.fit <- kimma_nyy_rstr_unrelated$lmekin.fit %>% 
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
kimma_nnn_rstr_related$lm.fit <- kimma_nnn_rstr_related$lm.fit %>% 
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
kimma_nny_rstr_related$lm.fit <- kimma_nny_rstr_related$lm.fit %>% 
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
kimma_nyn_rstr_related$lmekin.fit <- kimma_nyn_rstr_related$lmekin.fit %>% 
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
kimma_nyy_rstr_related$lmekin.fit <- kimma_nyy_rstr_related$lmekin.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "weights",
         subset="related")

#### Save ####
#Save indiv results for easier updating in future
save(kimma_nnn_rstr_unrelated,kimma_nny_rstr_unrelated,
     kimma_nyn_rstr_unrelated,kimma_nyy_rstr_unrelated,
     kimma_nnn_rstr_related,kimma_nny_rstr_related,
     kimma_nyn_rstr_related,kimma_nyy_rstr_related,
     file="results/model_fit/sample_group/Sample_Group_subset_all.RData")

#Combine and format 1 df
rstr_subset_result <- bind_rows(kimma_nnn_rstr_unrelated$lm,
                                kimma_nny_rstr_unrelated$lm,
                                kimma_nyn_rstr_unrelated$lmekin,
                                kimma_nyy_rstr_unrelated$lmekin,
                                kimma_nnn_rstr_related$lm,
                                kimma_nny_rstr_related$lm,
                                kimma_nyn_rstr_related$lmekin,
                                kimma_nyy_rstr_related$lmekin) %>% 
  mutate(variable = recode(variable, "Sample_GroupRSTR"="Sample_Group")) %>% 
  filter(variable == "Sample_Group") %>% 
  select(subset, software:weights,gene:FDR)

rownames(rstr_subset_result) <- NULL

#Combine fit into 1 df
rstr_subset_metric <- bind_rows(kimma_nnn_rstr_unrelated$lm.fit,
                                kimma_nny_rstr_unrelated$lm.fit,
                                kimma_nyn_rstr_unrelated$lmekin.fit,
                                kimma_nyy_rstr_unrelated$lmekin.fit,
                                kimma_nnn_rstr_related$lm.fit,
                                kimma_nny_rstr_related$lm.fit,
                                kimma_nyn_rstr_related$lmekin.fit,
                                kimma_nyy_rstr_related$lmekin.fit) %>% 
  select(subset, software:weights,gene:adj_Rsq)

#Save
save(rstr_subset_result, rstr_subset_metric, 
     file="results/rstr_subset_fit.RData")
