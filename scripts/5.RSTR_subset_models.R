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
#Subset to RSTR vs LTBI in TB samples
datV.tb.unrelated <- datV
datV.tb.unrelated$targets <- datV.tb.unrelated$targets %>% 
  filter(ptID %in% unrelated & condition == "TB")
datV.tb.unrelated$E <- as.data.frame(datV.tb.unrelated$E) %>% 
  select(all_of(datV.tb.unrelated$targets$libID))
# identical(colnames(datV.tb.unrelated$E), datV.tb.unrelated$targets$libID)
datV.tb.unrelated$weights <- as.data.frame(datV.tb.unrelated$weights) %>% 
  select(all_of(datV.tb.unrelated$targets$libID)) %>% 
  as.matrix()

kin.unrelated <- kin[unrelated,unrelated]

#Subset related
#Subset to RSTR vs LTBI in TB samples
datV.tb.related <- datV
datV.tb.related$targets <- datV.tb.related$targets %>% 
  filter(ptID %in% related & condition == "TB")
datV.tb.related$E <- as.data.frame(datV.tb.related$E) %>% 
  select(all_of(datV.tb.related$targets$libID))
# identical(colnames(datV.tb.related$E), datV.tb.related$targets$libID)
datV.tb.related$weights <- as.data.frame(datV.tb.related$weights) %>% 
  select(all_of(datV.tb.related$targets$libID)) %>% 
  as.matrix()

kin.related <- kin[related,related]

#### unrelated ####
# Paired = N, kinship = N, weights = N
kimma_nnn_rstr_unrelated <- kmFit(dat = datV.tb.unrelated, 
                                  model = "~ Sample_Group", run.lm = TRUE, 
                                  use.weights = FALSE, metrics = TRUE)
kimma_nnn_rstr_unrelated$lm <- kimma_nnn_rstr_unrelated$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights",
         subset="unrelated")
kimma_nnn_rstr_unrelated$lm.fit <- kimma_nnn_rstr_unrelated$lm.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights",
         subset="unrelated")

# Paired = N, kinship = N, weights = Y
kimma_nny_rstr_unrelated <- kmFit(dat = datV.tb.unrelated, 
                                  model = "~ Sample_Group", run.lm = TRUE, 
                                  use.weights = TRUE, metrics = TRUE)
kimma_nny_rstr_unrelated$lm <- kimma_nny_rstr_unrelated$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "voom weights",
         subset="unrelated")
kimma_nny_rstr_unrelated$lm.fit <- kimma_nny_rstr_unrelated$lm.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "voom weights",
         subset="unrelated")

# Paired = N, kinship = Y, weights = N
kimma_nyn_rstr_unrelated <- kmFit(dat = datV.tb.unrelated, kin=kin.unrelated,
                                       model = "~ Sample_Group + (1|ptID)", run.lmerel = TRUE, 
                                       use.weights = FALSE, metrics = TRUE)
kimma_nyn_rstr_unrelated$lmerel <- kimma_nyn_rstr_unrelated$lmerel %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "no weights",
         subset="unrelated")
kimma_nyn_rstr_unrelated$lmerel.fit <- kimma_nyn_rstr_unrelated$lmerel.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "no weights",
         subset="unrelated")

# Paired = Y, kinship = Y, weights = Y
kimma_nyy_rstr_unrelated <- kmFit(dat = datV.tb.unrelated, kin=kin.unrelated,
                                       model = "~ Sample_Group + (1|ptID)", run.lmerel = TRUE, 
                                       use.weights = TRUE, metrics = TRUE)
kimma_nyy_rstr_unrelated$lmerel <- kimma_nyy_rstr_unrelated$lmerel %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "voom weights",
         subset="unrelated")
kimma_nyy_rstr_unrelated$lmerel.fit <- kimma_nyy_rstr_unrelated$lmerel.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "voom weights",
         subset="unrelated")

#### related ####
# Paired = N, kinship = N, weights = N
kimma_nnn_rstr_related <- kmFit(dat = datV.tb.related, 
                                  model = "~ Sample_Group", run.lm = TRUE, 
                                  use.weights = FALSE, metrics = TRUE)
kimma_nnn_rstr_related$lm <- kimma_nnn_rstr_related$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights",
         subset="related")
kimma_nnn_rstr_related$lm.fit <- kimma_nnn_rstr_related$lm.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights",
         subset="related")

# Paired = N, kinship = N, weights = Y
kimma_nny_rstr_related <- kmFit(dat = datV.tb.related, 
                                  model = "~ Sample_Group", run.lm = TRUE, 
                                  use.weights = TRUE, metrics = TRUE)
kimma_nny_rstr_related$lm <- kimma_nny_rstr_related$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "voom weights",
         subset="related")
kimma_nny_rstr_related$lm.fit <- kimma_nny_rstr_related$lm.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "voom weights",
         subset="related")

# Paired = N, kinship = Y, weights = N
kimma_nyn_rstr_related <- kmFit(dat = datV.tb.related, kin=kin.related,
                                     model = "~ Sample_Group + (1|ptID)", run.lmerel = TRUE, 
                                     use.weights = FALSE, metrics = TRUE)
kimma_nyn_rstr_related$lmerel <- kimma_nyn_rstr_related$lmerel %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "no weights",
         subset="related")
kimma_nyn_rstr_related$lmerel.fit <- kimma_nyn_rstr_related$lmerel.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "no weights",
         subset="related")

# Paired = N, kinship = Y, weights = Y
kimma_nyy_rstr_related <- kmFit(dat = datV.tb.related, kin=kin.related,
                                     model = "~ Sample_Group + (1|ptID)", run.lmerel = TRUE, 
                                     use.weights = TRUE, metrics = TRUE)
kimma_nyy_rstr_related$lmerel <- kimma_nyy_rstr_related$lmerel %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "voom weights",
         subset="related")
kimma_nyy_rstr_related$lmerel.fit <- kimma_nyy_rstr_related$lmerel.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "voom weights",
         subset="related")

#### Save ####
#Save indiv results for easier updating in future
save(kimma_nnn_rstr_unrelated,kimma_nny_rstr_unrelated,
     kimma_nyn_rstr_unrelated,kimma_nyy_rstr_unrelated,
     kimma_nnn_rstr_related,kimma_nny_rstr_related,
     kimma_nyn_rstr_related,kimma_nyy_rstr_related,
     file="results/model_fit/sample_group/Sample_Group_subset_all_tb.RData")

#Combine and format 1 df
rstr_subset_result <- bind_rows(kimma_nnn_rstr_unrelated$lm,
                                kimma_nny_rstr_unrelated$lm,
                                kimma_nyn_rstr_unrelated$lmerel,
                                kimma_nyy_rstr_unrelated$lmerel,
                                kimma_nnn_rstr_related$lm,
                                kimma_nny_rstr_related$lm,
                                kimma_nyn_rstr_related$lmerel,
                                kimma_nyy_rstr_related$lmerel) %>% 
  mutate(variable = recode(variable, "Sample_GroupRSTR"="Sample_Group")) %>% 
  filter(variable == "Sample_Group") %>% 
  select(subset, software:weights,gene:FDR)

rownames(rstr_subset_result) <- NULL

#Combine fit into 1 df
rstr_subset_metric <- bind_rows(kimma_nnn_rstr_unrelated$lm.fit,
                                kimma_nny_rstr_unrelated$lm.fit,
                                kimma_nyn_rstr_unrelated$lmerel.fit,
                                kimma_nyy_rstr_unrelated$lmerel.fit,
                                kimma_nnn_rstr_related$lm.fit,
                                kimma_nny_rstr_related$lm.fit,
                                kimma_nyn_rstr_related$lmerel.fit,
                                kimma_nyy_rstr_related$lmerel.fit) %>% 
  select(subset, software:weights,gene:adj_Rsq)

#Save
save(rstr_subset_result, rstr_subset_metric, 
     file="results/rstr_tb_subset_fit.RData")
