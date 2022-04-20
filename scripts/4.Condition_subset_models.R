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
datV.ltbi.unrelated <- datV
datV.ltbi.unrelated$targets <- datV.ltbi.unrelated$targets %>% 
  filter(ptID %in% unrelated & Sample_Group == "LTBI")
datV.ltbi.unrelated$E <- as.data.frame(datV.ltbi.unrelated$E) %>% 
  select(all_of(datV.ltbi.unrelated$targets$libID))
# identical(colnames(datV.ltbi.unrelated$E), datV.ltbi.unrelated$targets$libID)
datV.ltbi.unrelated$weights <- as.data.frame(datV.ltbi.unrelated$weights) %>% 
  select(all_of(datV.ltbi.unrelated$targets$libID)) %>% 
  as.matrix()

kin.unrelated <- kin[unrelated,unrelated]

#Subset related
#Subset to MEDIA vs TB in LTBI samples
datV.ltbi.related <- datV
datV.ltbi.related$targets <- datV.ltbi.related$targets %>% 
  filter(ptID %in% related & Sample_Group == "LTBI")
datV.ltbi.related$E <- as.data.frame(datV.ltbi.related$E) %>% 
  select(all_of(datV.ltbi.related$targets$libID))
# identical(colnames(datV.ltbi.related$E), datV.ltbi.related$targets$libID)
datV.ltbi.related$weights <- as.data.frame(datV.ltbi.related$weights) %>% 
  select(all_of(datV.ltbi.related$targets$libID)) %>% 
  as.matrix()

kin.related <- kin[related,related]

#### unrelated ####
# Paired = Y, kinship = N, weights = N
kimma_ynn_condition_unrelated <- kmFit(dat = datV.ltbi.unrelated, 
                                       model = "~ condition + (1|ptID)", run.lme = TRUE, 
                                       use.weights = FALSE, metrics = TRUE)
kimma_ynn_condition_unrelated$lme <- kimma_ynn_condition_unrelated$lme %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "no kinship", weights = "no weights",
         subset = "unrelated")
kimma_ynn_condition_unrelated$lme.fit <- kimma_ynn_condition_unrelated$lme.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "no kinship", weights = "no weights",
         subset = "unrelated")

# Paired = Y, kinship = N, weights = Y
kimma_yny_condition_unrelated <- kmFit(dat = datV.ltbi.unrelated, 
                                       model = "~ condition + (1|ptID)", run.lme = TRUE, 
                                       use.weights = TRUE, metrics = TRUE)
kimma_yny_condition_unrelated$lme <- kimma_yny_condition_unrelated$lme %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "no kinship", weights = "weights",
         subset="unrelated")
kimma_yny_condition_unrelated$lme.fit <- kimma_yny_condition_unrelated$lme.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "no kinship", weights = "weights",
         subset="unrelated")

# Paired = Y, kinship = Y, weights = N
kimma_yyn_condition_unrelated <- kmFit(dat = datV.ltbi.unrelated, kin=kin.unrelated,
                             model = "~ condition + (1|ptID)", run.lmerel = TRUE, 
                             use.weights = FALSE, metrics = TRUE)
kimma_yyn_condition_unrelated$lmerel <- kimma_yyn_condition_unrelated$lmerel %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "no weights",
         subset = "unrelated")
kimma_yyn_condition_unrelated$lmerel.fit <- kimma_yyn_condition_unrelated$lmerel.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "no weights",
         subset = "unrelated")

# Paired = Y, kinship = Y, weights = Y
kimma_yyy_condition_unrelated <- kmFit(dat = datV.ltbi.unrelated, kin=kin.unrelated,
                             model = "~ condition + (1|ptID)", run.lmerel = TRUE, 
                             use.weights = TRUE, metrics = TRUE)
kimma_yyy_condition_unrelated$lmerel <- kimma_yyy_condition_unrelated$lmerel %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "weights",
         subset="unrelated")
kimma_yyy_condition_unrelated$lmerel.fit <- kimma_yyy_condition_unrelated$lmerel.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "weights",
         subset="unrelated")

#### related ####
# Paired = Y, kinship = N, weights = N
kimma_ynn_condition_related <- kmFit(dat = datV.ltbi.related, 
                                     model = "~ condition + (1|ptID)", run.lme = TRUE, 
                                     use.weights = FALSE, metrics = TRUE)
kimma_ynn_condition_related$lme <- kimma_ynn_condition_related$lme %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "no kinship", weights = "no weights",
         subset="related")
kimma_ynn_condition_related$lme.fit <- kimma_ynn_condition_related$lme.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "no kinship", weights = "no weights",
         subset="related")

# Paired = Y, kinship = N, weights = Y
kimma_yny_condition_related <- kmFit(dat = datV.ltbi.related, 
                                     model = "~ condition + (1|ptID)", run.lme = TRUE, 
                                     use.weights = TRUE, metrics = TRUE)
kimma_yny_condition_related$lme <- kimma_yny_condition_related$lme %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "no kinship", weights = "weights",
         subset="related")
kimma_yny_condition_related$lme.fit <- kimma_yny_condition_related$lme.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "no kinship", weights = "weights",
         subset="related")

# Paired = Y, kinship = Y, weights = N
kimma_yyn_condition_related <- kmFit(dat = datV.ltbi.related, kin=kin.related,
                                       model = "~ condition + (1|ptID)", run.lmerel = TRUE, 
                                       use.weights = FALSE, metrics = TRUE)
kimma_yyn_condition_related$lmerel <- kimma_yyn_condition_related$lmerel %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "no weights",
         subset="related")
kimma_yyn_condition_related$lmerel.fit <- kimma_yyn_condition_related$lmerel.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "no weights",
         subset="related")

# Paired = Y, kinship = Y, weights = Y
kimma_yyy_condition_related <- kmFit(dat = datV.ltbi.related, kin=kin.related,
                                       model = "~ condition + (1|ptID)", run.lmerel = TRUE, 
                                       use.weights = TRUE, metrics = TRUE)
kimma_yyy_condition_related$lmerel <- kimma_yyy_condition_related$lmerel %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "weights",
         subset="related")
kimma_yyy_condition_related$lmerel.fit <- kimma_yyy_condition_related$lmerel.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "weights",
         subset="related")

#### Save ####
#Save indiv results for easier updating in future
save(kimma_ynn_condition_unrelated, kimma_yny_condition_unrelated,
     kimma_yyn_condition_unrelated, kimma_yyy_condition_unrelated,
     kimma_ynn_condition_related, kimma_yny_condition_related,
     kimma_yyn_condition_related, kimma_yyy_condition_related,
     file="results/model_fit/condition/condition_subset_all.RData")

#Combine and format 1 df
condition_subset_result <- bind_rows(kimma_ynn_condition_unrelated$lme,
                                     kimma_yny_condition_unrelated$lme,
                                     kimma_yyn_condition_unrelated$lmerel,
                                     kimma_yyy_condition_unrelated$lmerel,
                                     kimma_ynn_condition_related$lme,
                                     kimma_yny_condition_related$lme,
                                     kimma_yyn_condition_related$lmerel,
                                     kimma_yyy_condition_related$lmerel) %>% 
  mutate(variable = recode(variable, "conditionTB"="condition")) %>% 
  filter(variable == "condition") %>% 
  select(subset, software:weights,gene:FDR)

rownames(condition_subset_result) <- NULL

#Combine fit into 1 df
condition_subset_metric <- bind_rows(kimma_ynn_condition_unrelated$lme.fit,
                                    kimma_yny_condition_unrelated$lme.fit,
                                    kimma_yyn_condition_unrelated$lmerel.fit,
                                    kimma_yyy_condition_unrelated$lmerel.fit,
                                    kimma_ynn_condition_related$lme.fit,
                                    kimma_yny_condition_related$lme.fit,
                                    kimma_yyn_condition_related$lmerel.fit,
                                    kimma_yyy_condition_related$lmerel.fit) %>% 
  select(subset, software:weights,gene:adj_Rsq)

#Save
save(condition_subset_result, condition_subset_metric,
     file="results/condition_subset_fit.RData")
