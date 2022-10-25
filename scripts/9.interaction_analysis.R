library(tidyverse)
library(kimma)

#### Data ####
attach("data/RSTR_RNAseq_data_combined_clean.RData")
dat.combined.voom <- dat.combined.voom
dat.combined.voom$weights <- as.data.frame(dat.combined.voom$weights)
rownames(dat.combined.voom$weights) <- rownames(dat.combined.voom$E)
colnames(dat.combined.voom$weights) <- colnames(dat.combined.voom$E)

kin <- load("data/RSTR_data_all.RData")

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
  filter(group == "related" & ptID %in% dat.combined.voom$targets$FULLIDNO) %>% 
  pull(ptID)
## Indiv with no related persons in data set
unrelated <- kin.summ %>% 
  filter(group == "unrelated" & ptID %in% dat.combined.voom$targets$FULLIDNO) %>% 
  pull(ptID)

#Subset unrelated
#Subset to RSTR vs LTBI in TB samples
datV.tb.unrelated <- dat.combined.voom
datV.tb.unrelated$targets <- datV.tb.unrelated$targets %>% 
  filter(FULLIDNO %in% unrelated)
datV.tb.unrelated$E <- as.data.frame(datV.tb.unrelated$E) %>% 
  select(all_of(datV.tb.unrelated$targets$libID))
# identical(colnames(datV.tb.unrelated$E), datV.tb.unrelated$targets$libID)
datV.tb.unrelated$weights <- as.data.frame(datV.tb.unrelated$weights) %>% 
  select(all_of(datV.tb.unrelated$targets$libID)) %>% 
  as.matrix()

kin.unrelated <- kin[unrelated,unrelated]

#Subset related
#Subset to TB vs TB in LTBI samples
datV.tb.related <- dat.combined.voom
datV.tb.related$targets <- datV.tb.related$targets %>% 
  filter(FULLIDNO %in% related)
datV.tb.related$E <- as.data.frame(datV.tb.related$E) %>% 
  select(all_of(datV.tb.related$targets$libID))
# identical(colnames(datV.tb.related$E), datV.tb.related$targets$libID)
datV.tb.related$weights <- as.data.frame(datV.tb.related$weights) %>% 
  select(all_of(datV.tb.related$targets$libID)) %>% 
  as.matrix()

kin.related <- kin[related,related]

#### lm ####
unrelated.lmerel <- kmFit(dat = datV.tb.unrelated, kin=kin.unrelated,
                          patientID = "FULLIDNO",
                          model = "~ condition*Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + experiment + (1|FULLIDNO)", run.lmerel = TRUE, run.contrast = TRUE,
                          use.weights = TRUE, metrics = TRUE,
                          contrast.var = "condition:Sample_Group")

unrelated.lmerel$lmerel <- unrelated.lmerel$lmerel %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "voom weights",
         subset="unrelated")
unrelated.lmerel$lmerel.fit <- unrelated.lmerel$lmerel.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "voom weights",
         subset="unrelated")

save(unrelated.lmerel, file = "results/rstr_interaction_fit.RData")

###
related.lmerel <- kmFit(dat = datV.tb.related, kin=kin.related,
                        patientID = "FULLIDNO",
                        model = "~ condition*Sample_Group + KCHCA_AGE_YR_CURRENT + M0_KCVSEX + experiment + (1|FULLIDNO)", run.lmerel = TRUE, run.contrast = TRUE,
                        use.weights = TRUE, metrics = TRUE,
                        contrast.var = "condition:Sample_Group")

related.lmerel$lmerel <- related.lmerel$lmerel %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "voom weights",
         subset="related")
related.lmerel$lmerel.fit <- related.lmerel$lmerel.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "voom weights",
         subset="related")

save(unrelated.lmerel, related.lmerel, file = "results/rstr_interaction_fit.RData")
