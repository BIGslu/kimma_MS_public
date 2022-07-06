library(tidyverse)
library(limma)
library(kimma)

#### Data ####
attach("data/RSTR_data_all.RData")
load("data/RSTR_snp.RData")

#SNP data
dat.snp <- geno.final %>% 
  dplyr::select(-gene) %>% 
  distinct() %>% 
  column_to_rownames("refsnp_id") %>% 
  t() %>% as.data.frame()

dat.map <- geno.final %>% 
  distinct(refsnp_id, gene)

#Subset related individuals
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

#Subset related
#Subset to TB vs TB in LTBI samples
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

#### eQTL model ####
eqtl <- kmFit_eQTL(dat.snp=dat.snp, dat.map = dat.map,
                   genotypeID="refsnp_id",
                   dat = datV.tb.related, kin=kin.related,
                   model = "~ genotype + (1|ptID)", run.lmerel = TRUE, 
                   use.weights = TRUE, metrics = TRUE)

eqtl$lmerel <- eqtl$lmerel %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "voom weights",
         subset="related")
eqtl$lmerel.fit <- eqtl$lmerel.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "voom weights",
         subset="related")

save(eqtl, file="results/rstr_subset_eqtl.RData")
