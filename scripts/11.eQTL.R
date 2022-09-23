library(tidyverse)
library(limma)
library(kimma)

#### Data ####
attach("data/RSTR_data_all.RData")
load("data/RSTR_snp.RData")

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

#### eQTL model ####
dat.snp <- geno.final %>% 
  select(-gene) %>% 
  distinct() %>% 
  column_to_rownames("refsnp_id") %>% 
  t() %>% as.data.frame
dat.map <- geno.final %>% 
  select(refsnp_id, gene)

eqtl <- kmFit_eQTL(dat.snp=dat.snp, dat.map=dat.map,
                   genotypeID="refsnp_id",
                   dat = datV.ltbi.related, kin=kin.related,
                   model = "~ condition*genotype + (1|ptID)", 
                   run.lmerel = TRUE, 
                   use.weights = TRUE, metrics = TRUE)

eqtl$lmerel <- eqtl$lmerel %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "voom weights",
         subset="related") %>% 
  mutate(FDR = p.adjust(pval, method = "BH"))
eqtl$lmerel.fit <- eqtl$lmerel.fit %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "voom weights",
         subset="related") 
save(eqtl, file="results/rstr_subset_eqtl.RData")

#### Summary numbers ####
#DEG
length(deg1)
#eQTL tested
nrow(dat.map)
#eQTL genes
length(unique(dat.map$gene))

#Signif eQTL
eqtl$lmerel %>% 
  filter(grepl(":",variable) & FDR < 0.05) %>% 
  nrow()
eqtl$lmerel %>% 
  filter(grepl(":",variable) & FDR < 0.05) %>% 
  pull(gene) %>% unique() %>% length()


eqtl$lmerel %>% 
  filter(grepl(":",variable)) %>% 
  slice_min(FDR)
