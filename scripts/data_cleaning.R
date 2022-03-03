library(tidyverse)
#voom object
library(edgeR)
library(limma)
#dream object
library(variancePartition)

#### Load data ####
attach("data/RSTR_RNAseq_data_combined_clean.RData")
dat.raw <- dat.combined
datV.raw <- dat.combined.voom
datV.raw$design <- NULL

load("data/RSTR_RNAseq_kinship.RData")

#### Calculate dream weights ####
datD.raw <- voomWithDreamWeights(dat.raw, formula = ~condition + (1|FULLIDNO),
                             data = dat.raw$samples)

#### Remove duplicate FULLIDNO ####
# List duplicates
dups2 <- dat.raw$samples %>% 
  dplyr::count(FULLIDNO, condition) %>% 
  filter(n>1)

# list RSID with the most media sequences per FULLIDNO donor
to.keep <- dat.raw$samples %>% 
  #media samples
  filter(condition == "MEDIA") %>% 
  #Keep library with most sequences
  group_by(FULLIDNO) %>% 
  slice_max(lib.size) %>% 
  ungroup() %>% 
  select(RSID) %>% unlist(use.names = FALSE)

# subset dat objects
## edgeR
dat.dedup <- dat.raw

dat.dedup$samples <- dat.dedup$samples %>% 
  filter(RSID %in% to.keep)

dat.dedup$counts <- as.data.frame(dat.dedup$counts) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(dat.dedup$samples$libID)) %>% 
  column_to_rownames()

## voom
datV.dedup <- datV.raw

datV.dedup$targets <- datV.dedup$targets %>% 
  filter(RSID %in% to.keep)

datV.dedup$E <- as.data.frame(datV.dedup$E) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datV.dedup$targets$libID)) %>% 
  column_to_rownames()

rownames(datV.dedup$weights) <- rownames(datV.raw$E)
colnames(datV.dedup$weights) <- colnames(datV.raw$E)

datV.dedup$weights <- as.data.frame(datV.dedup$weights) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datV.dedup$targets$libID)) %>% 
  column_to_rownames()

## dream
datD.dedup <- datD.raw

datD.dedup$targets <- datD.dedup$targets %>% 
  filter(RSID %in% to.keep)

datD.dedup$E <- as.data.frame(datD.dedup$E) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datD.dedup$targets$libID)) %>% 
  column_to_rownames()

rownames(datD.dedup$weights) <- rownames(datD.raw$E)
colnames(datD.dedup$weights) <- colnames(datD.raw$E)
datD.dedup$weights <- as.data.frame(datD.dedup$weights) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datD.dedup$targets$libID)) %>% 
  column_to_rownames()

#### Remove missing kinship ####
## rename dat to FULLIDNO
## edgeR
dat.dedup.rename <- dat.dedup

dat.dedup.rename$samples <- dat.dedup.rename$samples %>% 
  select(group, lib.size, norm.factors, libID, FULLIDNO, condition, Sample_Group) %>% 
  rename(ptID=FULLIDNO, libID2=libID) %>% 
  mutate(libID = paste(condition, ptID, sep="_"))
rownames(dat.dedup.rename$samples) <- dat.dedup.rename$samples$libID

identical(dat.dedup.rename$samples$libID2, colnames(dat.dedup.rename$counts))
colnames(dat.dedup.rename$counts) <- dat.dedup.rename$samples$libID

### voom
datV.dedup.rename <- datV.dedup

datV.dedup.rename$targets <- datV.dedup.rename$targets %>% 
  select(group, lib.size, norm.factors, libID, FULLIDNO, condition, Sample_Group) %>% 
  rename(ptID=FULLIDNO, libID2=libID) %>% 
  mutate(libID = paste(condition, ptID, sep="_"))
rownames(datV.dedup.rename$targets) <- datV.dedup.rename$targets$libID

identical(datV.dedup.rename$targets$libID2, colnames(datV.dedup.rename$E))
colnames(datV.dedup.rename$E) <- datV.dedup.rename$targets$libID

identical(datV.dedup.rename$targets$libID2, colnames(datV.dedup.rename$weights))
colnames(datV.dedup.rename$weights) <- datV.dedup.rename$targets$libID

### dream
datD.dedup.rename <- datD.dedup

datD.dedup.rename$targets <- datD.dedup.rename$targets %>% 
  select(group, lib.size, norm.factors, libID, FULLIDNO, condition, Sample_Group) %>% 
  rename(ptID=FULLIDNO, libID2=libID) %>% 
  mutate(libID = paste(condition, ptID, sep="_"))
rownames(datD.dedup.rename$targets) <- datD.dedup.rename$targets$libID

identical(datD.dedup.rename$targets$libID2, colnames(datD.dedup.rename$E))
colnames(datD.dedup.rename$E) <- datD.dedup.rename$targets$libID

identical(datD.dedup.rename$targets$libID2, colnames(datD.dedup.rename$weights))
colnames(datD.dedup.rename$weights) <- datD.dedup.rename$targets$libID

## filter samples with kinship
overlap <- intersect(colnames(kin.matrix), datV.dedup.rename$targets$ptID)

## edgeR
dat.dedup.rename.kin <- dat.dedup.rename

dat.dedup.rename.kin$samples <- dat.dedup.rename.kin$samples %>% 
  filter(ptID %in% overlap) %>% 
  select(-libID2)

dat.dedup.rename.kin$counts <- as.data.frame(dat.dedup.rename.kin$counts) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(dat.dedup.rename.kin$samples$libID)) %>% 
  column_to_rownames()

## voom
datV.dedup.rename.kin <- datV.dedup.rename

datV.dedup.rename.kin$targets <- datV.dedup.rename.kin$targets %>% 
  filter(ptID %in% overlap) %>% 
  select(-libID2)

datV.dedup.rename.kin$E <- as.data.frame(datV.dedup.rename.kin$E) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datV.dedup.rename.kin$targets$libID)) %>% 
  column_to_rownames()

datV.dedup.rename.kin$weights <- as.data.frame(datV.dedup.rename.kin$weights) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datV.dedup.rename.kin$targets$libID)) %>% 
  column_to_rownames()

## dream
datD.dedup.rename.kin <- datD.dedup.rename

datD.dedup.rename.kin$targets <- datD.dedup.rename.kin$targets %>% 
  filter(ptID %in% overlap) %>% 
  select(-libID2)

datD.dedup.rename.kin$E <- as.data.frame(datD.dedup.rename.kin$E) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datD.dedup.rename.kin$targets$libID)) %>% 
  column_to_rownames()

datD.dedup.rename.kin$weights <- as.data.frame(datD.dedup.rename.kin$weights) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datD.dedup.rename.kin$targets$libID)) %>% 
  column_to_rownames()

#### Return to matrices as needed #####
datV.dedup.rename.kin$weights <- as.matrix(datV.dedup.rename.kin$weights)
datD.dedup.rename.kin$weights <- as.matrix(datD.dedup.rename.kin$weights)

#### Rename and save ####
dat <- dat.dedup.rename.kin
save(dat, file="data/RSTR_kimma_dat.RData")

datV <- datV.dedup.rename.kin
save(datV, file="data/RSTR_kimma_voom.RData")

datD <- datD.dedup.rename.kin
save(datD, file="data/RSTR_kimma_dream.RData")

# Checks
dim(dat$counts)
dim(dat$samples)
dim(dat$genes)

dim(datV$E)
dim(datV$targets)
dim(datV$genes)
dim(datV$weights)

dim(datD$E)
dim(datD$targets)
dim(datD$genes)
dim(datD$weights)
