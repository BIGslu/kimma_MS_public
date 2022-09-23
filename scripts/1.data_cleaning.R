library(tidyverse)
library(stringr)
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
## filter samples with kinship
overlap <- sort(intersect(colnames(kin.matrix), datV.dedup$targets$FULLIDNO))

## edgeR
dat.dedup.kin <- dat.dedup

dat.dedup.kin$samples <- dat.dedup.kin$samples %>% 
  filter(FULLIDNO %in% overlap)

dat.dedup.kin$counts <- as.data.frame(dat.dedup.kin$counts) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(dat.dedup.kin$samples$libID)) %>% 
  column_to_rownames()

## voom
datV.dedup.kin <- datV.dedup

datV.dedup.kin$targets <- datV.dedup.kin$targets %>% 
  filter(FULLIDNO %in% overlap)

datV.dedup.kin$E <- as.data.frame(datV.dedup.kin$E) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datV.dedup.kin$targets$libID)) %>% 
  column_to_rownames()

datV.dedup.kin$weights <- as.data.frame(datV.dedup.kin$weights) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datV.dedup.kin$targets$libID)) %>% 
  column_to_rownames()

## dream
datD.dedup.kin <- datD.dedup

datD.dedup.kin$targets <- datD.dedup.kin$targets %>% 
  filter(FULLIDNO %in% overlap) 

datD.dedup.kin$E <- as.data.frame(datD.dedup.kin$E) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datD.dedup.kin$targets$libID)) %>% 
  column_to_rownames()

datD.dedup.kin$weights <- as.data.frame(datD.dedup.kin$weights) %>% 
  rownames_to_column() %>% 
  select(rowname, all_of(datD.dedup.kin$targets$libID)) %>% 
  column_to_rownames()

#### De-identify ####
## Create fake ptID
all.ID <- unique(dat.dedup.kin$samples$FULLIDNO)
ID <- data.frame(FULLIDNO = all.ID,
                 ptID = paste0("pt",
                               str_pad(sample.int(length(all.ID), 
                                                  length(all.ID)),2, pad="0")))
## Save ID key
write_csv(file = "data/RSTR_ID_key.csv", ID)

## edgeR
dat.dedup.kin.rename <- dat.dedup.kin

dat.dedup.kin.rename$samples <- dat.dedup.kin.rename$samples %>% 
  select(group, lib.size, norm.factors, libID, FULLIDNO, condition, Sample_Group) %>% 
  full_join(ID) %>% 
  rename(libID2=libID) %>% 
  mutate(libID = paste(condition, ptID, sep="_"))
rownames(dat.dedup.kin.rename$samples) <- dat.dedup.kin.rename$samples$libID

identical(dat.dedup.kin.rename$samples$libID2, colnames(dat.dedup.kin.rename$counts))
colnames(dat.dedup.kin.rename$counts) <- dat.dedup.kin.rename$samples$libID

### voom
datV.dedup.kin.rename <- datV.dedup.kin

datV.dedup.kin.rename$targets <- datV.dedup.kin.rename$targets %>% 
  select(group, lib.size, norm.factors, libID, FULLIDNO, condition, Sample_Group) %>% 
  full_join(ID) %>% 
  rename(libID2=libID) %>% 
  mutate(libID = paste(condition, ptID, sep="_"))
rownames(datV.dedup.kin.rename$targets) <- datV.dedup.kin.rename$targets$libID

identical(datV.dedup.kin.rename$targets$libID2, colnames(datV.dedup.kin.rename$E))
colnames(datV.dedup.kin.rename$E) <- datV.dedup.kin.rename$targets$libID

identical(datV.dedup.kin.rename$targets$libID2, colnames(datV.dedup.kin.rename$weights))
colnames(datV.dedup.kin.rename$weights) <- datV.dedup.kin.rename$targets$libID

### dream
datD.dedup.kin.rename <- datD.dedup.kin

datD.dedup.kin.rename$targets <- datD.dedup.kin.rename$targets %>% 
  select(group, lib.size, norm.factors, libID, FULLIDNO, condition, Sample_Group) %>% 
  full_join(ID) %>% 
  rename(libID2=libID) %>% 
  mutate(libID = paste(condition, ptID, sep="_"))
rownames(datD.dedup.kin.rename$targets) <- datD.dedup.kin.rename$targets$libID

identical(datD.dedup.kin.rename$targets$libID2, colnames(datD.dedup.kin.rename$E))
colnames(datD.dedup.kin.rename$E) <- datD.dedup.kin.rename$targets$libID

identical(datD.dedup.kin.rename$targets$libID2, colnames(datD.dedup.kin.rename$weights))
colnames(datD.dedup.kin.rename$weights) <- datD.dedup.kin.rename$targets$libID

### kinship
kin.matrix.rename <- as.data.frame(kin.matrix) %>% 
  rownames_to_column("FULLIDNO") %>% 
  full_join(ID) %>% 
  select(-FULLIDNO) %>% 
  pivot_longer(-ptID, names_to = "FULLIDNO") %>% 
  full_join(ID, by="FULLIDNO") %>% 
  select(-FULLIDNO) %>%
  arrange(ptID.y) %>% 
  pivot_wider(names_from = ptID.y) %>% 
  arrange(ptID.x) %>% 
  column_to_rownames("ptID.x") %>% 
  as.matrix()

#### Remove old libID ####
dat.dedup.kin.rename$samples <- dat.dedup.kin.rename$samples %>% 
  select(-libID2, -FULLIDNO) %>% 
  arrange(libID)

datV.dedup.kin.rename$targets <- datV.dedup.kin.rename$targets %>% 
  select(-libID2, -FULLIDNO) %>% 
  arrange(libID)

datD.dedup.kin.rename$targets <- datD.dedup.kin.rename$targets %>% 
  select(-libID2, -FULLIDNO) %>% 
  arrange(libID)

#### Return to matrices as needed #####
datV.dedup.kin.rename$weights <- as.matrix(datV.dedup.kin.rename$weights)
datD.dedup.kin.rename$weights <- as.matrix(datD.dedup.kin.rename$weights)

#### Rename and save ####
dat <- dat.dedup.kin.rename
save(dat, file="data/RSTR_kimma_dat.RData")

datV <- datV.dedup.kin.rename
save(datV, file="data/RSTR_kimma_voom.RData")

datD <- datD.dedup.kin.rename
save(datD, file="data/RSTR_kimma_dream.RData")

kin <- kin.matrix.rename
save(kin, file="data/RSTR_kin.RData")

save(dat, datV, datD, kin, file="data/RSTR_data_all.RData")

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

identical(rownames(kin), unique(dat$samples$ptID))
identical(rownames(kin), unique(datV$targets$ptID))
identical(rownames(kin), unique(datD$targets$ptID))
