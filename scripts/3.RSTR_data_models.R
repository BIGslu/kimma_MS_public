library(tidyverse)
library(edgeR)
library(limma)
library(kimma)
library(DESeq2)
library(variancePartition)

#### Data ####
attach("data/RSTR_data_all.RData")
kin <- kin

#Subset to RSTR vs LTBI in media samples
dat.media <- dat
dat.media$samples <- dat.media$samples %>% 
  filter(condition == "MEDIA")
dat.media$counts <- as.data.frame(dat.media$counts) %>% 
  select(all_of(dat.media$samples$libID))
# identical(colnames(dat.media$counts), dat.media$samples$libID)

datV.media <- datV
datV.media$targets <- datV.media$targets %>% 
  filter(condition == "MEDIA")
datV.media$E <- as.data.frame(datV.media$E) %>% 
  select(all_of(datV.media$targets$libID))
# identical(colnames(datV.media$E), datV.media$targets$libID)
datV.media$weights <- as.data.frame(datV.media$weights) %>% 
  select(all_of(datV.media$targets$libID)) %>% 
  as.matrix()

#### kimma ####
#### Unpaired ####
# Paired = N, kinship = N, weights = N
kimma_nnn_rstr <- kmFit(dat = datV.media,
                        model = "~ Sample_Group", run.lm = TRUE, 
                        use.weights = FALSE)

kimma_nnn_rstr$lm <- kimma_nnn_rstr$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights")

# Paired = N, kinship = N, weights = Y
kimma_nny_rstr <- kmFit(dat = datV.media,
                        model = "~ Sample_Group", run.lm = TRUE, 
                        use.weights = TRUE)
kimma_nny_rstr$lm <- kimma_nnn_rstr$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "weights")

# Paired = N, kinship = Y, weights = N
kimma_nyn_rstr <- kmFit(dat = datV.media, kin=kin,
                        model = "~ Sample_Group + (1|ptID)", run.lmekin = TRUE, 
                        use.weights = FALSE, processors = 4)
kimma_nyn_rstr$lmekin <- kimma_nyn_rstr$lmekin %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "no weights")

# Paired = N, kinship = Y, weights = Y
kimma_nyy_rstr <- kmFit(dat = datV.media, kin=kin,
                        model = "~ Sample_Group + (1|ptID)", run.lmekin = TRUE, 
                        use.weights = TRUE, processors = 4)
kimma_nyy_rstr$lmekin <- kimma_nyy_rstr$lmekin %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "weights")

#### limma ####
mm <- model.matrix(~ Sample_Group, data=datV.media$targets)

# Paired = N, kinship = N, weights = N
datV.media.noW <- datV.media
datV.media.noW$weights <- NULL

fit <- lmFit(object = datV.media.noW, design = mm)
efit <- eBayes(fit)
limma_nnn_rstr  <- extract_lmFit(design=mm, fit=efit) %>% 
  mutate(software = "limma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights")

# Paired = N, kinship = N, weights = Y
fit <- lmFit(object = datV.media, design = mm,
             weights =  datV.media$weights)
efit <- eBayes(fit)
limma_nny_rstr  <- extract_lmFit(design=mm, fit=efit) %>% 
  mutate(software = "limma", paired = "unpaired",
         kinship = "no kinship", weights = "weights")

#### DESeq2 #### 
# Paired = N, kinship = N, weights = N
datS <- DESeqDataSetFromMatrix(countData = round(dat.media$counts,digits=0),
                              colData = dat.media$samples,
                              design = ~ Sample_Group)

fit <- DESeq(datS, parallel =  TRUE, BPPARAM = MulticoreParam(6))
deseq2_nnn_rstr <- as.data.frame(results(fit)) %>% 
  dplyr::rename(pval=pvalue, FDR=padj) %>% 
  rownames_to_column("gene") %>% 
  mutate(variable = "Sample_Group") %>% 
  mutate(software = "DESeq2", paired = "unpaired",
         kinship = "no kinship", weights = "no weights")

#### Save ####
rstr.ls <- list()

for(d in c("kimma_nnn_rstr", "kimma_nny_rstr", 
           "kimma_nyn_rstr", "kimma_nyy_rstr",
           "limma_nnn_rstr", "limma_nny_rstr",
           "deseq2_nnn_rstr")){
  print(d)
  rstr.ls[[d]] <- get(d)
}

save(rstr.ls, file="results/rstr_fit.RData")
