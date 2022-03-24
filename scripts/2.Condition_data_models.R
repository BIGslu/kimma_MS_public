library(tidyverse)
library(edgeR)
library(limma)
library(kimma)
library(DESeq2)
library(variancePartition)

#### Data ####
attach("data/RSTR_data_all.RData")

#Subset to MEDIA vs TB in LTBI samples
dat.ltbi <- dat
dat.ltbi$samples <- dat.ltbi$samples %>% 
  filter(Sample_Group == "LTBI")
dat.ltbi$counts <- as.data.frame(dat.ltbi$counts) %>% 
  select(all_of(dat.ltbi$samples$libID))
# identical(colnames(dat.ltbi$counts), dat.ltbi$samples$libID)

datV.ltbi <- datV
datV.ltbi$targets <- datV.ltbi$targets %>% 
  filter(Sample_Group == "LTBI")
datV.ltbi$E <- as.data.frame(datV.ltbi$E) %>% 
  select(all_of(datV.ltbi$targets$libID))
# identical(colnames(datV.ltbi$E), datV.ltbi$targets$libID)
datV.ltbi$weights <- as.data.frame(datV.ltbi$weights) %>% 
  select(all_of(datV.ltbi$targets$libID)) %>% 
  as.matrix()

datD.ltbi <- datD
datD.ltbi$targets <- datD.ltbi$targets %>% 
  filter(Sample_Group == "LTBI")
datD.ltbi$E <- as.data.frame(datD.ltbi$E) %>% 
  select(all_of(datD.ltbi$targets$libID))
# identical(colnames(datD.ltbi$E), datD.ltbi$targets$libID)
datD.ltbi$weights <- as.data.frame(datD.ltbi$weights) %>% 
  select(all_of(datD.ltbi$targets$libID)) %>% 
  as.matrix()

#### kimma ####
#### Unpaired ####
# Paired = N, kinship = N, weights = N
kimma_nnn_condition <- kmFit(dat = datV.ltbi,
                             model = "~ condition", run.lm = TRUE, 
                             use.weights = FALSE)

kimma_nnn_condition$lm <- kimma_nnn_condition$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights")

# Paired = N, kinship = N, weights = Y
kimma_nny_condition <- kmFit(dat = datV.ltbi,
                             model = "~ condition", run.lm = TRUE, 
                             use.weights = TRUE)
kimma_nny_condition$lm <- kimma_nnn_condition$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "weights")

#### Paired ####
# Paired = Y, kinship = N, weights = N
kimma_ynn_condition <- kmFit(dat = datV.ltbi,
                             model = "~ condition + (1|ptID)", run.lme = TRUE, 
                             use.weights = FALSE)

kimma_ynn_condition$lme <- kimma_ynn_condition$lme %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "no kinship", weights = "no weights")

# Paired = Y, kinship = N, weights = Y
kimma_yny_condition <- kmFit(dat = datV.ltbi,
                             model = "~ condition + (1|ptID)", run.lme = TRUE, 
                             use.weights = TRUE)
kimma_yny_condition$lme <- kimma_ynn_condition$lme %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "no kinship", weights = "weights")

# Paired = Y, kinship = Y, weights = N
kimma_yyn_condition <- kmFit(dat = datV.ltbi, kin=kin,
                             model = "~ condition + (1|ptID)", run.lmekin = TRUE, 
                             use.weights = FALSE)
kimma_yyn_condition$lmekin <- kimma_yyn_condition$lmekin %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "no weights")

# Paired = Y, kinship = Y, weights = Y
kimma_yyy_condition <- kmFit(dat = datV.ltbi, kin=kin,
                             model = "~ condition + (1|ptID)", run.lmekin = TRUE, 
                             use.weights = TRUE)
kimma_yyy_condition$lmekin <- kimma_yyy_condition$lmekin %>% 
  mutate(software = "kimma", paired = "paired",
         kinship = "kinship", weights = "weights")

#### limma ####
mm <- model.matrix(~ condition, data=datV.ltbi$targets)
#### unpaired ####
# Paired = N, kinship = N, weights = N
datV.ltbi.noW <- datV.ltbi
datV.ltbi.noW$weights <- NULL

fit <- lmFit(object = datV.ltbi.noW, design = mm)
efit <- eBayes(fit)
limma_nnn_condition  <- extract_lmFit(design=mm, fit=efit) %>% 
  mutate(software = "limma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights")

# Paired = N, kinship = N, weights = Y
fit <- lmFit(object = datV.ltbi, design = mm,
             weights =  datV.ltbi$weights)
efit <- eBayes(fit)
limma_nny_condition <- extract_lmFit(design=mm, fit=efit) %>% 
  mutate(software = "limma", paired = "unpaired",
         kinship = "no kinship", weights = "weights")

#### paired ####
# Paired = Y, kinship = N, weights = N
dupCor <- duplicateCorrelation(object = datV.ltbi.noW$E,
                               design = mm,
                               block = datV.ltbi.noW$targets$ptID)
efit <- eBayes(lmFit(object = datV.ltbi.noW$E, design = mm,
                     block=datV.ltbi.noW$targets$ptID,
                     correlation=dupCor$consensus.correlation))
limma_ynn_condition <- extract_lmFit(design = mm, fit = efit)  %>% 
  mutate(software = "limma", paired = "paired",
         kinship = "no kinship", weights = "no weights")

# Paired = Y, kinship = N, weights = Y
dupCor <- duplicateCorrelation(object = datV.ltbi$E,
                               design = mm,
                               block = datV.ltbi$targets$ptID)
efit <- eBayes(lmFit(object = datV.ltbi$E, design = mm,
                     block=datV.ltbi$targets$ptID,
                     correlation=dupCor$consensus.correlation,
                     weights = datV.ltbi$weights))
limma_yny_condition <- extract_lmFit(design = mm, fit = efit)  %>% 
  mutate(software = "limma", paired = "paired",
         kinship = "no kinship", weights = "weights")

#### DESeq2 #### 
# Paired = N, kinship = N, weights = N
datS <- DESeqDataSetFromMatrix(countData = round(dat.ltbi$counts,digits=0),
                               colData = dat.ltbi$samples,
                               design = ~ condition)

fit <- DESeq(datS, parallel =  TRUE, BPPARAM = MulticoreParam(6))
deseq2_nnn_condition <- as.data.frame(results(fit)) %>% 
  dplyr::rename(pval=pvalue, FDR=padj) %>% 
  rownames_to_column("gene") %>% 
  mutate(variable = "condition") %>% 
  mutate(software = "DESeq2", paired = "unpaired",
         kinship = "no kinship", weights = "no weights")

# Paired = Y, kinship = N, weights = N
datS <- DESeqDataSetFromMatrix(countData = round(dat.ltbi$counts,digits=0),
                               colData = dat.ltbi$samples,
                               design = ~ condition + ptID)

fit <- DESeq(datS, parallel =  TRUE, BPPARAM = MulticoreParam(6))
deseq2_ynn_condition <- as.data.frame(results(fit)) %>% 
  dplyr::rename(pval=pvalue, FDR=padj) %>% 
  rownames_to_column("gene") %>% 
  mutate(variable = "condition") %>% 
  mutate(software = "DESeq2", paired = "paired",
         kinship = "no kinship", weights = "no weights")

#### dream ####
# Paired = Y, kinship = N, weights = N
datD.ltbi.noW <- datD.ltbi
datD.ltbi.noW$weights <- NULL

fit <- dream(datD.ltbi.noW, ~condition + (1|ptID),
                  useWeights = FALSE,
                  BPPARAM = SnowParam(6),
             datD.ltbi.noW$targets, REML=TRUE)
efit <- eBayes(fit)
dream_ynn_condition <- extract_lmFit(design = fit$design, fit=efit) %>% 
  mutate(geneName=symbol) %>% 
  mutate(software = "dream", paired = "paired",
         kinship = "no kinship", weights = "no weights")

# Paired = Y, kinship = N, weights = Y
fit <- dream(datD.ltbi, ~condition + (1|ptID),
             useWeights = TRUE,
             BPPARAM = SnowParam(6),
             datD.ltbi$targets, REML=TRUE)
efit <- eBayes(fit)
dream_yny_condition <- extract_lmFit(design = fit$design, fit=efit) %>% 
  mutate(geneName=symbol) %>% 
  mutate(software = "dream", paired = "paired",
         kinship = "no kinship", weights = "weights")

#### Save ####
condition.ls <- list()

for(d in c("kimma_nnn_condition", "kimma_nny_condition", 
           "kimma_ynn_condition", "kimma_yny_condition",
           "kimma_yyn_condition", "kimma_yyy_condition",
           "limma_nnn_condition", "limma_nny_condition",
           "limma_ynn_condition", "limma_yny_condition",
           "deseq2_nnn_condition", "deseq2_ynn_condition",
           "dream_ynn_condition", "dream_yny_condition")){
  print(d)
  condition.ls[[d]] <- get(d)
}

save(condition.ls, file="results/condition_fit.RData")
