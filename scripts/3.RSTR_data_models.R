library(tidyverse)
library(edgeR)
library(limma)
library(kimma)
library(DESeq2)
library(variancePartition)

#### Data ####
attach("data/RSTR_data_all.RData")
kin <- kin

#Subset to RSTR vs LTBI in tb samples
dat.tb <- dat
dat.tb$samples <- dat.tb$samples %>% 
  filter(condition == "TB")
dat.tb$counts <- as.data.frame(dat.tb$counts) %>% 
  select(all_of(dat.tb$samples$libID))
# identical(colnames(dat.tb$counts), dat.tb$samples$libID)

datV.tb <- datV
datV.tb$targets <- datV.tb$targets %>% 
  filter(condition == "TB")
datV.tb$E <- as.data.frame(datV.tb$E) %>% 
  select(all_of(datV.tb$targets$libID))
# identical(colnames(datV.tb$E), datV.tb$targets$libID)
datV.tb$weights <- as.data.frame(datV.tb$weights) %>% 
  select(all_of(datV.tb$targets$libID)) %>% 
  as.matrix()

#### kimma ####
# Paired = N, kinship = N, weights = N
kimma_nnn_rstr <- kmFit(dat = datV.tb,
                        model = "~ Sample_Group", run.lm = TRUE, 
                        use.weights = FALSE, metrics = TRUE)

kimma_nnn_rstr$lm <- kimma_nnn_rstr$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights")
kimma_nnn_rstr$lm.fit <- kimma_nnn_rstr$lm.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights")

# Paired = N, kinship = N, weights = Y
kimma_nny_rstr <- kmFit(dat = datV.tb,
                        model = "~ Sample_Group", run.lm = TRUE, 
                        use.weights = TRUE, metrics = TRUE)
kimma_nny_rstr$lm <- kimma_nny_rstr$lm %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "voom weights")
kimma_nny_rstr$lm.fit <- kimma_nny_rstr$lm.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "no kinship", weights = "voom weights")

# Paired = N, kinship = Y, weights = N
kimma_nyn_rstr <- kmFit(dat = datV.tb, kin=kin,
                        model = "~ Sample_Group + (1|ptID)", run.lmerel = TRUE, 
                        use.weights = FALSE, processors = 4, metrics = TRUE)
kimma_nyn_rstr$lmerel <- kimma_nyn_rstr$lmerel %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "no weights")
kimma_nyn_rstr$lmerel.fit <- kimma_nyn_rstr$lmerel.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "no weights")

# Paired = N, kinship = Y, weights = Y
kimma_nyy_rstr <- kmFit(dat = datV.tb, kin=kin,
                        model = "~ Sample_Group + (1|ptID)", run.lmerel = TRUE, 
                        use.weights = TRUE, processors = 4, metrics = TRUE)
kimma_nyy_rstr$lmerel <- kimma_nyy_rstr$lmerel %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "voom weights")
kimma_nyy_rstr$lmerel.fit <- kimma_nyy_rstr$lmerel.fit %>% 
  mutate(software = "kimma", paired = "unpaired",
         kinship = "kinship", weights = "voom weights")

#### limma ####
mm <- model.matrix(~ Sample_Group, data=datV.tb$targets)

# Paired = N, kinship = N, weights = N
datV.tb.noW <- datV.tb
datV.tb.noW$weights <- NULL

fit <- lmFit(object = datV.tb.noW, design = mm)
efit <- eBayes(fit)
limma_nnn_rstr  <- extract_lmFit(design=mm, fit=efit) %>% 
  mutate(software = "limma", paired = "unpaired",
         kinship = "no kinship", weights = "no weights")

# Paired = N, kinship = N, weights = Y
fit <- lmFit(object = datV.tb, design = mm,
             weights =  datV.tb$weights)
efit <- eBayes(fit)
limma_nny_rstr  <- extract_lmFit(design=mm, fit=efit) %>% 
  mutate(software = "limma", paired = "unpaired",
         kinship = "no kinship", weights = "voom weights")

#### DESeq2 #### 
# Paired = N, kinship = N, weights = N
datS <- DESeqDataSetFromMatrix(countData = round(dat.tb$counts,digits=0),
                              colData = dat.tb$samples,
                              design = ~ Sample_Group)

fit <- DESeq(datS, parallel =  TRUE, BPPARAM = MulticoreParam(6))
deseq2_nnn_rstr <- as.data.frame(results(fit)) %>% 
  dplyr::rename(pval=pvalue, FDR=padj) %>% 
  rownames_to_column("gene") %>% 
  mutate(variable = "Sample_Group") %>% 
  mutate(software = "DESeq2", paired = "unpaired",
         kinship = "no kinship", weights = "no weights")

#### Save ####
#Save indiv results for easier updating in future
save(kimma_nnn_rstr,kimma_nny_rstr,kimma_nyn_rstr,kimma_nyy_rstr,
     limma_nnn_rstr,limma_nny_rstr,deseq2_nnn_rstr,
     file="results/model_fit/sample_group/Sample_Group_all_tb.RData")

#Combine and format 1 df
resultK <- bind_rows(kimma_nnn_rstr$lm, kimma_nny_rstr$lm,
                     kimma_nyn_rstr$lmerel, kimma_nyy_rstr$lmerel) %>% 
  mutate(variable = recode(variable, "Sample_GroupRSTR"="Sample_Group")) %>% 
  filter(variable == "Sample_Group") %>% 
  mutate(subset = "all") %>% 
  select(subset, software:weights,gene:FDR)

resultL <- bind_rows(limma_nnn_rstr,limma_nny_rstr) %>% 
  mutate(variable = recode(variable, "Sample_GroupRSTR"="Sample_Group")) %>% 
  filter(variable == "Sample_Group") %>% 
  mutate(subset = "all") %>% 
  #rename to match kimma results
  select(-gene) %>% 
  dplyr::rename(gene=symbol) %>% 
  select(all_of(colnames(resultK)))

resultS <- deseq2_nnn_rstr %>% 
  mutate(subset = "all") %>% 
  #rename to match kimma results
  dplyr::rename(estimate=log2FoldChange) %>% 
  select(all_of(colnames(resultK)))

rstr_result <- bind_rows(resultK, resultL, resultS)
rownames(rstr_result) <- NULL

#Combine fit into 1 df
rstr_metric <- bind_rows(kimma_nnn_rstr$lm.fit,kimma_nny_rstr$lm.fit,
                         kimma_nyn_rstr$lmerel.fit,kimma_nyy_rstr$lmerel.fit) %>% 
  mutate(subset = "all") %>% 
  select(subset, software:weights,gene:adj_Rsq)

#Save
save(rstr_result, rstr_metric, file="results/rstr_tb_fit.RData")
