library(tidyverse)
library(BIGpicture)

#### LIMMA ####
attach("results/model_fit/LTBI_limma.RData")
LTBI.efit.noW <- LTBI.efit
LTBI.results.limma.noW <- LTBI.results.limma
attach("results/model_fit/LTBI_limma_weights.RData")
LTBI.efit.W <- efit
LTBI.results.limma.W <- results_limma_LTBI_weights

## Venn
venn.ls <- list()
venn.ls[["noW"]] <- LTBI.results.limma.noW %>% 
  filter(variable == "conditionTB" & adj.P.Val < 0.05) %>% 
  distinct(symbol) %>% unlist()
venn.ls[["W"]] <- LTBI.results.limma.W %>% 
  filter(variable == "conditionTB" & adj.P.Val < 0.05) %>% 
  distinct(symbol) %>% unlist()

#### KIMMA ####
attach("results/model_fit/LTBI_kimma.RData")
LTBI.efit.k.noW <- kimma.LTBI
LTBI.results.kimma.noW <- LTBI.results.kimma
attach("results/model_fit/LTBI_weight_kimma.RData")
LTBI.efit.k.W <- kimma.LTBI.wt
LTBI.results.kimma.W <- LTBI.wt.results.kimma

## Venn
venn.ls2 <- list()
venn.ls2[["noW"]] <- LTBI.results.kimma.noW %>% 
  filter(variable == "conditionTB" & FDR < 0.05) %>% 
  distinct(gene) %>% unlist()
venn.ls2[["W"]] <- LTBI.results.kimma.W %>% 
  filter(variable == "conditionTB" & FDR < 0.05) %>% 
  distinct(gene) %>% unlist()


## Weights should have an effect on a minority of genes
#Overlap within software, between models
par(mfrow=c(2,1))
venn::venn(venn.ls)
title(sub="limma", line = -1)
venn::venn(venn.ls2)
title(sub="kimma", line = -1)
dev.off()

##If everything is the same, 100% overlap
## But we have slightly different pval estimaters in kimma and limma
#Overlap within the same model, between software
par(mfrow=c(2,1))
venn::venn(list(venn.ls[[1]],venn.ls2[[1]]),
           snames = c("limma","kimma"))
title(sub="noW", line = -1)
venn::venn(list(venn.ls[[2]],venn.ls2[[2]]),
           snames = c("limma","kimma"))
title(sub="W", line = -1)
dev.off()

##If only kimma can use BIGpicture fxn
plot_venn_genes(LTBI.efit.k.noW, model="lm")

##### Paired samples #####
load("data/RSTR_RNAseq_data_combined_uniqueFULLID_weights.RData")
library(limma)
library(kimma)

#Limma
mm <- model.matrix(~condition, data = dat.combined.voom$targets)
#Block by donor - paired sample calc
dupCor <- duplicateCorrelation(object = dat.combined.voom$E,
                               design = mm,
                               block = dat.combined.voom$targets$FULLIDNO)
#Mean correlation of expression of all genes within a sample
dupCor$consensus.correlation

#run model
fit <- lmFit(dat.combined.voom$E, 
             mm, 
             block=dat.combined.voom$targets$FULLIDNO,
             correlation=dupCor$consensus.correlation)
efit <- eBayes(fit)

#kimma
k.fit <- kmFit(dat.combined.voom, patientID = "FULLIDNO",
               use.weights = TRUE/FALSE,
               model = "~ condition + (1|FULLIDNO)", run.lme = TRUE)

##### Kinship #####
k.fit <- kmFit(dat.combined.voom, patientID = "FULLIDNO",
               kin=kinship.mat,
               use.weights = TRUE/FALSE,
               model = "~ condition + (1|FULLIDNO)", run.lmekin = TRUE)