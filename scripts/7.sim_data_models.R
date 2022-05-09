library(tidyverse)
library(limma)
library(kimma)
library(variancePartition)
library(DESeq2)
library(foreach)
library(parallel)
library(doParallel)
library(BiocParallel)
`%notin%` <- Negate(`%in%`)

#### Sensitivity and specificity fxn ####
ss_calc <- function(i, dat, fdr, fdr.name){
  deg_real <- dat$genes %>% 
    filter(sim.group != "none") %>% 
    select(symbol) %>% unlist(use.names = FALSE)
  
  df <- data.frame()
  for(p.cutoff in unique(c(seq(0,0.01, 0.001), seq(0,1, 0.01)))){
    if("P.Value" %in% colnames(fdr)){
    deg_test <- fdr %>% 
      filter(grepl("condition",variable) & P.Value < p.cutoff) %>% 
      select(geneName) %>% unlist(use.names = FALSE)
    deg_test2 <- fdr %>% 
      filter(grepl("condition",variable) & adj.P.Val < p.cutoff) %>% 
      select(geneName) %>% unlist(use.names = FALSE)
    } else{
      deg_test <- fdr %>% 
        filter(grepl("condition",variable) & pval < p.cutoff) %>% 
        select(gene) %>% unlist(use.names = FALSE)
      deg_test2 <- fdr %>% 
        filter(grepl("condition",variable) & FDR < p.cutoff) %>% 
        select(gene) %>% unlist(use.names = FALSE)
    }
    
    df <- data.frame(model = fdr.name,
                     simulation = i,
                     cutoff = p.cutoff,
                     TP_p = length(deg_test[deg_test %in% deg_real]),
                     FP_p = length(deg_test[deg_test %notin% deg_real]),
                     FN_p = length(deg_real[deg_real %notin% deg_test])) %>% 
      mutate(TN_p = nrow(dat$genes) - TP_p - FP_p - FN_p) %>% 
      mutate(TP_fdr = length(deg_test2[deg_test2 %in% deg_real]),
             FP_fdr = length(deg_test2[deg_test2 %notin% deg_real]),
             FN_fdr = length(deg_real[deg_real %notin% deg_test2]),
             TN_fdr = nrow(dat$genes) - TP_fdr - FP_fdr - FN_fdr) %>% 
      mutate(sensitivity_p = TP_p/(TP_p+FN_p),
             specificity_p = TN_p/(TN_p+FP_p),
             sensitivity_fdr = TP_fdr/(TP_fdr+FN_fdr),
             specificity_fdr = TN_fdr/(TN_fdr+FP_fdr),
             genes_p=list(deg_test),
             genes_fdr=list(deg_test2)) %>% 
      bind_rows(df)
  }
  return(df)
}

#### Combine lists fxn ####
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

#### Kinship data ####
load("data_sim/simulated_kinship.RData")

#### Limma and kimma ####
cl <- makeCluster(30)
doParallel::registerDoParallel(cl)
result <- list()

result <- foreach(i=1:100, .combine='comb', .multicombine=TRUE,
                  .init=list(list(), list()),
                  .packages = c("tidyverse","limma","kimma",
                                "foreach","doParallel")
                  ) %dopar% {
  print(i)
  #### Data ####
  load(paste0("data_sim/simulated",i,".RData"))
  dat.sim$counts <- NULL
  
  # Data without weights
  dat.sim.noW <- dat.sim
  dat.sim.noW$weights <- NULL
  # Data with voom weights
  dat.sim.voomW <- dat.sim
  dat.sim.voomW$weights_dream <- NULL
  # Data with dream weights
  dat.sim.voomD <- dat.sim
  dat.sim.voomD$weights <- dat.sim.voomD$weights_dream
  dat.sim.voomW$weights_dream <- NULL
  
  # Set to NULL to capture errors
  fdrKnnn <- fdrKnny <- fdrKynn <- fdrKyny <- fdrKyyn <- fdrKyyy <- fdrLnnn <- fdrLnny <- fdrLynn <- fdrLyny <- NULL
  
  # Run models
  #### Paired = N, kinship = N, weights = N ####
  print("NNN")
  ## limma
  start <- Sys.time()
  mm <- model.matrix(~condition, data = dat.sim.noW$targets)
  fitLnnn <- eBayes(lmFit(object = dat.sim.noW$E, design = mm))
  fdrLnnn <- extract_lmFit(design = mm, fit = fitLnnn)
  endLnnn <- Sys.time()

  ## kimma
  fdrKnnn <- kmFit(dat.sim.noW, use.weights = FALSE,
                 model = "~ condition", run.lm = TRUE,
                 processors=1)[[1]]
  endKnnn <- Sys.time()

  #### Paired = N, kinship = N, weights = Y ####
  print("NNY")
  ## limma
  mm <- model.matrix(~condition, data = dat.sim.voomW$targets)
  fitLnny <- eBayes(lmFit(object = dat.sim.voomW$E, design = mm,
                          weights=dat.sim.voomW$weights))
  fdrLnny <- extract_lmFit(design = mm, fit = fitLnny)
  endLnny <- Sys.time()

  ## kimma
  fdrKnny <- kmFit(dat.sim.voomW, use.weights = TRUE,
                   model = "~ condition", run.lm = TRUE,
                   processors=1)[[1]]
  endKnny <- Sys.time()
  
  fdrKnny2 <- kmFit(dat.sim.voomD, use.weights = TRUE,
                   model = "~ condition", run.lm = TRUE,
                   processors=1)[[1]]
  endKnny2 <- Sys.time()
  
  #### Paired = Y, kinship = N, weights = N ####
  print("YNN")
  ## limma
  mm <- model.matrix(~condition, data = dat.sim.noW$targets)
  dupCor <- duplicateCorrelation(object = dat.sim.noW$E,
                                 design = mm,
                                 block = dat.sim.noW$targets$ptID)
  fitLynn <- eBayes(lmFit(object = dat.sim.noW$E, design = mm,
                          block=dat.sim.noW$targets$ptID,
                          correlation=dupCor$consensus.correlation))
  fdrLynn <- extract_lmFit(design = mm, fit = fitLynn)
  endLynn <- Sys.time()

  ## kimma
  fdrKynn <- kmFit(dat.sim.noW, use.weights = FALSE,
                   model = "~ condition + (1|ptID)", run.lme = TRUE, 
                   processors=1)[[1]]
  endKynn <- Sys.time()
  
  #### Paired = Y, kinship = N, weights = Y ####
  print("YNY")
  ## limma
  mm <- model.matrix(~condition, data = dat.sim.voomW$targets)
  dupCor <- duplicateCorrelation(object = dat.sim.voomW$E,
                                 design = mm,
                                 block = dat.sim.voomW$targets$ptID)
  fitLyny <- eBayes(lmFit(object = dat.sim.voomW$E, design = mm,
                          block=dat.sim.voomW$targets$ptID,
                          correlation=dupCor$consensus.correlation,
                          weights=dat.sim.voomW$weights))
  fdrLyny <- extract_lmFit(design = mm, fit = fitLyny)
  endLyny <- Sys.time()
  
  ## kimma
  fdrKyny <- kmFit(dat.sim.voomW, use.weights = TRUE,
                   model = "~ condition + (1|ptID)", run.lme = TRUE, 
                   processors=1)[[1]]
  endKyny <- Sys.time()
  
  fdrKyny2 <- kmFit(dat.sim.voomD, use.weights = TRUE,
                   model = "~ condition + (1|ptID)", run.lme = TRUE, 
                   processors=1)[[1]]
  endKyny2 <- Sys.time()
  
  #### kinship = Y ####
  print("XYX")
  ## kimma only
  ### Replace kinship matrix to remove paired design
  # Paired = N, kinship = Y, weights = N
  fdrKnyn <- kmFit(dat.sim.noW, kin=kin.sim.unpair, use.weights = FALSE,
                   patientID = "libID",
                   model = "~ condition + (1|libID)", run.lmerel = TRUE,
                   processors=1)[[1]]
  endKnyn <- Sys.time()

  # Paired = N, kinship = Y, weights = Y
  fdrKnyy <- kmFit(dat.sim.voomW, kin=kin.sim.unpair, use.weights = TRUE,
                   patientID = "libID",
                   model = "~ condition + (1|libID)", run.lmerel = TRUE,
                   processors=1)[[1]]
  endKnyy <- Sys.time()

  fdrKnyy2 <- kmFit(dat.sim.voomD, kin=kin.sim.unpair, use.weights = TRUE,
                   patientID = "libID",
                   model = "~ condition + (1|libID)", run.lmerel = TRUE,
                   processors=1)[[1]]
  endKnyy2 <- Sys.time()
  
  # Paired = Y, kinship = Y, weights = N
  fdrKyyn <- kmFit(dat.sim.noW, kin=kin.sim, use.weights = FALSE,
                   model = "~ condition + (1|ptID)", run.lmerel = TRUE,
                   processors=1)[[1]]
  endKyyn <- Sys.time()

  # Paired = Y, kinship = Y, weights = Y
  fdrKyyy <- kmFit(dat.sim.voomW, kin=kin.sim, use.weights = TRUE,
                   model = "~ condition + (1|ptID)", run.lmerel = TRUE,
                   processors=1)[[1]]
  endKyyy <- Sys.time()
  
  fdrKyyy2 <- kmFit(dat.sim.voomD, kin=kin.sim, use.weights = TRUE,
                   model = "~ condition + (1|ptID)", run.lmerel = TRUE,
                   processors=1)[[1]]
  endKyyy2 <- Sys.time()
  
  #### sensitivity and specificity ####
  ss_result <- data.frame()
  for(m in ls(pattern = "^fdr")){
    if(!is.null(get(m))){
      ss_result <- ss_calc(i, dat.sim, get(m), m) %>% 
        bind_rows(ss_result)
    }
  }
  
  #### time ####
  time_result <- data.frame(model=c("Lnnn", "Knnn", "Lnny", "Knny", "Knny2",
                                    "Lynn", "Kynn", "Lyny", "Kyny", "Kyny2",
                                    "Knyn", "Knyy", "Knyy2",
                                    "Kyyn", "Kyyy", "Kyyy2"),
                            simulation=i,
                            time_s=c(endLnnn-start, endKnnn-endLnnn,
                                     endLnny-endKnnn, endKnny-endLnny,
                                     endKnny2-endKnny, endLynn-endKnny2, 
                                     endKynn-endLynn, endLyny-endKynn, 
                                     endKyny-endLyny, endKyny2-endKyny, 
                                     endKnyn-endKyny2, endKnyy-endKnyn,
                                     endKnyy2-endKnyy, endKyyn-endKnyy2,
                                     endKyyy-endKyyn, endKyyy2-endKyyy))
  
  #### Final results ####
  list(ss_result, time_result)
                  }

stopCluster(cl)
save(result, file="results/model_fit/simulated_limma_kimma.RData")

#### Multi-processor kimma ####
time_resultKM <- data.frame()
# Just for time
# ss_resultKM <- data.frame()
p = 6

for(i in 1:100){
  #### Data ####
  print(i)
  load(paste0("data_sim/simulated",i,".RData"))
  dat.sim$counts <- NULL
  
  # Data without weights
  dat.sim.noW <- dat.sim
  dat.sim.noW$weights <- NULL
  # Data with voom weights
  dat.sim.voomW <- dat.sim
  dat.sim.voomW$weights_dream <- NULL
  # Data with dream weights
  dat.sim.voomD <- dat.sim
  dat.sim.voomD$weights <- dat.sim.voomD$weights_dream
  dat.sim.voomW$weights_dream <- NULL
  
  # Set to NULL to capture errors
  fdrKnnn <- fdrKnny <- fdrKynn <- fdrKyny <- fdrKyyn <- fdrKyyy <- NULL
  
  # Run models
  #### Paired = N, kinship = N, weights = N ####
  print("NNN")
  ## kimma
  start <- Sys.time()
  fdrKnnn <- kmFit(dat.sim.noW, use.weights = FALSE,
                   model = "~ condition", run.lm = TRUE, 
                   processors=p)[[1]]
  endKnnn <- Sys.time()
  
  #### Paired = N, kinship = N, weights = Y ####
  ## kimma
  fdrKnny <- kmFit(dat.sim.voomW, use.weights = TRUE,
                   model = "~ condition", run.lm = TRUE, 
                   processors=p)[[1]]
  endKnny <- Sys.time()
  
  fdrKnny2 <- kmFit(dat.sim.voomD, use.weights = TRUE,
                   model = "~ condition", run.lm = TRUE, 
                   processors=p)[[1]]
  endKnny2 <- Sys.time()
  
  #### Paired = Y, kinship = N, weights = N ####
  print("YNN")
  ## kimma
  fdrKynn <- kmFit(dat.sim.noW, use.weights = FALSE,
                   model = "~ condition + (1|ptID)", run.lme = TRUE, 
                   processors=p)[[1]]
  endKynn <- Sys.time()
  
  #### Paired = Y, kinship = N, weights = Y ####
  print("YNY")
  ## kimma
  fdrKyny <- kmFit(dat.sim.voomW, use.weights = TRUE,
                   model = "~ condition + (1|ptID)", run.lme = TRUE, 
                   processors=p)[[1]]
  endKyny <- Sys.time()
  
  fdrKyny2 <- kmFit(dat.sim.voomD, use.weights = TRUE,
                   model = "~ condition + (1|ptID)", run.lme = TRUE, 
                   processors=p)[[1]]
  endKyny2 <- Sys.time()
  
  #### kinship = Y ####
  print("XYX")
  ## kimma only
  ### Replace kinship matrix to remove paired design
  # Paired = N, kinship = Y, weights = N
  fdrKnyn <- kmFit(dat.sim.noW, kin=kin.sim.unpair, use.weights = FALSE,
                   patientID = "libID",
                   model = "~ condition + (1|libID)", run.lmerel = TRUE, 
                   processors=p)[[1]]
  endKnyn <- Sys.time()
  
  # Paired = N, kinship = Y, weights = Y
  fdrKnyy <- kmFit(dat.sim.voomW, kin=kin.sim.unpair, use.weights = TRUE,
                   patientID = "libID",
                   model = "~ condition + (1|libID)", run.lmerel = TRUE, 
                   processors=p)[[1]]
  endKnyy <- Sys.time()
  
  fdrKnyy2 <- kmFit(dat.sim.voomD, kin=kin.sim.unpair, use.weights = TRUE,
                   patientID = "libID",
                   model = "~ condition + (1|libID)", run.lmerel = TRUE, 
                   processors=p)[[1]]
  endKnyy2 <- Sys.time()
  
  # Paired = Y, kinship = Y, weights = N
  fdrKyyn <- kmFit(dat.sim.noW, kin=kin.sim, use.weights = FALSE,
                   model = "~ condition + (1|ptID)", run.lmerel = TRUE, 
                   processors=p)[[1]]
  endKyyn <- Sys.time()
  
  # Paired = Y, kinship = Y, weights = Y
  fdrKyyy <- kmFit(dat.sim.voomW, kin=kin.sim, use.weights = TRUE,
                   model = "~ condition + (1|ptID)", run.lmerel = TRUE, 
                   processors=p)[[1]]
  endKyyy <- Sys.time()
  
  fdrKyyy2 <- kmFit(dat.sim.voomD, kin=kin.sim, use.weights = TRUE,
                   model = "~ condition + (1|ptID)", run.lmerel = TRUE, 
                   processors=p)[[1]]
  endKyyy2 <- Sys.time()
  
  #### sensitivity and specificity ####
  #Dont run on multi-processor results. Same as single processor
  # for(m in ls(pattern = "^fdr")){
  #   if(!is.null(get(m))){
  #     ss_resultKM <- ss_calc(i, dat.sim, get(m), m) %>% 
  #       bind_rows(ss_resultKM)
  #   }
  # }
  #### time ####
  time_resultKM <- data.frame(model=c("KnnnM", "KnnyM", "Knny2M", 
                                      "KynnM", "KynyM", "Kyny2M",
                                      "KnynM", "KnyyM", "Knyy2M", 
                                      "KyynM", "KyyyM", "Kyyy2M"),
                              simulation=i,
                              time_s=c(endKnnn-start, endKnny-endKnnn, 
                                       endKnny2-endKnny,
                                       endKynn-endKnny2, endKyny-endKynn, 
                                       endKyny2-endKyny,
                                       endKnyn-endKyny2, endKnyy-endKnyn, 
                                       endKnyy2-endKnyy,
                                       endKyyn-endKnyy2, endKyyy-endKyyn,
                                       endKyyy2-endKyyy)) %>% 
    bind_rows(time_resultKM)
}

save(time_resultKM, #ss_resultKM, 
     file="results/model_fit/simulated_kimma_multi.RData")

#### dream ####
# Cannot be run in foreach loop as dream parallel processing conflicts with dopar parallel
# Restart R session to remove all parallel setup from previous foreach loop

ss_resultD <- data.frame()
time_resultD <- data.frame()
p = 6

for(i in 1:100) {
  #### Data ####
  print(i)
  load(paste0("data_sim/simulated",i,".RData"))
  dat.sim$counts <- NULL
  
  # Data with dream weights
  dat.sim.dreamW <- dat.sim
  dat.sim.dreamW$weights <- dat.sim.dreamW$weights_dream
  dat.sim.dreamW$weights_dream <- NULL
  rownames(dat.sim.dreamW$targets) <- dat.sim.dreamW$targets$libID
  
  # Data no weights
  dat.sim.noW <- dat.sim
  dat.sim.noW$weights <- NULL
  dat.sim.noW$weights_dream <- NULL
  rownames(dat.sim.noW$targets) <- dat.sim.noW$targets$libID
  
  # Set to NULL to capture errors
  fdrDyny <- fitDynyM <- fitDynn <- fitDynn <- NULL
  
  # Run models
  #### Paired = Y, kinship = N, weights = Y ####
  start <- Sys.time()
  fitDyny <- dream(dat.sim.dreamW, ~condition + (1|ptID),
                   useWeights = TRUE,
                   BPPARAM = SnowParam(1),
                   dat.sim.dreamW$targets, REML=TRUE)
  efitDyny <- eBayes(fitDyny)
  fdrDyny <- extract_lmFit(design = fitDyny$design, fit=efitDyny) %>% 
    select(-gene) %>% 
    mutate(gene=symbol)
  
  endDyny <- Sys.time()
  
  fitDynyM <- dream(dat.sim.dreamW, ~condition + (1|ptID),
                   useWeights = TRUE,
                   BPPARAM = SnowParam(p),
                   dat.sim.dreamW$targets, REML=TRUE)
  efitDynyM <- eBayes(fitDynyM)
  fdrDynyM <- extract_lmFit(design = fitDynyM$design, fit=efitDynyM) %>% 
    select(-gene) %>% 
    mutate(gene=symbol)
  
  endDynyM <- Sys.time()
  
  
  #### Paired = Y, kinship = N, weights = N ####
  fitDynn <- dream(dat.sim.noW, ~condition + (1|ptID),
                   useWeights = FALSE,
                   BPPARAM = SnowParam(1),
                   dat.sim.noW$targets, REML=TRUE)
  efitDynn <- eBayes(fitDynn)
  fdrDynn <- extract_lmFit(design = fitDynn$design, fit=efitDynn) %>% 
    select(-gene) %>% 
    mutate(gene=symbol)
  
  endDynn <- Sys.time()
  
  fitDynnM <- dream(dat.sim.noW, ~condition + (1|ptID),
                   useWeights = FALSE,
                   BPPARAM = SnowParam(p),
                   dat.sim.noW$targets, REML=TRUE)
  efitDynnM <- eBayes(fitDynnM)
  fdrDynnM <- extract_lmFit(design = fitDynnM$design, fit=efitDynnM) %>% 
    select(-gene) %>% 
    mutate(gene=symbol)
  
  endDynnM <- Sys.time()
  
  #### sensitivity and specificity ####
  #Dont run on multi-processor results. Same as single processor
  for(m in c("fdrDyny","fdrDynn")){
    if(!is.null(get(m))){
      ss_resultD <- ss_calc(i, dat.sim, get(m), m) %>% 
        bind_rows(ss_resultD)
    }
  }
  
  #### time ####
  time_resultD <- data.frame(model=c("Dyny", "DynyM", "Dynn", "DynnM"),
                             simulation=i,
                             time_s=c(endDyny-start, endDynyM-endDyny, 
                                      endDynn-endDynyM, endDynnM-endDynn)) %>% 
    bind_rows(time_resultD)
}

save(time_resultD, ss_resultD, 
     file="results/model_fit/simulated_dream.RData")

#### DESeq2 ####
p = 6
ss_resultS <- data.frame()
time_resultS <- data.frame()

for(i in 1:100){
  print(i)
  load(paste0("data_sim/simulated",i,".RData"))
  
  dat.sim$targets <- dat.sim$targets %>% 
    mutate(condition = factor(condition),
           ptID = factor(ptID))
  
  start <- Sys.time()
  #### Paired = N, kinship = N, weights = N ####
  print("NNN")
  dat.sim.nnn <- DESeqDataSetFromMatrix(countData = round(dat.sim$counts, digits=0),
                                        colData = dat.sim$targets,
                                        design= ~ condition)
  fitSnnn <- DESeq(dat.sim.nnn)
  fdrSnnn <- as.data.frame(results(fitSnnn)) %>% 
    dplyr::rename(pval=pvalue, FDR=padj) %>% 
    rownames_to_column("gene") %>% 
    mutate(variable = "condition")
  endSnnn <- Sys.time()
  
  dat.sim.nnnM <- DESeqDataSetFromMatrix(countData = round(dat.sim$counts, 
                                                           digits=0),
                                         colData = dat.sim$targets,
                                         design= ~ condition)
  fitSnnnM <- DESeq(dat.sim.nnnM, 
                    parallel =  TRUE, BPPARAM = MulticoreParam(p))
  fdrSnnnM <- as.data.frame(results(fitSnnnM)) %>% 
    dplyr::rename(pval=pvalue, FDR=padj) %>% 
    rownames_to_column("gene") %>% 
    mutate(variable = "condition")
  endSnnnM <- Sys.time()
  
  #### Paired = Y, kinship = N, weights = N ####
  print("YNN")
  dat.sim.ynn <- DESeqDataSetFromMatrix(countData = round(dat.sim$counts, 
                                                          digits=0),
                                        colData = dat.sim$targets,
                                        design= ~ ptID + condition)
  fitSynn <- DESeq(dat.sim.ynn)
  fdrSynn <- as.data.frame(results(fitSynn)) %>% 
    dplyr::rename(pval=pvalue, FDR=padj) %>% 
    rownames_to_column("gene") %>% 
    mutate(variable = "condition")
  endSynn <- Sys.time()
  
  dat.sim.ynnM <- DESeqDataSetFromMatrix(countData = round(dat.sim$counts, 
                                                           digits=0),
                                         colData = dat.sim$targets,
                                         design= ~ ptID + condition)
  fitSynnM <- DESeq(dat.sim.ynnM, 
                    parallel =  TRUE, BPPARAM = MulticoreParam(p))
  fdrSynnM <- as.data.frame(results(fitSynnM)) %>% 
    dplyr::rename(pval=pvalue, FDR=padj) %>% 
    rownames_to_column("gene") %>% 
    mutate(variable = "condition")
  endSynnM <- Sys.time()
  
  #### sensitivity and specificity ####
  #Don't run on multi-processor result. Same as single processor
  for(m in c("fdrSnnn","fdrSynn")){
    if(!is.null(get(m))){
      ss_resultS <- ss_calc(i, dat.sim, get(m), m) %>% 
        bind_rows(ss_resultS)
    }
  }
  
  #### time ####
  time_resultS <- data.frame(model=c("Snnn", "SnnnM",
                                     "Synn", "SynnM"),
                             simulation=i,
                             time_s=c(endSnnn-start, endSnnnM-endSnnn, 
                                      endSynn-endSnnnM, endSynnM-endSynn)) %>% 
    bind_rows(time_resultS)
}

save(time_resultS, ss_resultS, 
     file="results/model_fit/simulated_deseq.RData")

#### Combine and format results ####

load("results/model_fit/simulated/simulated_limma_kimma.RData")
load("results/model_fit/simulated/simulated_kimma_multi.RData")
load("results/model_fit/simulated/simulated_dream.RData")
load("results/model_fit/simulated/simulated_deseq.RData")

time_result <- bind_rows(result[[2]]) %>% 
  bind_rows(bind_rows(time_resultKM)) %>% 
  bind_rows(bind_rows(time_resultD)) %>%
  bind_rows(bind_rows(time_resultS)) %>%
  #format labels
  separate(model, into=c("trash","software","paired",
                         "kinship","weights1","weights_proc"),
           sep="", remove = FALSE, extra = "merge") %>% 
  mutate(software = recode_factor(factor(software),
                                  "K"="kimma", "L"="limma",
                                  "D"="dream", "S"="DESeq2"),
         paired = recode_factor(factor(paired),
                                "n"="unpaired","y"="paired"),
         kinship = recode_factor(factor(kinship),
                                 "n"="no kinship","y"="kinship"),
         weights = recode_factor(factor(weights1),
                                 "n"="no weights","y"="weights"),
         weights_type = paste0(weights1, weights_proc),
         weights_type = ifelse(software == "dream" & weights != "no weights",
                               "dream\nweights",
                               weights_type),
         weights_type = recode_factor(factor(weights_type),
                                 "nNA"="no\nweights",
                                 "nM"="no\nweights",
                                 "yM"="voom\nweights",
                                 "yNA"="voom\nweights",
                                 "y2M"="dream\nweights",
                                 "y2"="dream\nweights"),
         processors = ifelse(is.na(weights_proc) | weights_proc == 2,
                             "1 processor", "6 processors"),
         processors = factor(processors, c("1 processor","6 processors"))) %>% 
  select(model, kinship, paired, weights, weights_type, software,
         processors, simulation, time_s)

ss_result <- bind_rows(result[[1]]) %>% 
  bind_rows(bind_rows(ss_resultD)) %>%
  bind_rows(bind_rows(ss_resultS)) %>%
  #format labels
  mutate(model = gsub("fdr","", model)) %>% 
  separate(model, into=c("trash","software","paired",
                         "kinship","weights1","weights_proc"),
           sep="", remove = FALSE) %>% 
  mutate(software = recode_factor(factor(software),
                                  "K"="kimma", "L"="limma",
                                  "D"="dream", "S"="DESeq2"),
         paired = recode_factor(factor(paired),
                                "n"="unpaired","y"="paired"),
         kinship = recode_factor(factor(kinship),
                                 "n"="no kinship","y"="kinship"),
         weights = recode_factor(factor(weights1),
                                 "n"="no weights","y"="weights"),
         weights_type = paste0(weights1, weights_proc),
         weights_type = ifelse(software == "dream" & weights != "no weights", 
                               "dream\nweights",
                               weights_type),
         weights_type = recode_factor(factor(weights_type),
                                      "nNA"="no\nweights",
                                      "nM"="no\nweights",
                                      "yM"="voom\nweights",
                                      "yNA"="voom\nweights",
                                      "y2M"="dream\nweights",
                                      "y2"="dream\nweights"),
         processors = ifelse(is.na(weights_proc) | weights_proc == 2, 
                             "1 processor", "6 processors"),
         processors = factor(processors, c("1 processor","6 processors"))) %>% 
  select(model, kinship, paired, weights, weights_type, software,
         processors, simulation, cutoff:genes_fdr)

save(ss_result, time_result, 
     file="results/simulated_fit.RData")
