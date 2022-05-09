library(tidyverse)
library(kimma)

#### Data ####
attach("data/RSTR_data_all.RData")

#Subset to MEDIA vs TB in LTBI samples
datV.ltbi <- datV
datV.ltbi$targets <- datV.ltbi$targets %>% 
  filter(Sample_Group == "LTBI") %>% 
  #fake sex group
  mutate(sex = sample(c("M","F"), replace=TRUE, size=98))
datV.ltbi$E <- as.data.frame(datV.ltbi$E) %>% 
  select(all_of(datV.ltbi$targets$libID))
# identical(colnames(datV.ltbi$E), datV.ltbi$targets$libID)
datV.ltbi$weights <- as.data.frame(datV.ltbi$weights) %>% 
  select(all_of(datV.ltbi$targets$libID)) %>% 
  as.matrix()

#### models ####
time_trial <- Sys.time()
temp <- kmFit(dat = datV.ltbi,
              model = "~ condition", run.lm = TRUE, 
              use.weights = FALSE, metrics = FALSE,
              processors=6)
time_trial <- c(time_trial, Sys.time())
#metrics
temp <- kmFit(dat = datV.ltbi,
              model = "~ condition", run.lm = TRUE, 
              use.weights = FALSE, metrics = TRUE,
              processors=6)
time_trial <- c(time_trial, Sys.time())
#weights
temp <- kmFit(dat = datV.ltbi,
              model = "~ condition", run.lm = TRUE, 
              use.weights = TRUE, metrics = FALSE,
              processors=6)
time_trial <- c(time_trial, Sys.time())
#paired
temp <- kmFit(dat = datV.ltbi,
              model = "~ condition+(1|ptID)", run.lme = TRUE, 
              use.weights = FALSE, metrics = FALSE,
              processors=6)
time_trial <- c(time_trial, Sys.time())
#covar
temp <- kmFit(dat = datV.ltbi,
              model = "~ condition+sex", run.lm = TRUE, 
              use.weights = FALSE, metrics = FALSE,
              processors=6)
time_trial <- c(time_trial, Sys.time())
#interaction
temp <- kmFit(dat = datV.ltbi,
              model = "~ condition*sex", run.lm = TRUE, 
              use.weights = FALSE, metrics = FALSE,
              processors=6)
time_trial <- c(time_trial, Sys.time())
#contrasts
temp <- kmFit(dat = datV.ltbi,
              model = "~ condition*sex", run.lm = TRUE, 
              run.contrast=TRUE, contrast.var="condition:sex",
              use.weights = FALSE, metrics = FALSE,
              processors=6)
time_trial <- c(time_trial, Sys.time())
#kin
temp <- kmFit(dat = datV.ltbi, kin=kin,
              model = "~ condition+(1|ptID)", run.lmerel = TRUE, 
              use.weights = FALSE, metrics = FALSE,
              processors=6)
time_trial <- c(time_trial, Sys.time())
#all
temp <- kmFit(dat = datV.ltbi, kin=kin,
              model = "~ condition*sex+(1|ptID)", 
              run.lmerel = TRUE, 
              run.contrast=TRUE, contrast.var="condition:sex",
              use.weights = TRUE, metrics = TRUE,
              processors=6)
time_trial <- c(time_trial, Sys.time())
#multi-models
temp <- kmFit(dat = datV.ltbi, kin=kin,
              model = "~ condition*sex+(1|ptID)", 
              run.lm = TRUE, run.lme = TRUE, run.lmerel = TRUE, 
              run.contrast=TRUE, contrast.var="condition:sex",
              use.weights = TRUE, metrics = TRUE,
              processors=6)
time_trial <- c(time_trial, Sys.time())

#### combine results ####
time_trial_df <- data.frame(
  model = c("none","base","metrics","weights","paired","covariate","interaction",
            "contrast","kin","all","multi"),
  time = time_trial) %>% 
  mutate(diff = time-lag(time),
         diff_from_base = diff - diff[model=="base"])  

save(time_trial_df, file="results/kimma_time_trial.RData")
