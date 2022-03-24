library(tidyverse)
attach("data/RSTR_kimma_dat.RData")
attach("data/RSTR_kimma_voom.RData")
attach("data/RSTR_kimma_dream.RData")
attach("data/RSTR_RNAseq_kinship.RData")

#### Media only samples ####
dat.m <- dat
dat.m$samples <- dat.m$samples %>% 
  filter(condition=="MEDIA")
dat.m$counts <- dat.m$counts %>% 
  select(all_of(dat.m$samples$libID))

datV.m <- datV
datV.m$targets <- datV.m$targets %>% 
  filter(condition=="MEDIA")
datV.m$E <- datV.m$E %>% 
  select(all_of(datV.m$targets$libID))
datV.m$weights <- as.data.frame(datV.m$weights) %>% 
  select(all_of(datV.m$targets$libID))

datD.m <- datD
datD.m$targets <- datD.m$targets %>% 
  filter(condition=="MEDIA")
datD.m$E <- datD.m$E %>% 
  select(all_of(datD.m$targets$libID))
datD.m$weights <- as.data.frame(datD.m$weights) %>% 
  select(all_of(datD.m$targets$libID))

#### Simulated samples base ####
N <- nrow(dat.m$samples)

metaS <- data.frame(libID = paste0("lib", 1:(N*2)),
                          ptID = paste0("pt", 1:N),
                          condition = c(rep("baseline", N),
                                        rep("simulated", N)))
## Repeat weights for simulated data
weightsVS <- cbind(datV.m$weights, datV.m$weights)
weightsDS <- cbind(datD.m$weights, datD.m$weights)

# Make 1 data object
datS <- datV
datS$counts <- dat.m$counts
datS$E <- datV.m$E
datS$targets <- metaS
datS$genes <- datV.m$genes
datS$weights <- as.matrix(weightsVS)
datS$weights_dream <- as.matrix(weightsDS)

## De-identify
colnames(datS$counts) <- paste0("lib", 1:N)
colnames(datS$E) <- paste0("lib", 1:N)
colnames(datS$weights) <- paste0("lib", 1:(N*2))
colnames(datS$weights_dream) <- paste0("lib", 1:(N*2))

#### Format kinship to match ####
meta.sub <- datV$targets %>% 
  filter(condition=="MEDIA") %>% 
  distinct(ptID, libID)
kin.sim <- as.data.frame(kin.matrix) %>% 
  rownames_to_column("ptID") %>% 
  left_join(meta.sub) %>% 
  arrange(libID) %>% 
  mutate(rowname = paste0("pt", 1:length(libID))) %>% 
  select(rowname, everything(), -libID, -ptID) %>% 
  column_to_rownames() %>% t() %>% as.data.frame() %>% 
  rownames_to_column("ptID") %>% 
  left_join(meta.sub) %>% 
  arrange(libID) %>% 
  mutate(rowname = paste0("pt", 1:length(libID))) %>% 
  select(rowname, everything(), -libID, -ptID) %>% 
  column_to_rownames() %>% t()

#create unpaired version
#Create fake kinship
set.seed(33)
kin.sim.unpair <- matrix(nrow=N*2, ncol=N*2,
                         runif((N*2)^2, min=0, max = 0.01))
#Copy lower triangle to upper
kin.sim.unpair[upper.tri(kin.sim.unpair)] <- NA
kin.sim.unpair[upper.tri(kin.sim.unpair)] <- t(kin.sim.unpair)[upper.tri(t(kin.sim.unpair))]
#Fill in diagonal
diag(kin.sim.unpair) <- 1
#Name
rownames(kin.sim.unpair) <- paste0("lib", 1:(N*2))
colnames(kin.sim.unpair) <- paste0("lib", 1:(N*2))

save(kin.sim, kin.sim.unpair, file="data_sim/simulated_kinship.RData")

#### Create simulated data ####
set.seed(33)
for(i in 1:100){
  # Subset to random 1000 genes
  dat.sub <- datS
  index <- sample(1:nrow(dat$genes), 1000, replace=FALSE)
  dat.sub$genes <- dat.sub$genes[index, ]
  dat.sub$E <- dat.sub$E[index, ] 
  dat.sub$counts <- dat.sub$counts[index, ] 
  dat.sub$weights <- dat.sub$weights[index, ]
  dat.sub$weights_dream <- dat.sub$weights_dream[index, ]
  
  # Introduce 50 genes each with fold change 0.1, 0.25, 0.5, 1, 2
  ## raw counts
  
  dat.sim.01 <- dat.sub$counts[c(1:50), ]*1.1
  dat.sim.05 <- dat.sub$counts[c(51:100), ]*1.5
  dat.sim.1 <- dat.sub$counts[c(101:150), ]*2
  dat.sim.005 <- dat.sub$counts[c(151:200), ]*1.05
  dat.sim.001 <- dat.sub$counts[c(201:250), ]*1.01
  dat.sim.base <- dat.sub$counts[c(251:1000), ]
  
  # Combine and rename with new libID
  dat.sim <- rbind(dat.sim.01, dat.sim.05, dat.sim.1,
                    dat.sim.005, dat.sim.001, dat.sim.base)
  colnames(dat.sim) <- paste0("lib", (N+1):(N*2))
  
  # Introduce simulated error by multiplying by a random number between 0.95 and 1.05
  rand.mat <- matrix(nrow=nrow(dat.sim), ncol=ncol(dat.sim),
                     runif(nrow(dat.sim)*ncol(dat.sim), min=0.95, max=1.05))
  dat.sim.er <- dat.sim*rand.mat
  
  # Combine original and simulated expression
  dat.sub$counts <- cbind(dat.sub$counts, dat.sim.er)
  
  ## Voom data
  dat.sim.01b <- dat.sub$E[c(1:50), ]*1.1
  dat.sim.05b <- dat.sub$E[c(51:100), ]*1.5
  dat.sim.1b <- dat.sub$E[c(101:150), ]*2
  dat.sim.005b <- dat.sub$E[c(151:200), ]*1.05
  dat.sim.001b <- dat.sub$E[c(201:250), ]*1.01
  dat.sim.baseb <- dat.sub$E[c(251:1000), ]
  
  # Combine and rename with new libID
  dat.simb <- rbind(dat.sim.01b, dat.sim.05b, dat.sim.1b,
                    dat.sim.005b, dat.sim.001b, dat.sim.baseb)
  colnames(dat.simb) <- paste0("lib", (N+1):(N*2))
  
  # Introduce simulated error by multiplying by a random number between 0.95 and 1.05
  dat.sim.erb <- dat.simb*rand.mat
  
  # Combine original and simulated expression
  dat.sub$E <- cbind(dat.sub$E, dat.sim.erb)
  
  # Add simulated data group to metadata
  dat.sub$genes <- dat.sub$genes %>% 
    mutate(sim.group = c(rep("FC0.1", 50), rep("FC0.25", 50),
                         rep("FC0.5", 50), rep("FC1", 50),
                         rep("FC2", 50), rep("none", 750)))

  dat.sim <- dat.sub
  # Save output
  save(dat.sim, file=paste0("data_sim/simulated",i,".RData"))
}
