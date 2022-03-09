library(tidyverse)
library(venn)

attach("results/model_fit/Sample_Group_kimma_kinship.RData")
fit_kimma_kin <- fit
fdr_kimma_kin <- results_kimma_LTBI_vs_RSTR_kinship
attach("results/model_fit/Sample_Group_kimma_kinship_weights.RData")
fit_kimma_kin_w <- fit
fdr_kimma_kin_w <- results_kimma_LTBI_vs_RSTR_kinship_weights

#Venn
# other packages: venneuler, VennDiagram

venn.ls1 <- list()
venn.ls1[["w/o weights"]] <- fdr_kimma_kin %>% 
  filter(FDR < 0.1 & variable == "Sample_GroupRSTR") %>% 
  distinct(gene) %>% unlist()
venn.ls1[["w/ weights"]] <- fdr_kimma_kin_w %>% 
  filter(FDR < 0.1 & variable == "Sample_GroupRSTR") %>% 
  distinct(gene) %>% unlist()
venn(venn.ls1)

# Compare P-value and FDR
fdr1 <- fdr_kimma_kin %>% 
  filter(variable == "Sample_GroupRSTR") %>% 
  select(gene, variable, pval, FDR) %>% 
  mutate(group = "w/o weights")
fdr2 <- fdr_kimma_kin_w %>% 
  filter(variable == "Sample_GroupRSTR") %>% 
  select(gene, variable, pval, FDR) %>% 
  mutate(group = "w/ weights")

fdr_bind <- bind_rows(fdr1, fdr2) %>% 
  pivot_longer(pval:FDR) %>% 
  mutate(newname = paste(group, name)) %>% 
  select(newname, value, gene, variable) %>% 
  pivot_wider(names_from = newname)



####

attach("results/model_fit/LTBI_limma_paired_no_weights.RData")
fit_limma <- paired.efit
fdr_limma <- LTBI.paired.results.limma
attach("results/model_fit/LTBI_limma_paired_weights.RData")
fit_limma_w <- paired.efit.w
fdr_limma_w <- LTBI.paired.results.limma.W

#Venn
venn.ls2 <- list()
venn.ls2[["w/o weights"]] <- fdr_limma %>% 
  filter(adj.P.Val < 0.01 & variable == "conditionTB") %>% 
  distinct(geneName) %>% unlist()
venn.ls2[["w/ weights"]] <- fdr_limma_w %>% 
  filter(adj.P.Val < 0.01 & variable == "conditionTB") %>% 
  distinct(geneName) %>% unlist()
venn(venn.ls2)

# Compare P-value and FDR
fdr1 <- fdr_limma %>% 
  filter(variable == "conditionTB") %>% 
  select(geneName, variable, P.Value, adj.P.Val) %>% 
  mutate(group = "w/o weights")
fdr2 <- fdr_limma_w %>% 
  filter(variable == "conditionTB") %>% 
  select(geneName, variable, P.Value, adj.P.Val) %>% 
  mutate(group = "w/ weights")

fdr_bind <- bind_rows(fdr1, fdr2) %>% 
  pivot_longer(P.Value:adj.P.Val) %>% 
  mutate(newname = paste(group, name)) %>% 
  select(newname, value, geneName, variable) %>% 
  pivot_wider(names_from = newname)

# Combining venns
par(mfrow=c(1,2))
venn(venn.ls1)
title(main="kimma with kinship", sub = "FDR < 0.1", line = -5)
venn(venn.ls2)
title(main="limma paired", sub = "FDR < 0.01", line = -5)
dev.off()

# Upset plots
# alternative to venn
library(UpSetR)

upset(fromList(venn.ls1))
upset(fromList(venn.ls2))
