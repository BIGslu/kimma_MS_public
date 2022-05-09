library(tidyverse)
library(kimma)
library(BIGpicture)
library(SEARchways)

#### Data ####
attach("data/RSTR_data_all.RData")

#Subset to MEDIA vs TB in LTBI samples
datV.ltbi <- datV
datV.ltbi$targets <- datV.ltbi$targets %>% 
  filter(Sample_Group == "LTBI")
datV.ltbi$E <- as.data.frame(datV.ltbi$E) %>% 
  select(all_of(datV.ltbi$targets$libID))
# identical(colnames(datV.ltbi$E), datV.ltbi$targets$libID)
datV.ltbi$weights <- as.data.frame(datV.ltbi$weights) %>% 
  select(all_of(datV.ltbi$targets$libID)) %>% 
  as.matrix()

#### Model ####
result <- kmFit(dat = datV.ltbi, kin=kin,
                model = "~ condition + (1|ptID)", 
                run.lme = TRUE, run.lme = TRUE, run.lmerel = TRUE, 
                use.weights = FALSE, metrics = TRUE)
summarise_kmFit(result$lme)

#### AIC ####
aic <- plot_fit(result, result, "lm","lme", metrics = "AIC")
aic
ggsave("figs/intro_diagram/AIC.pdf", aic, width=2.5, height=3)

#### PCA ####
pca <- plot_pca(datV.ltbi, vars = "condition")[[1]] +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
pca
ggsave("figs/intro_diagram/PCA.pdf", pca, width=3, height=2)

#### Venn ####
vd <- plot_venn_genes(result,"lme", fdr.cutoff = 0.05)
vd
ggsave("figs/intro_diagram/venn.pdf", vd, width=3, height=3.5)

#### STRING ####
genes.OI <- result$lme %>% 
  slice_min(order_by = FDR, n=10) %>% 
  pull(gene)

#Enrichment
enrich <- BIGprofiler(list(top10=genes.OI), category = "H")
  
#network
map <- map_string(genes.OI, score_threshold = 700)
nw <- plot_string(map, enrichment = enrich, fdr.cutoff = 0.01)
nw

ggsave("figs/intro_diagram/string.pdf", nw, width=8, height=4)

#### GSEA ####
fc <- result$lme %>% 
  filter(variable=="condition") %>% 
  pull(estimate, name=gene)

gs <- BIGsea(list(condition=fc), category = "H")

gs_plot <- gs %>% 
  filter(pathway %in% c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                        "HALLMARK_ALLOGRAFT_REJECTION",
                        "HALLMARK_APOPTOSIS")) %>% 
  mutate(pathway = gsub("HALLMARK_","",pathway),
         pathway = gsub("_"," ", pathway)) %>% 
  plot_gsea()
gs_plot

ggsave("figs/intro_diagram/gsea.pdf", gs_plot, width=5, height=1.5)
