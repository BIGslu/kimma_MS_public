library(tidyverse)
library(ggvenn)
library(SEARchways)

#### Data ####
load("results/rstr_interaction_fit.RData")

deg_rel <- related.lmerel$lmerel %>%
  filter(variable == "condition:Sample_Group" & FDR < 0.2) %>% 
  pull(gene)

deg_unrel <- unrelated.lmerel$lmerel %>%
  filter(variable == "condition:Sample_Group" & FDR < 0.2) %>% 
  pull(gene)

#### Overlap DEG with Mtb-induced DEG ####

#RSTR-Mtb paper 
deg_all <- read_csv("https://raw.githubusercontent.com/hawn-lab/RSTR_RNAseq_Mtb_public/main/Uganda/results/gene_level/RSTR.Mtb.model.results.anno.csv") %>% 
  filter(variable == "conditionTB:Sample_GroupRSTR" & FDR < 0.2) %>% 
  pull(gene)

ggvenn(list("unrelated"=deg_unrel, "related"=deg_rel, 
            "all"=deg_all), show_percentage = FALSE,
       text_size = 3, set_name_size = 3, stroke_size = 0.5)

#### Enrichment ####
H <- BIGprofiler(gene_list=list("related"=deg_rel), category = "H")

#### Enrichment plot ####

plot2 <- H %>% 
  filter(FDR < 0.1) %>% 
  mutate(pathway = gsub("HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway)) %>% 
  mutate(Significance = ifelse(FDR <= 0.1, paste("FDR <", 0.1), "NS")) %>%
  
  ggplot(aes(x=reorder(pathway,`k/K`), y=`k/K`, fill=Significance)) +
  geom_col() +
  scale_fill_manual(values=c("FDR < 0.1"="#fdae61")) +
  coord_flip() +
  theme_classic() +
  labs(x="", y="Proportion enriched (k/K)")
plot2
ggsave("figs/RSTR_DEG_enrich.png", plot2, width=6, height=2)

#### GSEA ####
# recalc FDR within interaction DEG
related.contrast <- related.lmerel$lmerel.contrast %>% 
  filter(variable == "condition*Sample_Group" & gene %in% deg_rel) %>% 
  group_by(variable, contrast_ref, contrast_lvl) %>% 
  mutate(FDR.new = p.adjust(pval, method="BH")) %>% 
  ungroup()

FC.ls <- list()

FC.ls[["RSTRinMEDIA"]] <- related.contrast %>% 
  filter(contrast_ref == "MEDIA LTBI" & contrast_lvl == "MEDIA RSTR") %>% 
  pull(estimate, name=gene)

FC.ls[["RSTRinTB"]] <- related.contrast %>% 
  filter(contrast_ref == "TB LTBI" & contrast_lvl == "TB RSTR") %>% 
  pull(estimate, name=gene)

FC.ls[["TBinLTBI"]] <- related.contrast %>% 
  filter(contrast_ref == "MEDIA LTBI" & contrast_lvl == "TB LTBI") %>% 
  pull(estimate, name=gene)

FC.ls[["TBinRSTR"]] <- related.contrast %>% 
  filter(contrast_ref == "MEDIA RSTR" & contrast_lvl == "TB RSTR") %>% 
  pull(estimate, name=gene)

H.gsea <- BIGsea(gene_list = FC.ls, category = "H")

#### GSEA plot ####
#Mtb-induced gsea + signif media gsea

H.gsea.signifTB <- H.gsea %>% 
  filter(FDR < 0.1 & group %in% c("TBinLTBI","TBinRSTR")) %>% 
  pull(pathway)
H.gsea.signifRSTR <- H.gsea %>% 
  filter(FDR < 0.1 & group %in% c("RSTRinMEDIA","RSTRinTB")) %>% 
  pull(pathway)
H.gsea.signif <- intersect(H.gsea.signifTB, H.gsea.signifRSTR)

load(url("https://github.com/hawn-lab/RSTR_RNAseq_Mtb_public/blob/main/publication/GSEA.signif.RData?raw=true"))

H.gsea.signif.all <- unique(c(H.gsea.signif, h.signif$pathway))

#Format and order signif terms 
H.gsea.format <- H.gsea %>% 
  #color by significance
  mutate(Significance = ifelse(FDR <= 0.1, paste("FDR <", 0.1), "NS")) %>%
  filter(pathway %in% H.gsea.signif.all) %>% 
  #Beautify labels
  mutate(pathway = gsub("HALLMARK_","",pathway)) %>% 
  ##lowercase and then correction
  mutate(pathway = gsub("_"," ", tolower(pathway)),
         pathway = gsub("tnfa","TNFA",pathway),
         pathway = gsub("nfkb","NF-kB",pathway),
         pathway = gsub("interferon","IFN",pathway),
         pathway = gsub("il2","IL2",pathway),
         pathway = gsub("stat5","STAT5",pathway),
         pathway = gsub("kras","KRAS",pathway),
         pathway = gsub(" dn"," down",pathway),
         pathway = gsub("myc ","MYC ",pathway),
         pathway = gsub("mtorc1","MTORC1",pathway)) %>% 
  ## order
  mutate(pathway = factor(pathway,
                          levels=c(
                            #New in related sub-analysis
                            "apoptosis",
                            "complement",
                            "KRAS signaling up",

                            #Split to up R
                            "coagulation",

                            #Split R
                            "epithelial mesenchymal transition",
                            "hypoxia",
                            "KRAS signaling down",
                            "IL2 STAT5 signaling",
                            #All R down
                            "IFN alpha response",
                            "allograft rejection",
                            "IFN gamma response",

                            #Down with TB
                            "adipogenesis",
                            "oxidative phosphorylation",
                            #Down to up R
                            "inflammatory response",
                            "TNFA signaling via NF-kB"))) %>%
  mutate(group2 = ifelse(pathway %in% c("adipogenesis",
                                        "oxidative phosphorylation"), "ii",
                         ifelse(pathway %in% c("IFN alpha response",
                                               "allograft rejection",
                                               "IFN gamma response"), "iii",
                                ifelse(pathway %in% c("inflammatory response",
                                                      "TNFA signaling via NF-kB"), "i",
                                       ifelse(pathway %in% c("IL2 STAT5 signaling",
                                                             "KRAS signaling down",
                                                             "hypoxia","coagulation",
                                                             "epithelial mesenchymal transition"), "iv", "v")))))

plot.lim <- max(abs(H.gsea.format$NES))+0.15

#facet labels
facet_lab <- c('RSTRinMEDIA'="Down in RSTR <-  -> Up in RSTR",
               'RSTRinTB'="Down in RSTR <-  -> Up in RSTR",
               'TBinLTBI'="Down in +Mtb <-  -> Up in +Mtb",
               'TBinRSTR'="Down in +Mtb <-  -> Up in +Mtb",
               'i'='i','ii'='ii','iii'='iii','iv'='iv', 'v'='v')

plot1 <- H.gsea.format %>% 
  filter(group %in% c("RSTRinMEDIA","RSTRinTB")) %>%
  ggplot() +
  geom_segment(aes(x=pathway, xend=pathway, 
                   y=0, yend=NES),
               size=0.5) +
  geom_point(aes(x=pathway, y=NES, fill = Significance),
             stroke=0.7, shape=24) +
  geom_hline(yintercept = 0) +
  
  scale_fill_manual(values=c("FDR < 0.1"="#fdae61",
                             "NS"="grey")) +
  lims(y=c(-plot.lim,plot.lim)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized enrichment score (NES)",
       fill = "FGSEA significance", shape="",
       title="a) RSTR vs LTBI in media   b) RSTR vs LTBI in +Mtb") + 
  facet_grid(group2 ~ group, scales="free_y", 
             labeller = as_labeller(facet_lab),
             space = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "vertical",
        strip.background = element_rect(fill="white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  #force legend color/shape
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "black")),
         size = "none") 

plot1
ggsave("figs/RSTR_DEG_GSEA.png", plot1, width=6, height=5)
