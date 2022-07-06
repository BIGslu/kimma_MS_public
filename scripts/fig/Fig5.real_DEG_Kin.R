library(tidyverse)
library(ggvenn)
library(patchwork)

#colors
clr_vec <- c("no kinship"="#FFB000","none"="#DC267F",
             "kinship"="#648FFF")

#### Data ####
load("results/condition_subset_fit.RData")
load("results/rstr_subset_fit.RData")

#### Mtb infection ####
signif.ls1 <- list()

signif.ls1[["no kinship"]] <- condition_subset_result %>% 
  filter(kinship == "no kinship" & paired == "paired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "unrelated") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls1[["kinship"]] <- condition_subset_result %>% 
  filter(kinship == "kinship" & paired == "paired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "unrelated") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p1 <- ggvenn(signif.ls1, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5) +
  labs(title="A  Unrelated", y = "Mtb-infected vs media in LTBI\nPaired") +
  theme(axis.title.y = element_text(angle=90), text = element_text(size=8))
# p1
#
signif.ls2 <- list()

signif.ls2[["no kinship"]] <- condition_subset_result %>% 
  filter(kinship == "no kinship" & paired == "paired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls2[["kinship"]] <- condition_subset_result %>% 
  filter(kinship == "kinship" & paired == "paired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p2 <- ggvenn(signif.ls2, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5) +
  labs(title="B  Related") +
  theme(text = element_text(size=8))
# p2

#### RSTR v LTBI ####
signif.ls3 <- list()

signif.ls3[["no kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "no kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "weights" & subset == "unrelated") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls3[["kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "weights" & subset == "unrelated") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p3 <- ggvenn(signif.ls3, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5) +
  labs(title="C  Unrelated", y="RSTR vs LTBI in Mtb-infected\nUnpaired") +
  theme(axis.title.y = element_text(angle=90), text = element_text(size=8))
# p3

signif.ls4 <- list()

signif.ls4[["no kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "no kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls4[["kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p4 <- ggvenn(signif.ls4, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5) +
  labs(title="D  Related") +
  theme(text = element_text(size=8))
# p4

#### eQTL ####
# attach("data/RSTR_data_all.RData")
load("data/RSTR_snp.RData")
load("data/related.unrelated.RData")

dat.eqtl <- eqtl$lmerel %>% 
  filter(variable != '(1 | ptID)') %>% 
  mutate(FDR = p.adjust(pval, method = "BH"))

gene.OI <- dat.eqtl %>% 
  slice_min(pval)

temp <- datV$E %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% gene.OI$gene) %>% 
  pivot_longer(-gene, names_to = "libID", values_to = "E") %>% 
  left_join(datV$targets) %>% 
  filter(condition == "TB")

temp2 <- geno.final %>% 
  filter(refsnp_id %in% gene.OI$genotype) %>% 
  pivot_longer(-c(refsnp_id:gene), names_to = "ptID", values_to = "genotype") %>% 
  left_join(temp)

p5 <- temp2 %>% 
  mutate(group = ifelse(ptID %in% related,"Related","Unrelated"),
         facet = paste("P =", formatC(gene.OI$pval,digits=2,format="E"))) %>% 

  ggplot(aes(x=as.character(genotype), y=E))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=group, color=group), width=0.2, height=0) +
  theme_classic() +
  labs(x=gene.OI$genotype, title="E", 
       y=paste("Normalized log2 expression", gene.OI$gene, sep="\n")) +
  geom_smooth(method = "lm", se=FALSE, aes(group=1), color="red") +
  theme(legend.position = "top", legend.direction = "vertical",
        text = element_text(size=8)) +
  facet_wrap(~facet) +
  scale_colour_manual(name = "",
                      labels = c("Related","Unrelated"),
                      values = c("grey","black")) +   
  scale_shape_manual(name = "",
                     labels = c("Related","Unrelated"),
                     values = c(19, 17))
# p5

#### Save plots #####

p <- p1+p2+p3+p4
# p


ggsave("figs/Fig5.real_DEG_Kin.png", p, height=5, width = 5)
# ggsave("figs/Fig5.real_DEG_Kin.pdf", p, height=5, width = 5)

ggsave("figs/Fig5.eQTL.png", p5, height=5, width = 2)
# ggsave("figs/Fig5.eQTL.pdf", p, height=5, width = 2)
