library(tidyverse)
library(ggvenn)
library(patchwork)

#colors
clr_vec <- c("no kinship"="#FFB000","none"="#DC267F",
             "kinship"="#648FFF")

#### Data ####
load("results/condition_subset_fit.RData")

#### AIC ####
aic <- condition_subset_metric %>%
  select(subset, weights, kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) %>% 
  mutate(subset = str_to_title(subset),
         subset = factor(subset, levels=c("Unrelated","Related"))) %>% 
  #final models
  filter(weights == "voom weights")  %>% 
  #Change in AIC
  mutate(aic.delta = `no kinship`-kinship) 
  

#### AIC ####
p1 <- aic %>% 
  
  ggplot(aes(x=subset, y=aic.delta)) +
  geom_violin(fill="grey90") +
  geom_hline(yintercept = c(-2,2), lty="dashed", color="grey60") + 
  stat_summary(fun.y=mean, geom="point", size=2, color="black", shape="square") +
  theme_classic() +
  labs(y = "Per gene change in AIC\nBetter fit without kinship <--   --> Better fit with kinship       ", 
       x="") 

# p1

#### Venn ####

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

p2 <- ggvenn(signif.ls1, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5)+
  annotate(geom = "text", label = "Unrelated", x = 0, y = -1.3, size=3)
# p2

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

p3 <- ggvenn(signif.ls2, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5) +
  annotate(geom = "text", label = "Related", x = 0, y = -1.3, size=3)
# p3

#### eQTL ####
attach("data/RSTR_data_all.RData")
load("data/RSTR_snp.RData")
load("data/related.unrelated.RData")
load("results/rstr_subset_eqtl.RData")

dat.eqtl <- eqtl$lmerel %>% 
  filter(grepl(":", variable))

gene.OI <- dat.eqtl %>% 
  slice_min(FDR)

#Get allele info
library(SNPRelate)
gds <- snpgdsOpen("data/merged_Omni_MegaEx_randomFounder_Hawn_filter_LD.gds")
alleles <- data.frame(
  snpID = read.gdsn(index.gdsn(gds, "snp.id")),
  allele = read.gdsn(index.gdsn(gds, "snp.allele")),
  chr = read.gdsn(index.gdsn(gds, "snp.chromosome")),
  pos = read.gdsn(index.gdsn(gds, "snp.position"))
) %>% 
  filter(snpID == gene.OI$genotype) %>% 
  separate(allele, into=c("A","B"))
#0 is two B allele 
# https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html

#Expression data
temp <- datV$E %>% 
  rownames_to_column("gene") %>% 
  filter(gene %in% gene.OI$gene) %>% 
  pivot_longer(-gene, names_to = "libID", values_to = "E") %>% 
  left_join(datV$targets) %>% 
  filter(Sample_Group == "LTBI" & ptID %in% related)

temp2 <- geno.final %>% 
  filter(refsnp_id %in% gene.OI$genotype) %>% 
  pivot_longer(-c(refsnp_id:gene), names_to = "ptID", values_to = "genotype") %>% 
  inner_join(temp)

p4 <- temp2 %>% 
  mutate(group = ifelse(ptID %in% related,"Related","Unrelated"),
         facet = paste("P =", formatC(gene.OI$pval,digits=2,format="E")),
         condition = recode(condition, "MEDIA"="Media", "TB"="Mtb"),
         allele = case_when(genotype == 0 ~ paste0(alleles$B, alleles$B),
                            genotype == 1 ~ paste0(alleles$B, alleles$A),
                            genotype == 2 ~ paste0(alleles$A, alleles$A))) %>% 
  
  ggplot(aes(x=reorder(allele, genotype), y=E))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, height=0) +
  theme_classic() +
  labs(x=gene.OI$genotype, 
       y=paste("Normalized log2 expression", gene.OI$gene, sep="\n")) +
  geom_smooth(method = "lm", se=FALSE, aes(group=condition),
              color="grey60") +
  
  theme(text = element_text(size=8)) +
  facet_wrap(~condition, ncol=2, scales="free_x")


#### Combine ####
design <- "
AB
AC
DD
"
p_final <- p1+p2+p3+p4 +
  plot_annotation(tag_levels = "A",
                  title = "Paired\nMtb vs media in LTBI") +
  plot_layout(design = design) &
  theme(plot.tag = element_text(size=13))
# p_final

#### Save ####
ggsave("figs/Fig5.real_DEG_AIC_mtb.png", p_final, height=7, width = 4.5)
ggsave("figs/Fig5.real_DEG_AIC_mtb.pdf", p_final, height=7, width = 4.5)


#### Check fit new DEG with kinship ####
outersect <- function(x, y) {
  unique(c(setdiff(x, y),
         setdiff(y, x)))
}

outside1 <- setdiff(signif.ls1[[1]], signif.ls1[[2]])
outside2 <- setdiff(signif.ls1[[2]], signif.ls1[[1]])

condition_subset_metric %>% 
  filter(subset == "unrelated" & weights =="voom weights" &
           gene %in% c(outside1,outside2)) %>% 
  select(kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) %>% 
  mutate(diff = kinship-`no kinship`,
         best = ifelse(diff > 2, "no kinship",
                       ifelse(diff < -2, "kinship", "none")),
         group = ifelse(gene %in% outside1, "no kinship", "kinship")) %>% 
  count(group, best)

outside1 <- setdiff(signif.ls2[[1]], signif.ls2[[2]])
outside2 <- setdiff(signif.ls2[[2]], signif.ls2[[1]])

condition_subset_metric %>% 
  filter(subset == "unrelated" & weights =="voom weights" &
           gene %in% c(outside1,outside2)) %>% 
  select(kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) %>% 
  mutate(diff = kinship-`no kinship`,
         best = ifelse(diff > 0, "no kinship",
                       ifelse(diff < 0, "kinship", "none")),
         group = ifelse(gene %in% outside1, "no kinship", "kinship")) %>% 
  count(group, best)
