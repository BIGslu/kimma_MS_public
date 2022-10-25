library(tidyverse)

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

p5 <- temp2 %>% 
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
  labs(x=gene.OI$genotype, title="G", 
       y=paste("Normalized log2 expression", gene.OI$gene, sep="\n")) +
  geom_smooth(method = "lm", se=FALSE, aes(group=condition),
              color="grey60") +
  
  theme(text = element_text(size=8)) +
  facet_wrap(~condition, ncol=1, scales="free_x")
p5

#### Save plots #####
ggsave("figs/Fig5.eQTL.png", p5, height=6, width = 1.5)
ggsave("figs/Fig5.eQTL.pdf", p5, height=6, width = 1.5)

