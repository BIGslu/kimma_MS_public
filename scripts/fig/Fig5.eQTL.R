library(tidyverse)

attach("results/rstr_subset_eqtl.RData")
attach("data/RSTR_data_all.RData")
load("data/RSTR_snp.RData")

dat.eqtl <- eqtl$lmerel %>% 
  filter(variable != '(1 | ptID)') %>% 
  mutate(FDR = p.adjust(pval, method = "BH"))

dat.eqtl %>% 
  pivot_longer(c(pval,FDR)) %>%
  filter(value < 0.05) %>% 
  count(name)

gene.OI <- dat.eqtl %>% 
  filter(FDR < 0.05)

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

temp2 %>% 
  ggplot(aes(x=as.character(genotype), y=E))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.2, height=0) +
  theme_classic() +
  labs(x=gene.OI$genotype, 
       y=paste("Normalized log2 expression", gene.OI$gene, sep="\n")) +
  geom_smooth(method = "lm", se=FALSE, aes(group=1), color="red")
