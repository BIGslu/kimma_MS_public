library(tidyverse)
library(SNPRelate)
library(biomaRt)

#### SNP data ####
gds <- snpgdsOpen("data/merged_Omni_MegaEx_randomFounder_Hawn_filter_LD.gds")
geno <- snpgdsGetGeno(gdsobj = gds, with.id=TRUE)
colnames(geno$genotype) <- geno$snp.id
rownames(geno$genotype) <- geno$sample.id
geno.df <- as.data.frame(geno$genotype)
geno.df[1:3,1:3]

snpgdsClose(gds)

#### Filter genes of interest ####
#Get genes significant for Mtb in LTBI samples
load("results/condition_subset_fit.RData")

deg1 <- condition_subset_result %>% 
  filter(kinship == "kinship" & paired == "paired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)


#### ENSEMBL ####
#Get Ensembl ID for genes of interest
grch38 <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

gene.anno <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters = "hgnc_symbol", 
                   values = deg1, 
                   mart = grch38)

#Annotate SNPs to Ensembl gene ID
grch38.snp <- useMart(biomart="ENSEMBL_MART_SNP", dataset="hsapiens_snp")

#### SNP anno ####
## Note getting all SNP at once times out from ENSEMBL server
## Instead get 1000 SNP chunks at a time

snp.anno <- data.frame()
for(i in seq(1, ncol(geno.df), by=1000)){
  print(i)
  
  #Recode last chunk to not go beyond data
  i2 <- i+999
  if(i2 > ncol(geno.df)){ i2 <- ncol(geno.df) }
  
  snp.temp <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id"),
                    filters = "snp_filter",
                    values = colnames(geno.df)[i:i2],
                    mart = grch38.snp) 
  
  if(nrow(snp.temp) > 0){ snp.anno <- bind_rows(snp.temp, snp.anno)}
}

save(snp.anno, file="data/snp.anno.RData")

#### Data cleaning ####
# Combine gene HGNC and SNP anno
snp.anno.all <- snp.anno %>%
  filter(ensembl_gene_stable_id %in% gene.anno$ensembl_gene_id) %>% 
  distinct() %>% 
  inner_join(gene.anno, by=c("ensembl_gene_stable_id"="ensembl_gene_id")) %>% 
  rename(ensembl=ensembl_gene_stable_id, gene=hgnc_symbol) %>% 
  dplyr::select(gene, refsnp_id)

#Filter genotypes data
geno.df.filter <- geno.df %>% 
  dplyr::select(all_of(snp.anno.all$refsnp_id)) %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("refsnp_id") %>% 
  left_join(snp.anno.all) %>% 
  dplyr::select(refsnp_id, gene, everything())

#Rename FULLIDNO to ptID
key <- read_csv("data/RSTR_ID_key.csv") %>% 
  filter(FULLIDNO %in% colnames(geno.df.filter)[-c(1:2)]) %>% 
  mutate(ID_sort = factor(FULLIDNO, levels=colnames(geno.df.filter)[-c(1:2)])) %>% 
  arrange(ID_sort)

geno.df.filter.rename <- geno.df.filter %>% 
  dplyr::select(refsnp_id, gene, all_of(key$ID_sort))
colnames(geno.df.filter.rename) <- c("refsnp_id","gene",key$ptID)

#### Save ####
geno.final <- geno.df.filter.rename

save(geno.final, deg1, file="data/RSTR_snp.RData")
