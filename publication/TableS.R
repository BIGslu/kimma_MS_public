library(tidyverse)

#### AUC ####
load("results/simulated_fit.RData")

ss_compare <- ss_result %>% 
  filter(processors == "1 processor" &  kinship == "no kinship" & cutoff == 0.05) %>% 
  select(kinship:processors, cutoff,
         TP_fdr:TN_fdr, sensitivity_fdr, specificity_fdr) %>% 
  group_by(kinship, paired, weights, weights_type, 
           software, processors, cutoff) %>% 
  summarise_all(~mean(.)) %>% 
  rename_all(~gsub("_fdr", "", .))

write_csv(ss_compare, "publication/TableS1.AUC.csv")
