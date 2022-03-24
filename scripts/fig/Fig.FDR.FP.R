library(tidyverse)

load("results/model_fit/simulated_models_all.RData")

#Summarise false positives at P < 0.01 and 0.001

ss_result %>% 
  filter(p == 0.01 | p == 0.001) %>%
  mutate(p.group = paste("P <", p)) %>% 
  group_by(model, software, paired, kinship, weights, processors, p.group) %>% 
  summarise(meanFP = mean(FP),
            sdFP = sd(FP), .groups = "drop") %>% 
  
ggplot(aes(x=software, y=meanFP, group=p.group)) +
  geom_bar(aes(fill=p.group), stat="identity",
             position = position_dodge()) +
  geom_errorbar(aes(ymin=meanFP-sdFP, ymax=meanFP+sdFP), width=0.2,
                position=position_dodge(0.9)) +
  theme_classic() +
  facet_nested(~kinship+paired+weights, scales="free_x", space = "free") +
  labs(y="False positives", color="") +
  theme(panel.border = element_rect(fill=NA, "black"),
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size=4)))
