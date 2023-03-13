library(tidyverse)
library(DescTools)
library(ggh4x)
library(ggpubr)
library(multcompView)

#colors
clr_vec <- c("kimma"="#AA4499","limma"="#117733",
             "dream"="#88CCEE","DESeq2"="#DDCC77", "edgeR"="#332288")

#### Data ####
load("results/simulated_fit.RData")

ss_result_summ <- ss_result %>% 
  group_by(model, software, paired, kinship, weights, 
           weights_type, processors, cutoff) %>% 
  summarise(sen.mean = mean(sensitivity_p, na.rm=TRUE),
            sen.sd = sd(sensitivity_p, na.rm=TRUE),
            spec.mean = mean(specificity_p, na.rm=TRUE),
            spec.sd = sd(specificity_p, na.rm=TRUE), .groups="drop")

#### Calculate AUC ####
auc <- ss_result %>% 
  group_by(model, software, paired, kinship, weights, 
           weights_type, processors, simulation) %>% 
  summarise(auc = DescTools::AUC(x=1-specificity_p, y=sensitivity_p), 
            .groups = "drop") %>% 
  mutate(model=gsub("fdr","", model)) %>% 
  filter(processors == "1 processor")

#### Plot AUC ####
p3 <- auc %>% 
  filter(kinship == "no kinship") %>% 
  ggplot(aes(x=software, y=auc, group=software)) +
  geom_point(aes(color=software), alpha=0.3, size=0.5, 
             position = position_jitterdodge(jitter.width = 2, 
                                             dodge.width = 1)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0,
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights_type, scales="free_x", space = "free") +
  labs(y="AUC", color="") +
  theme(panel.border = element_rect(fill=NA, "black"), 
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size=4))) +
  scale_color_manual(values=clr_vec)
p3

#### Stats AUC ####
auc_summ <- auc %>% 
  group_by(model, software, paired, kinship, weights, weights_type) %>% 
  summarise(meanAUC=mean(auc), .groups="drop")

aov.result <- aov(auc ~ model, data = filter(auc, kinship=="no kinship"))
group.letter <- data.frame(multcompLetters(tukey.result$model[,4])$Letters) %>% 
  rownames_to_column("model") %>%
  dplyr::rename(let=2)  %>% 
  left_join(auc_summ) %>% 
  mutate(auc = 1.01) %>%
  #set letters from shortest to longest time
  mutate(let_fct = fct_reorder(factor(let), -meanAUC),
         let_num = as.numeric(let_fct)) %>% 
  arrange(let_fct) %>% 
  mutate(let_new = recode(let_num,
                          "1"="a","2"="b","3"="c","4"="d","5"="e",
                          "6"="f","7"="g"))

p3Pb <- p3 +
  geom_text(data = group.letter, mapping = aes(label = let_new),
            vjust=0, size=3) +
  lims(y=c(0.4,1.03))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) 
p3Pb

ggsave("figs/Fig2.simulated_noK_AUC.png", p3Pb, width=6.5, height=4)
ggsave("figs/Fig2.simulated_noK_AUC.pdf", p3Pb, width=6.5, height=4)


