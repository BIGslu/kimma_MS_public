library(tidyverse)
library(DescTools)
library(ggh4x)
library(patchwork)

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

#### Plot curves ####
p1 <- ss_result_summ %>% 
  filter(kinship == "no kinship") %>% 
  ggplot(aes(x=1-spec.mean, y=sen.mean, color=software, fill=software)) +
  geom_ribbon(aes(ymax=sen.mean+sen.sd, ymin=sen.mean-sen.sd), 
              alpha=0.2, color=NA) +
  geom_path() +
  theme_classic() +
  facet_nested(~kinship+paired+weights_type, scales="free_x") +
  labs(y="Sensitivity", x="1 - Specificity") +
  theme(panel.border = element_rect(fill=NA, "black"))

p2 <- ss_result_summ %>% 
  filter(kinship == "kinship") %>% 
  ggplot(aes(x=1-spec.mean, y=sen.mean, color=software, fill=software)) +
  geom_ribbon(aes(x=1-spec.mean, ymax=sen.mean+sen.sd, ymin=sen.mean-sen.sd), 
              alpha=0.2, color=NA) +
  geom_path() +
  theme_classic() +
  facet_nested(~kinship+paired+weights_type, scales="free_x") +
  labs(y="Sensitivity", x="1 - Specificity") +
  theme(panel.border = element_rect(fill=NA, "black"))

#### Plot AUC ####
p3 <- auc %>% 
  filter(kinship == "no kinship") %>% 
  ggplot(aes(x=software, y=auc, group=software)) +
  geom_point(aes(color=software), alpha=0.2, size=0.5, 
             position = position_jitterdodge(jitter.width = 2, 
                                             dodge.width = 1)) +
  # geom_line() +
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
  guides(color = guide_legend(override.aes = list(size=4)))
# p3

p4 <- auc %>% 
  filter(kinship == "kinship") %>% 
  ggplot(aes(x=software, y=auc, group=software)) +
  geom_point(aes(color=software), alpha=0.2, size=0.5, 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             dodge.width = 0.8)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0,
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights_type, scales="free_x", space = "free") +
  labs(y="AUC", color="") +
  theme(panel.border = element_rect(fill=NA, "black"), 
        legend.position = "none") +
  guides(color = guide_legend(override.aes = list(size=4)))

# p4

#### Plot change in AUC ####
p5 <- auc %>% 
  filter(kinship == "no kinship") %>% 
  select(-model) %>% 
  pivot_wider(names_from = software, values_from = auc) %>% 
  mutate(`kimma - limma`=kimma-limma,
         `kimma - dream`=kimma-dream,
         `kimma - DESeq2`=kimma-DESeq2,
         `limma - dream`=limma-dream,
         `limma - DESeq2`=limma-DESeq2,
         `dream - DESeq2`=dream-DESeq2) %>% 
  pivot_longer(`kimma - limma`:`dream - DESeq2`) %>% 
  drop_na(value) %>% 
  
  ggplot(aes(x=name, y=value, group=name)) +
  geom_jitter(alpha=0.2, size=0.5, width=0.2, height=0) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0, color="red",
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", color="red",
               position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights, scales = "free_x") +
  labs(y="Difference in AUC", color="", x="") +
  theme(panel.border = element_rect(fill=NA, "black")) +
  geom_hline(yintercept = 0, color="red") +
  coord_flip()
# p5


#### Stats AUC ####
# Make groups without software in model name
aov.dat <- auc %>% mutate(x = gsub("^[A-Z]", "", model)) 

# Paired and unpaired t-tests
result <- data.frame()
for(m in c("nnn","nny","ynn","yny")){
  print(m)
  temp <- aov.dat %>% filter(x==m)
  
  fit1 <- broom::tidy(pairwise.t.test(x=temp$auc, g=temp$software, 
                                      p.adjust.method="none", 
                                      paired=TRUE)) %>% 
    mutate(paired=TRUE)
  fit2 <- broom::tidy(pairwise.t.test(x=temp$auc, g=temp$software, 
                                      p.adjust.method="none", 
                                      paired=FALSE)) %>% 
    mutate(paired=FALSE)
  
  #mean differences
  if(m == "yny") {
    diff2 <- temp %>% 
      select(-model) %>% 
      group_by(software) %>% 
      summarise(value = mean(auc)) %>% 
      pivot_wider(names_from = software) %>% 
      mutate(mean.diffKL=kimma-limma,
             mean.diffKD=kimma-dream,
             mean.diffLD=limma-dream)
    diff1 <- temp %>% 
      select(-model) %>% 
      pivot_wider(names_from = software, values_from = auc) %>% 
      mutate(diffKL=kimma-limma,
             diffKD=kimma-dream,
             diffLD=limma-dream) %>% 
      pivot_longer(diffKL:diffLD) %>% 
      group_by(x, name) %>% 
      summarise(mean.diff = mean(value),
                sd.diff = sd(value)) %>% 
      arrange(name)
    
    result <- bind_rows(fit1, fit2) %>% 
      mutate(model=m, 
             mean.diff=c(diff1[diff1$name == "diffKL",]$mean.diff,
                         diff1[diff1$name == "diffKD",]$mean.diff,
                         diff1[diff1$name == "diffLD",]$mean.diff,
                         diff2$mean.diffKL,diff2$mean.diffKD,diff2$mean.diffLD),
             sd.diff = c(diff1[diff1$name == "diffKL",]$sd.diff,
                         diff1[diff1$name == "diffKD",]$sd.diff,
                         diff1[diff1$name == "diffLD",]$sd.diff,
                         NA,NA,NA)) %>% 
      bind_rows(result)
  } else if(m == "nnn") {
    diff2 <- temp %>% 
      select(-model) %>% 
      group_by(software) %>% 
      summarise(value = mean(auc)) %>% 
      pivot_wider(names_from = software) %>% 
      mutate(mean.diffKL=kimma-limma,
             mean.diffKS=kimma-DESeq2,
             mean.diffLS=limma-DESeq2)
    diff1 <- temp %>% 
      select(-model) %>% 
      pivot_wider(names_from = software, values_from = auc) %>% 
      mutate(diffKL=kimma-limma,
             diffKS=kimma-DESeq2,
             diffLS=limma-DESeq2) %>% 
      pivot_longer(diffKL:diffLS) %>% 
      group_by(x, name) %>% 
      summarise(mean.diff = mean(value),
                sd.diff = sd(value)) %>% 
      arrange(name)
    
    result <- bind_rows(fit1, fit2) %>% 
      mutate(model=m, 
             mean.diff=c(diff1[diff1$name == "diffKL",]$mean.diff,
                         diff1[diff1$name == "diffKS",]$mean.diff,
                         diff1[diff1$name == "diffLS",]$mean.diff,
                         diff2$mean.diffKL,diff2$mean.diffKS,diff2$mean.diffLS),
             sd.diff = c(diff1[diff1$name == "diffKL",]$sd.diff,
                         diff1[diff1$name == "diffKS",]$sd.diff,
                         diff1[diff1$name == "diffLS",]$sd.diff,
                         NA,NA,NA)) %>% 
      bind_rows(result)
  } else if (m == "ynn") {
    diff2 <- temp %>% 
      select(-model) %>% 
      group_by(software) %>% 
      summarise(value = mean(auc)) %>% 
      pivot_wider(names_from = software) %>% 
      mutate(mean.diffKL=kimma-limma,
             mean.diffKD=kimma-dream,
             mean.diffKS=kimma-DESeq2,
             mean.diffLD=limma-dream,
             mean.diffLS=limma-DESeq2,
             mean.diffDS=dream-DESeq2)
    diff1 <- temp %>% 
      select(-model) %>% 
      pivot_wider(names_from = software, values_from = auc) %>% 
      mutate(diffKL=kimma-limma,
             diffKD=kimma-dream,
             diffKS=kimma-DESeq2,
             diffLD=limma-dream,
             diffLS=limma-DESeq2,
             diffDS=dream-DESeq2) %>% 
      pivot_longer(diffKL:diffDS) %>% 
      group_by(x, name) %>% 
      summarise(mean.diff = mean(value),
                sd.diff = sd(value)) %>% 
      arrange(name)
    
    result <- bind_rows(fit1, fit2) %>% 
      mutate(model=m, 
             mean.diff=c(diff1[diff1$name=="diffKL",]$mean.diff,
                         diff1[diff1$name=="diffKD",]$mean.diff,
                         diff1[diff1$name=="diffLD",]$mean.diff,
                         diff1[diff1$name=="diffKS",]$mean.diff,
                         diff1[diff1$name=="diffLS",]$mean.diff,
                         diff1[diff1$name=="diffDS",]$mean.diff,
                         diff2$mean.diffKL,diff2$mean.diffKD,diff2$mean.diffLD,
                         diff2$mean.diffKS,diff2$mean.diffLS,diff2$mean.diffDS),
             sd.diff = c(diff1[diff1$name=="diffKL",]$sd.diff,
                         diff1[diff1$name=="diffKD",]$sd.diff,
                         diff1[diff1$name=="diffLD",]$sd.diff,
                         diff1[diff1$name=="diffKS",]$sd.diff,
                         diff1[diff1$name=="diffLS",]$sd.diff,
                         diff1[diff1$name=="diffDS",]$sd.diff,
                         NA,NA,NA,NA,NA,NA)) %>% 
      bind_rows(result)
  } else {
    diff2 <- temp %>% 
      group_by(software) %>% 
      summarise(value = mean(auc)) %>% 
      pivot_wider(names_from = software) %>% 
      mutate(mean.diff=kimma-limma)
    diff1 <- temp %>% 
      select(-model) %>% 
      pivot_wider(names_from = software, values_from = auc) %>% 
      mutate(diff = kimma-limma) %>% 
      group_by(x) %>% 
      summarise(mean.diff = mean(diff),
                sd.diff = sd(diff))
    result <- bind_rows(fit1, fit2) %>% 
      mutate(model=m, 
             mean.diff=c(diff1$mean.diff,diff2$mean.diff),
             sd.diff = c(diff1$sd.diff,NA)) %>% 
      bind_rows(result)
  }
}

#### Save ####
ggsave(filename="figs/AUC/AUC_curve_full.png",
       p1/p2, width=10, height=5)
ggsave(filename="figs/AUC/AUC_curve_full.pdf",
       p1/p2, width=10, height=5)

p1b <- p1 + lims(y=c(0.7,1.02))
p2b <- p2 + lims(y=c(0.7,1.02))
ggsave(filename="figs/AUC/AUC_curve.png",
       p1b/p2b, width=10, height=5)
ggsave(filename="figs/AUC/AUC_curve.pdf",
       p1b/p2b, width=10, height=5)

ggsave(filename="figs/AUC/AUC.png",
       p3/p4, width=5, height=7)
ggsave(filename="figs/AUC/AUC.pdf",
       p3/p4, width=5, height=7)

ggsave(filename="figs/AUC/AUC_diff.png",
       p5, width=8, height=2.5)
ggsave(filename="figs/AUC/AUC_diff.pdf",
       p5, width=8, height=2.5)

#### Add Pval ####
library(ggpubr)

stat.dat <- result %>% 
  dplyr::rename(paired.stat=paired) %>% 
  separate(model, into=c("trash","paired","kinship","weights"),
           sep="", remove = FALSE) %>% 
  mutate(paired = recode_factor(factor(paired),
                                "n"="unpaired","y"="paired"),
         kinship = recode_factor(factor(kinship),
                                 "n"="no kinship","y"="kinship"),
         weights = recode_factor(factor(weights),
                                 "n"="no weights","y"="weights")) %>% 
  filter(paired.stat==FALSE & p.value < 0.05) %>% 
  mutate(p.label = ifelse(p.value<0.001,"***","*"),
         p.value = formatC(p.value, format = "e", digits = 2)) %>% 
  arrange(model, group1, group2) %>% 
  mutate(y.position = c(1,1.05,
                        1.05,1.2,1.15,1.1,1.05,
                        1.05,1.1)) %>% 
  mutate(y.max = y.position+0.1,
         software="kimma")


p3P <- p3 +
  stat_pvalue_manual(stat.dat, label="p.label") +
  geom_point(data=stat.dat, x=NA, aes(y=y.max)) +
  theme(legend.position = "none") 

p3P

stat.dat2 <- data.frame(
  kinship = "no kinship",
  paired = c(rep("unpaired",5), rep("paired",7)),
  weights = c(rep("no weights",3),rep("weights",2),
              rep("no weights",4),rep("weights",3)),
  software = c("kimma","limma","DESeq2",
               "kimma","limma",
               "kimma","limma","dream","DESeq2",
               "kimma","limma","dream"),
  let <- c("a","a","b",
           "a","a",
           "a","b","a","c",
           "a","b","a"),
  auc = 1.01
) %>% 
  mutate(paired = factor(paired, level=c("unpaired","paired")))
p3Pb <- p3 +
  geom_text(data = stat.dat2, mapping = aes(label = let),
            vjust=0, size=3) +
  lims(y=c(0.4,1.03))+
  theme(legend.position = "none") 
p3Pb

ggsave("figs/AUC/AUC_pval.png", p3Pb, width=7, height=3)
ggsave("figs/AUC/AUC_pval.pdf", p3Pb, width=7, height=3)

