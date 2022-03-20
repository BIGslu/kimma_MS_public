library(tidyverse)
library(DescTools)
library(ggh4x)
library(patchwork)

#### Data ####
load("results/model_fit/simulated_models_all.RData")

ss_result_summ <- ss_result %>% 
  group_by(model, software, paired, kinship, weights, processors, cutoff) %>% 
  summarise(sen.mean = mean(sensitivity_p, na.rm=TRUE),
            sen.sd = sd(sensitivity_p, na.rm=TRUE),
            spec.mean = mean(specificity_p, na.rm=TRUE),
            spec.sd = sd(specificity_p, na.rm=TRUE), .groups="drop")

#### Calculate AUC ####
auc <- ss_result %>% 
  group_by(model, software, paired, kinship, weights, processors, simulation) %>% 
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
  facet_nested(~kinship+paired+weights, scales="free_x") +
  labs(y="Sensitivity", x="1 - Specificity") +
  theme(panel.border = element_rect(fill=NA, "black"))

p2 <- ss_result_summ %>% 
  filter(kinship == "kinship") %>% 
  ggplot(aes(x=1-spec.mean, y=sen.mean, color=software, fill=software)) +
  geom_ribbon(aes(x=1-spec.mean, ymax=sen.mean+sen.sd, ymin=sen.mean-sen.sd), 
              alpha=0.2, color=NA) +
  geom_path() +
  theme_classic() +
  facet_nested(~kinship+paired+weights, scales="free_x") +
  labs(y="Sensitivity", x="1 - Specificity") +
  theme(panel.border = element_rect(fill=NA, "black"))

#### Plot AUC ####
p3 <- auc %>% 
  filter(kinship == "no kinship") %>% 
  ggplot(aes(x=software, y=auc, group=software)) +
  geom_point(aes(color=software), alpha=0.2, size=0.5, 
             position = position_jitterdodge(jitter.width = 0.3, 
                                             dodge.width = 0.8)) +
  # geom_line() +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0,
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights, scales="free_x", space = "free") +
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
  facet_nested(~kinship+paired+weights, scales="free_x", space = "free") +
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
         `limma - dream`=limma-dream) %>% 
  pivot_longer(`kimma - limma`:`limma - dream`) %>% 
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
  if(m %in% c("yny","ynn")) {
    diff2 <- temp %>% 
      group_by(software) %>% 
      summarise(value = mean(auc)) %>% 
      pivot_wider(names_from = software) %>% 
      mutate(mean.diffKL=kimma-limma,
             mean.diffKD=kimma-dream,
             mean.diffLD=limma-dream)
    diff1 <- temp %>% 
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
             mean.diff=c(diff1$mean.diff[2],diff1$mean.diff[1],
                         diff1$mean.diff[3],diff2$mean.diffKL,
                         diff2$mean.diffKD,diff2$mean.diffLD),
             sd.diff = c(diff1$sd.diff[2],diff1$sd.diff[1],
                         diff1$sd.diff[3],NA,NA,NA)) %>% 
      bind_rows(result)
    
  } else {
  diff2 <- temp %>% 
    group_by(software) %>% 
    summarise(value = mean(auc)) %>% 
    pivot_wider(names_from = software) %>% 
    mutate(mean.diff=kimma-limma)
  diff1 <- temp %>% 
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

