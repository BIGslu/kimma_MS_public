library(tidyverse)
library(ggh4x)
library(patchwork)
library(ggpubr)
library(multcompView)

#colors
clr_vec <- c("1 processor"="#0072B2","6 processors"="#E69F00")

# Data
load("results/simulated_fit.RData")

# summary
time_summ <- time_result %>% 
  group_by(model, software, paired, kinship, weights, weights_type, processors) %>% 
  summarise(meanT=mean(time_s), .groups="drop")

# Plot
p1a <- time_result %>% 
  filter(kinship == "no kinship" & paired=="unpaired" &
           weights_type != "dream weights") %>% 
  ggplot(aes(x=software, y=time_s, group=processors)) +
  geom_point(aes(color=processors), alpha=0.3, size=0.5, 
             position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0,
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights_type, scales="free_x", 
               space = "free") +
  labs(y="Time (sec)", color="") +
  theme(panel.border = element_rect(fill=NA, "black"), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")+
  guides(color = guide_legend(override.aes = list(size=4)))+
  scale_color_manual(values=clr_vec)
# p1a

p1b <- time_result %>% 
  filter(kinship == "no kinship"& paired=="paired") %>% 
  ggplot(aes(x=software, y=time_s/60, group=processors)) +
  geom_point(aes(color=processors), alpha=0.3, size=0.5, 
             position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0,
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights_type, scales="free_x", space = "free") +
  labs(y="Time (min)", color="") +
  theme(panel.border = element_rect(fill=NA, "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size=4))) +
  scale_color_manual(values=clr_vec)
# p1b

#### Add Pval ####
aov.result1 <- aov(as.numeric(time_s)~model, 
                   data=filter(time_result, kinship=="no kinship" &
                                 paired=="unpaired"))
tukey.result1 <- TukeyHSD(aov.result1)

group.letter1 <- data.frame(multcompLetters(tukey.result1$model[,4])$Letters,
                            reversed=TRUE) %>% 
  rownames_to_column("model") %>%
  rename(let=2)  %>% 
  left_join(time_summ) %>% 
  mutate(time_s = 30) %>% 
  #set letters from shortest to longest time
  mutate(let_fct = fct_reorder(factor(let), meanT),
         let_num = as.numeric(let_fct)) %>% 
  arrange(let_fct) %>% 
  mutate(let_new = recode(let_num,
                          "1"="a","2"="b","3"="bc","4"="c","5"="d"))

p1a_signif <- p1a +
  geom_text(data = group.letter1, 
            aes(label = let_new, fill = processors),
            vjust=0, size=3,
            position = position_jitterdodge(dodge.width = 0.8)) +
  lims(y=c(0,32))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
# p1a_signif

#####
aov.result2 <- aov(as.numeric(time_s)~model, 
                   data=filter(time_result, kinship=="no kinship" &
                                 paired=="paired"))
tukey.result2 <- TukeyHSD(aov.result2)

group.letter2 <- data.frame(multcompLetters(tukey.result2$model[,4],
                                            reversed=TRUE)$Letters) %>% 
  rownames_to_column("model") %>%
  rename(let=2)  %>% 
  left_join(time_summ) %>% 
  mutate(time_s = 5.2*60) %>% 
  #set letters from shortest to longest time
  mutate(let_fct = fct_reorder(factor(let), meanT),
         let_num = as.numeric(let_fct)) %>% 
  arrange(let_fct) %>% 
  mutate(let_new = recode(let_num,
                          "1"="a","2"="ab","3"="b","4"="c","5"="d","6"="e","7"="f","8"="g"))
  
p1b_signif <- p1b +
  geom_text(data = group.letter2, 
            aes(label = let_new, fill = processors),
            vjust=0, size=3,
            position = position_jitterdodge(dodge.width = 0.8)) +
  lims(y=c(0,5.6))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
 # p1b_signif

#### Save ####

ggsave(filename="figs/simulated_noK_runtime_unpaired.png", p1a_signif, width=3.5, height=5)
ggsave(filename="figs/simulated_noK_runtime_paired.png", p1b_signif, width=3.5, height=5)

