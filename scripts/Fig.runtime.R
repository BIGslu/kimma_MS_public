library(tidyverse)
library(ggh4x)
library(patchwork)

# Data
load("results/model_fit/simulated_models_all.RData")

time_summ <- time_result %>% 
  group_by(software, paired, kinship, weights, processors) %>% 
  summarise(.groups="drop",
            s = sd(time_s),
            time_s = mean(time_s))
# Plot
p1a <- time_result %>% 
  filter(kinship == "no kinship" & paired=="unpaired") %>% 
  ggplot(aes(x=software, y=time_s, group=processors)) +
  geom_point(aes(color=processors), alpha=0.2, size=0.5, 
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0,
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights, scales="free_x", space = "free") +
  labs(y="Time (sec)", color="") +
  theme(panel.border = element_rect(fill=NA, "black"), legend.position = "none")

p1b <- time_result %>% 
  filter(kinship == "no kinship"& paired=="paired") %>% 
  ggplot(aes(x=software, y=time_s, group=processors)) +
  geom_point(aes(color=processors), alpha=0.2, size=0.5, 
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", width=0,
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights, scales="free_x", space = "free") +
  labs(y="Time (sec)", color="") +
  theme(panel.border = element_rect(fill=NA, "black"), legend.position = "none")
# p1

p2 <- time_result %>% 
  filter(kinship == "kinship") %>% 
  ggplot(aes(x=software, y=time_s, group=processors)) +
  geom_point(aes(color=processors), alpha=0.2, size=0.5,
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0, 
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights, scales="free_x", space = "free") +
  labs(y="Time (sec)", color="") +
  theme(panel.border = element_rect(fill=NA, "black"), legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size=4)))

# p2

p1 <- (p1a+p1b) + plot_layout(widths=c(8,12))
p <- p1/p2 
# p

ggsave(filename="figs/time/simulated_time.pdf", p, width=7, height=6)
ggsave(filename="figs/time/simulated_time.png", p, width=7, height=6)

