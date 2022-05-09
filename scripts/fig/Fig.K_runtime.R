library(tidyverse)
library(ggh4x)
library(patchwork)

#colors
clr_vec <- c("1 processor"="#0072B2","6 processors"="#E69F00")

# Data
load("results/simulated_fit.RData")

# Plot
p2a <- time_result %>% 
  filter(kinship == "kinship" & paired == "unpaired") %>% 
  ggplot(aes(x=software, y=time_s/60, group=processors)) +
  geom_point(aes(color=processors), alpha=0.3, size=0.5,
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0, 
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights_type, scales="free_x", space = "free") +
  labs(y="Time (min)", color="", x="kimma") +
  theme(panel.border = element_rect(fill=NA, "black"), legend.position = "none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  scale_color_manual(values=clr_vec)
# p2a

p2b <- time_result %>% 
  filter(kinship == "kinship" & paired == "paired") %>% 
  ggplot(aes(x=software, y=time_s/60, group=processors)) +
  geom_point(aes(color=processors), alpha=0.3, size=0.5,
             position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0, 
               position=position_dodge(width = 0.8)) +
  stat_summary(fun=mean, geom="point", position=position_dodge(width = 0.8),
               size=1, shape="square") +
  theme_classic() +
  facet_nested(~kinship+paired+weights_type, scales="free_x", space = "free") +
  labs(y="Time (min)", color="", x="kimma") +
  theme(panel.border = element_rect(fill=NA, "black"), legend.position = "bottom",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  lims(y=c(0,3.6)) +
  scale_color_manual(values=clr_vec)

p2 <- (p2a+p2b)
p2

# ggsave(filename="figs/simulated_K_runtime.pdf", p2, width=7, height=3.5)
ggsave(filename="figs/simulated_K_runtime.png", p2, width=7, height=3.5)

