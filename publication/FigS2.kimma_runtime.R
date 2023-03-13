library(tidyverse)

load("results/kimma_time_trial.RData")

p1 <- time_trial_df %>%
  mutate(group = ifelse(model %in% c("none","base","metrics","weights",
                                     "covariate","interaction"),
                        "a","b")) %>% 
  mutate(model = fct_recode(factor(model),
                            "start"="none",
                            "none"="base",
                            "interaction +\ncontrast"="contrast",
                            "paired + kinship"="kin",
                            "all preceding"="all",
                            "3 models"="multi"),
         model = fct_relevel(model, c("start","none","metrics","weights",
                                      "covariate","interaction",
                                      "interaction +\ncontrast",
                                      "paired","paired + kinship",
                                      "all preceding",
                                      "3 models"))) %>% 
  filter(!(model %in% c("start"))) %>%
  ggplot(aes(x=model, y=diff)) +
  geom_col() +
  labs(y="Time (min)",
       x="Added to kimma model") +
  theme_classic(base_size = 10) +
  # facet_wrap(~group, scales="free") +
  geom_hline(yintercept = time_trial_df[time_trial_df$model=="base","diff"],
             lty="dashed", color="red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.text.x = element_blank())

p1

ggsave("figs/FigS2.kimma_runtime.png", p1, width=3.5, height = 2.5)
ggsave("figs/FigS2.kimma_runtime.pdf", p1, width=3.5, height = 2.5)
