library(tidyverse)
library(patchwork)

#colors
clr_vec <- c("no kinship"="#FFB000","none"="#DC267F",
             "kinship"="#648FFF")

#### Condition data ####
load("results/condition_subset_fit.RData")

aic <- condition_subset_metric %>%
  select(subset, weights, kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) 

#### Unrelated, voom weights ####

p1 <- aic %>% 
  filter(subset == "unrelated" & weights == "voom weights") %>% 
  mutate(best = ifelse(kinship<`no kinship`, "kinship",
                       ifelse(`no kinship`<kinship, "no kinship", "none")),
         best = factor(best, c("no kinship","none","kinship"))) %>% 
  arrange(best) %>% 
  
  ggplot(aes(x=`no kinship`, y=kinship, color=best)) +
  geom_point(alpha=0.3) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  coord_fixed() +
  theme(legend.position = "none") +
  scale_color_manual(values = clr_vec)

p1_inset <- p1 +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  geom_abline(intercept = 2, slope = 1, lty="dashed") +
  geom_abline(intercept = -2, slope = 1, lty="dashed") +
  scale_x_continuous(breaks=c(0,2,4), limits = c(0,4)) +
  scale_y_continuous(breaks=c(0,2,4), limits = c(0,4))

p1_all <- p1 + ggtitle("Mtb-infected vs media\nUnrelated AIC") +
  inset_element(p1_inset, 0.5,0.01,1,0.49) #lbrt
# p1_all

#### Related, voom weights ####

p2 <- aic %>% 
  filter(subset == "related" & weights == "voom weights") %>% 
  mutate(best = ifelse(kinship<`no kinship`, "kinship",
                       ifelse(`no kinship`<kinship, "no kinship", "none")),
         best = factor(best, c("no kinship","none","kinship"))) %>% 
  arrange(best) %>% 
  
  ggplot(aes(x=`no kinship`, y=kinship, color=best)) +
  geom_point(alpha=0.3) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  coord_fixed() +
  labs(color="Best fit\nwith weights") +
  scale_color_manual(values = clr_vec)

p2_inset <- p2 +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  geom_abline(intercept = 2, slope = 1, lty="dashed") +
  geom_abline(intercept = -2, slope = 1, lty="dashed") +
  scale_x_continuous(breaks=c(0,2,4), limits = c(0,4)) +
  scale_y_continuous(breaks=c(0,2,4), limits = c(0,4))

p2_all <- p2 + ggtitle("Mtb-infected vs media\nRelated AIC") +
  inset_element(p2_inset, 0.5,0.01,1,0.49) #lbrt
# p2_all

#### RSTR data ####
load("results/rstr_subset_fit.RData")

aic2 <- rstr_subset_metric %>%
  select(subset, weights, kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) 

#### Unrelated, voom weights ####

p3 <- aic2 %>% 
  filter(subset == "unrelated" & weights == "weights") %>% 
  mutate(best = ifelse(kinship<`no kinship`, "kinship",
                       ifelse(`no kinship`<kinship, "no kinship", "none")),
         best = factor(best, c("no kinship","none","kinship"))) %>% 
  arrange(best) %>% 
  
  ggplot(aes(x=`no kinship`, y=kinship, color=best)) +
  geom_point(alpha=0.3) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  coord_fixed() +
  theme(legend.position = "none") +
  scale_color_manual(values = clr_vec) +
  labs(title="RSTR vs LTBI\nUnrelated AIC")

#### Related, voom weights ####

p4 <- aic2 %>% 
  filter(subset == "related" & weights == "weights") %>% 
  mutate(best = ifelse(kinship<`no kinship`, "kinship",
                       ifelse(`no kinship`<kinship, "no kinship", "none")),
         best = factor(best, c("no kinship","none","kinship"))) %>% 
  arrange(best) %>% 
  
  ggplot(aes(x=`no kinship`, y=kinship, color=best)) +
  geom_point(alpha=0.3) +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic() +
  coord_fixed() +
  labs(color="Best fit\nwith weights", title = "RSTR vs LTBI\nRelated AIC") +
  scale_color_manual(values = clr_vec)+
  geom_abline(intercept = 0, slope = 1)

#### combine ####
p_combo <- p1_all + p2_all + p3 + p4 +
  plot_annotation(tag_levels = "A")
# p_combo

ggsave("figs/real_AIC_weights.png", p_combo, width=8, height=8)

#### Addtl metrics ####

aic %>% 
  mutate(diff = `no kinship`-kinship) %>% 
  group_by(subset, weights) %>% 
  summarise(meanD = mean(diff),
            sdD = sd(diff),
            minD = min(diff),
            maxD = max(diff))

aic %>% 
  mutate(diff = `no kinship`-kinship,
         sign = ifelse(diff<0, "neg", "pos")) %>% 
  filter(abs(diff)>2) %>% 
  count(weights, subset, sign)
# pos is best fit kinship, neg is best fit no kinship

length(unique(aic$gene))

###
aic2 %>% 
  mutate(diff = `no kinship`-kinship) %>% 
  group_by(subset, weights) %>% 
  summarise(meanD = mean(diff),
            sdD = sd(diff),
            minD = min(diff),
            maxD = max(diff))

aic2 %>% 
  mutate(diff = `no kinship`-kinship,
         sign = ifelse(diff<0, "neg", "pos")) %>% 
  filter(abs(diff)>2) %>% 
  count(weights, subset, sign)
# pos is best fit kinship, neg is best fit no kinship
