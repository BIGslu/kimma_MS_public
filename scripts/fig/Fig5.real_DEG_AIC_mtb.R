library(tidyverse)
library(ggvenn)
library(patchwork)

#colors
clr_vec <- c("no kinship"="#FFB000","none"="#DC267F",
             "kinship"="#648FFF")

#### Data ####
load("results/condition_subset_fit.RData")

#### AIC ####
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
  scale_color_manual(values = clr_vec) +
  labs(y = "Unrelated\n\nAIC kinship", x="AIC no kinship")

p1_inset <- p1 +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  geom_abline(intercept = 2, slope = 1, lty="dashed") +
  geom_abline(intercept = -2, slope = 1, lty="dashed") +
  scale_x_continuous(breaks=c(0,2,4), limits = c(0,4)) +
  scale_y_continuous(breaks=c(0,2,4), limits = c(0,4))

p1_all <- p1 + 
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
  labs(color="Best fit", y="Related\n\nAIC kinship",
       x="AIC no kinship") +
  scale_color_manual(values = clr_vec) +
  theme(legend.position = "bottom")

p2_inset <- p2 +
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  geom_abline(intercept = 2, slope = 1, lty="dashed") +
  geom_abline(intercept = -2, slope = 1, lty="dashed") +
  scale_x_continuous(breaks=c(0,2,4), limits = c(0,4)) +
  scale_y_continuous(breaks=c(0,2,4), limits = c(0,4))

p2_all <- p2 + 
  inset_element(p2_inset, 0.5,0.01,1,0.49) #lbrt
# p2_all

#### Venn ####

signif.ls1 <- list()

signif.ls1[["no kinship"]] <- condition_subset_result %>% 
  filter(kinship == "no kinship" & paired == "paired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "unrelated") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls1[["kinship"]] <- condition_subset_result %>% 
  filter(kinship == "kinship" & paired == "paired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "unrelated") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p3 <- ggvenn(signif.ls1, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5)
# p3

signif.ls2 <- list()

signif.ls2[["no kinship"]] <- condition_subset_result %>% 
  filter(kinship == "no kinship" & paired == "paired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls2[["kinship"]] <- condition_subset_result %>% 
  filter(kinship == "kinship" & paired == "paired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p4 <- ggvenn(signif.ls2, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5) 
# p4

#### Combine ####
p_final <- p1_all + p3 + p2_all + p4 +
  plot_annotation(tag_levels = "A",
                  title = "Mtb vs media in LTBI")
p_final

#### Save ####
ggsave("figs/Fig5.real_DEG_AIC_mtb.png", p_final, height=6, width = 6)
# ggsave("figs/Fig5.real_DEG_AIC_mtb.pdf", p_final, height=6, width = 6)

#### Check fit new DEG with kinship ####
outersect <- function(x, y) {
  unique(c(setdiff(x, y),
         setdiff(y, x)))
}

outside1 <- setdiff(signif.ls1[[1]], signif.ls1[[2]])
outside2 <- setdiff(signif.ls1[[2]], signif.ls1[[1]])

condition_subset_metric %>% 
  filter(subset == "unrelated" & weights =="voom weights" &
           gene %in% c(outside1,outside2)) %>% 
  select(kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) %>% 
  mutate(diff = kinship-`no kinship`,
         best = ifelse(diff > 2, "no kinship",
                       ifelse(diff < -2, "kinship", "none")),
         group = ifelse(gene %in% outside1, "no kinship", "kinship")) %>% 
  count(group, best)

outside1 <- setdiff(signif.ls2[[1]], signif.ls2[[2]])
outside2 <- setdiff(signif.ls2[[2]], signif.ls2[[1]])

condition_subset_metric %>% 
  filter(subset == "unrelated" & weights =="voom weights" &
           gene %in% c(outside1,outside2)) %>% 
  select(kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) %>% 
  mutate(diff = kinship-`no kinship`,
         best = ifelse(diff > 0, "no kinship",
                       ifelse(diff < 0, "kinship", "none")),
         group = ifelse(gene %in% outside1, "no kinship", "kinship")) %>% 
  count(group, best)
