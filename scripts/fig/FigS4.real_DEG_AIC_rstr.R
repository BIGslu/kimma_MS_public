library(tidyverse)
library(ggvenn)
library(patchwork)

#colors
clr_vec <- c("no kinship"="#FFB000","none"="#DC267F",
             "kinship"="#648FFF")

#### Data ####
load("results/rstr_subset_fit.RData")

#### AIC ####
aic <- rstr_subset_metric %>%
  select(subset, weights, kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) 

#### Unrelated, voom weights ####
p1 <- aic %>% 
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
  labs(y = "Unrelated\n\nAIC kinship", x="AIC no kinship")

#p1

#### Related, voom weights ####
p2 <- aic %>% 
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
  labs(color="Best fit", y="Related\n\nAIC kinship",
       x="AIC no kinship") +
  scale_color_manual(values = clr_vec) +
  theme(legend.position = "bottom")

#p2

#### Venn ####
signif.ls1 <- list()

signif.ls1[["no kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "no kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "weights" & subset == "unrelated") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls1[["kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "weights" & subset == "unrelated") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p3 <- ggvenn(signif.ls1, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5)
# p3

signif.ls2 <- list()

signif.ls2[["no kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "no kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls2[["kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p4 <- ggvenn(signif.ls2, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5) 
# p4

#### Combine ####
p_final <- p1 + p3 + p2 + p4 +
  plot_annotation(tag_levels = "A",
                  title = "RSTR vs LTBI in media")
p_final

#### Save ####
ggsave("figs/FigS4.real_DEG_AIC_rstr.png", p_final, height=6, width = 6)
# ggsave("figs/FigS4.real_DEG_AIC_rstr.pdf", p_final, height=6, width = 6)

#### Check fit new DEG ####
outersect <- function(x, y) {
  unique(c(setdiff(x, y),
           setdiff(y, x)))
}

outside1 <- setdiff(signif.ls2[[1]], signif.ls2[[2]])
outside2 <- setdiff(signif.ls2[[2]], signif.ls2[[1]])

rstr_subset_metric %>% 
  filter(subset == "unrelated" & weights =="weights" &
           gene %in% c(outside1,outside2)) %>% 
  select(kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) %>% 
  mutate(diff = kinship-`no kinship`,
         best = ifelse(diff > 0, "no kinship",
                       ifelse(diff < 0, "kinship", "none")),
         group = ifelse(gene %in% outside2, "kinship", "no kinship")) %>% 
  count(group, best)

