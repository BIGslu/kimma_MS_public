library(tidyverse)
library(ggvenn)
library(patchwork)

#colors
clr_vec <- c("no kinship"="#FFB000","none"="#DC267F",
             "kinship"="#648FFF")

#### Data ####
load("results/rstr_tb_subset_fit.RData")

#### AIC ####
aic <- rstr_subset_metric %>%
  select(subset, weights, kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) %>% 
  mutate(subset = str_to_title(subset),
         subset = factor(subset, levels=c("Unrelated","Related"))) %>% 
  #final models
  filter(weights == "voom weights")  %>% 
  #Change in AIC
  mutate(aic.delta = `no kinship`-kinship) 


#### AIC ####
p1 <- aic %>% 
  
  ggplot(aes(x=subset, y=aic.delta)) +
  geom_violin(fill="grey90") +
  geom_hline(yintercept = c(-2,2), lty="dashed", color="grey60") + 
  stat_summary(fun.y=mean, geom="point", size=2, color="black", shape="square") +
  theme_classic() +
  labs(y = "Per gene change in AIC\nBetter fit without kinship <--   --> Better fit with kinship       ", 
       x="") 
#p1

#### Venn ####
signif.ls1 <- list()

signif.ls1[["no kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "no kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "unrelated") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls1[["kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "unrelated") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p2 <- ggvenn(signif.ls1, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5)+
  annotate(geom = "text", label = "Unrelated", x = 0, y = -1.3, size=3)
# p2

signif.ls2 <- list()

signif.ls2[["no kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "no kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls2[["kinship"]] <- rstr_subset_result %>% 
  filter(kinship == "kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "voom weights" & subset == "related") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p3 <- ggvenn(signif.ls2, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 2.5, set_name_size = 3, stroke_size = 0.5) +
  annotate(geom = "text", label = "Related", x = 0, y = -1.3, size=3)
# p3

#### Combine ####
design <- "
AB
AC
"
p_final <- p1 + p2 + p3 +
  plot_annotation(tag_levels = "A",
                  title = "Unpaired\nRSTR vs LTBI in Mtb") +
  plot_layout(design = design)
p_final

#### Save ####
ggsave("figs/FigS3.real_DEG_AIC_rstr.png", p_final, height=4.5, width = 4.5)
ggsave("figs/FigS3.real_DEG_AIC_rstr.pdf", p_final, height=4.5, width = 4.5)

#### Check fit new DEG ####
outersect <- function(x, y) {
  unique(c(setdiff(x, y),
           setdiff(y, x)))
}

outside1 <- setdiff(signif.ls2[[1]], signif.ls2[[2]])
outside2 <- setdiff(signif.ls2[[2]], signif.ls2[[1]])

rstr_subset_metric %>% 
  filter(subset == "unrelated" & weights =="voom weights" &
           gene %in% c(outside1,outside2)) %>% 
  select(kinship, gene, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) %>% 
  mutate(diff = kinship-`no kinship`,
         best = ifelse(diff > 0, "no kinship",
                       ifelse(diff < 0, "kinship", "none")),
         group = ifelse(gene %in% outside2, "kinship", "no kinship")) %>% 
  count(group, best)

