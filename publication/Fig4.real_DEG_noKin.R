library(tidyverse)
library(ggvenn)
# library(patchwork)

#colors
clr_vec <- c("kimma"="#AA4499","limma"="#117733",
             "dream"="#88CCEE","DESeq2"="#DDCC77")

#### Condition data ####
load("results/condition_fit.RData")

#### Venn ####
signif.ls1 <- list()

signif.ls1[["kimma with\nno weights"]] <- condition_result %>% 
  filter(kinship == "no kinship" & paired == "paired" &
           software == "kimma" & 
           weights_type == "no weights") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls1[["dream with\nno weights"]] <- condition_result %>% 
  filter(kinship == "no kinship" & paired == "paired" &
           software == "dream" & 
           weights_type == "no weights") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls1[["kimma with\ndream weights"]] <- condition_result %>% 
  filter(kinship == "no kinship" & paired == "paired" &
           software == "kimma" & 
           weights_type == "dream weights") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls1[["dream with\ndream weights"]] <- condition_result %>% 
  filter(kinship == "no kinship" & paired == "paired" &
           software == "dream" & 
           weights_type == "dream weights") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

p1 <- ggvenn(signif.ls1, show_percentage = FALSE,
       fill_color = c("#AA4499","#88CCEE","#AA4499","#88CCEE"),
       text_size = 3, set_name_size = 3, stroke_size = 0.5) +
  labs(title="Paired\nMtb vs media in LTBI") +
  lims(x=c(-2.8,2.9)) +
  theme(plot.title = element_text(hjust = 0.5))
# p1

#### RSTR data ####
# load("results/rstr_tb_fit.RData")

#### Venn ####
# signif.ls2 <- list()
# 
# signif.ls2[["kimma with\nno weights"]] <- rstr_result %>% 
#   filter(kinship == "no kinship" & paired == "unpaired" &
#            software == "kimma" & 
#            weights == "no weights") %>%
#   filter(FDR < 0.05) %>% 
#   pull(gene)
# 
# signif.ls2[["limma with\nno weights"]] <- rstr_result %>% 
#   filter(kinship == "no kinship" & paired == "unpaired" &
#            software == "limma" & 
#            weights == "no weights") %>%
#   filter(FDR < 0.05) %>% 
#   pull(gene)
# 
# signif.ls2[["kimma with\nvoom weights"]] <- rstr_result %>% 
#   filter(kinship == "no kinship" & paired == "unpaired" &
#            software == "kimma" & 
#            weights == "voom weights") %>%
#   filter(FDR < 0.05) %>% 
#   pull(gene)
# 
# signif.ls2[["limma with\nvoom weights"]] <- rstr_result %>% 
#   filter(kinship == "no kinship" & paired == "unpaired" &
#            software == "limma" & 
#            weights == "voom weights") %>%
#   filter(FDR < 0.05) %>% 
#   pull(gene)
# 
# p3 <- ggvenn(signif.ls2, show_percentage = FALSE,
#              fill_color = c("#AA4499","#117733","#AA4499","#117733"),
#              text_size = 3, set_name_size = 3, stroke_size = 0.5) +
#   labs(title="B  Unpaired") +
#   lims(x=c(-2.7,2.8))
# p3

#### combine plots #####

ggsave("figs/Fig4.real_DEG_noKin.png", p1, height=3, width = 4)
ggsave("figs/Fig4.real_DEG_noKin.pdf", p1, height=3, width = 4)

