library(tidyverse)
library(ggvenn)
# library(venn)
# library(ggupset)
library(ggridges)
library(patchwork)

#colors
clr_vec <- c("no kinship"="#FFB000","none"="#DC267F",
             "kinship"="#648FFF")

#### Condition data ####
load("results/condition_subset_fit.RData")
load("results/rstr_subset_fit.RData")

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

p1 <- ggvenn(signif.ls1, show_percentage = FALSE,
       fill_color = c("#FFB000","#648FFF"),
       text_size = 3, set_name_size = 3, stroke_size = 0.5) +
  labs(title="A  Unrelated\nkimma, paired, weights") +
  lims(x=c(-2.8,2.9))
# p1

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

p2 <- ggvenn(signif.ls2, show_percentage = FALSE,
             fill_color = c("#FFB000","#648FFF"),
             text_size = 3, set_name_size = 3, stroke_size = 0.5) +
  labs(title="B  Related\nkimma, paired, weights") +
  lims(x=c(-2.8,2.9))
p2

#### combine plots #####

p <-p1+p2
# p

ggsave("figs/real_DEG_Kin.png", p, height=3, width = 6)
# ggsave("figs/real_DEG_Kin.pdf", p, height=6, width = 9)

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
             text_size = 3, set_name_size = 3, stroke_size = 0.5) +
  labs(title="C  Unrelated\nkimma, unpaired, weights") +
  lims(x=c(-2.8,2.9))
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
             text_size = 3, set_name_size = 3, stroke_size = 0.5) +
  labs(title="D  Related\nkimma, unpaired, weights") +
  lims(x=c(-2.8,2.9))
# p4

#### combine plots #####

p <-p1+p2+p3+p4
 p

ggsave("figs/real_DEG_Kin.png", p, height=6, width = 6)
# ggsave("figs/real_DEG_Kin.pdf", p, height=6, width = 9)
