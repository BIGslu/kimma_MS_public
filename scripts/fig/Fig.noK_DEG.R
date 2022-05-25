library(tidyverse)
library(ggvenn)
# library(venn)
# library(ggupset)
library(ggridges)
library(patchwork)

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

# p1 <- venn(signif.ls1, zcolor = "style",
#            box=FALSE, ggplot = TRUE) +
#   labs(title="A") +
#   theme(plot.title = element_text(size=14))

p1 <- ggvenn(signif.ls1, show_percentage = FALSE,
       fill_color = c("#AA4499","#88CCEE","#AA4499","#88CCEE"),
       text_size = 3, set_name_size = 3, stroke_size = 0.5) +
  labs(title="A  Paired") +
  lims(x=c(-2.8,2.9))
# p1

#### FDR distrib ####
#Select genes signif in all
all_method <- condition_result %>%
  filter(kinship == "no kinship" & paired == "paired" &
           weights_type != "voom weights",
           software %in% c("kimma","dream") & FDR < 0.05) %>%
  dplyr::count(gene) %>%
  filter(n==4)

#Select genes signif in at least 1
method1 <- condition_result %>% 
  #top performing models
  filter(kinship == "no kinship" & paired == "paired" &
           software %in% c("kimma","dream")) %>%
  filter(FDR < 0.05) %>% 
  distinct(gene) %>% pull(gene)

ridge_data <- condition_result %>% 
  #top performing models
  filter(kinship == "no kinship" & paired == "paired" &
           software %in% c("kimma","dream") &
           weights_type %in% c("no weights", "dream weights")) %>% 
  #Label
  mutate(group = paste(software, "with\n", weights_type),
         group = factor(group, c("dream with\n no weights",
                                 "kimma with\n no weights",
                                 "dream with\n dream weights",
                                 "kimma with\n dream weights"))) %>% 
  #non-overlap genes
  filter(!(gene %in% all_method) & gene %in% method1
         & FDR >= 0.05)

p2 <- ridge_data %>% 
  ggplot(aes(x = FDR, y = group, fill = software)) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.000001) +
  theme_classic() +
  theme(legend.position = "none")  +
  geom_vline(xintercept = c(0.05,0.1,0.25), lty="dashed") +
  labs(y="", x="FDR", title="B") +
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.5,0.75,1),
                    limits = c(0.05,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(drop=FALSE) +
  scale_fill_manual(values=clr_vec)
# p2

#### RSTR data ####
load("results/rstr_tb_fit.RData")

#### Venn ####
signif.ls2 <- list()

signif.ls2[["kimma with\nno weights"]] <- rstr_result %>% 
  filter(kinship == "no kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "no weights") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls2[["limma with\nno weights"]] <- rstr_result %>% 
  filter(kinship == "no kinship" & paired == "unpaired" &
           software == "limma" & 
           weights == "no weights") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls2[["kimma with\nvoom weights"]] <- rstr_result %>% 
  filter(kinship == "no kinship" & paired == "unpaired" &
           software == "kimma" & 
           weights == "voom weights") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

signif.ls2[["limma with\nvoom weights"]] <- rstr_result %>% 
  filter(kinship == "no kinship" & paired == "unpaired" &
           software == "limma" & 
           weights == "voom weights") %>%
  filter(FDR < 0.05) %>% 
  pull(gene)

# p3 <- venn(signif.ls2, zcolor = "style",
#            box=FALSE, ggplot = TRUE) +
#   labs(title="C") +
#   theme(plot.title = element_text(size=14))

p3 <- ggvenn(signif.ls2, show_percentage = FALSE,
             fill_color = c("#AA4499","#117733","#AA4499","#117733"),
             text_size = 3, set_name_size = 3, stroke_size = 0.5) +
  labs(title="C  Unpaired") +
  lims(x=c(-2.7,2.8))
# p3

#### combine plots #####

p <-p1+p2+p3+plot_spacer()
# p

ggsave("figs/real_DEG_noKin.png", p, height=5.5, width = 8)
# ggsave("figs/real_DEG_noKin.pdf", p, height=6, width = 9)

##### Addtl results for text #####
#how many hits between 0.05 and a higher cutoff
length(ridge_data$FDR[ridge_data$FDR < 0.1])
length(ridge_data$FDR[ridge_data$FDR < 0.25])/length(ridge_data$FDR)*100

