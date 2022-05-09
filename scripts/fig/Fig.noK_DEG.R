library(tidyverse)
library(ggupset)
library(ggridges)
library(patchwork)

#colors
clr_vec <- c("kimma"="#AA4499","limma"="#117733",
             "dream"="#88CCEE","DESeq2"="#DDCC77")

#### Condition data ####
load("results/condition_fit.RData")

#### ggplot upset ####

p1 <- condition_result %>% 
  #top performing models
  filter(kinship == "no kinship" & paired == "paired" &
           software %in% c("kimma","dream") & 
           weights_type %in% c("no weights","dream weights")) %>%
  filter(FDR < 0.05) %>% 
  mutate(group=paste(weights_type, software)) %>% 
  group_by(gene) %>% 
  summarise(models=list(group)) %>% 
  
  ggplot(aes(x=models)) +
  geom_bar() +
  geom_text(stat='count', aes(label=after_stat(count)), 
            vjust=-1, size=3) +
  scale_x_upset(sets = c("voom weights kimma",
                         "dream weights kimma",
                         "dream weights dream",
                         "no weights kimma",
                         "no weights dream")) +
  theme_classic() +
  labs(x="", y="Genes FDR < 0.05", 
       title="Mtb-infected vs media") +
  lims(y=c(0,10800))

p1

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
  mutate(group = paste(weights_type, software),
         group = factor(group, c("no weights dream","no weights kimma",
                                 "dream weights dream","dream weights kimma"))) %>% 
  #non-overlap genes
  filter(!(gene %in% all_method) & gene %in% method1
         & FDR >= 0.05)

p2 <- ridge_data %>% 
  ggplot(aes(x = FDR, y = group, fill = software)) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.000001) +
  theme_classic() +
  theme(legend.position = "none")  +
  geom_vline(xintercept = c(0.05,0.1,0.25), lty="dashed") +
  labs(y="", x="FDR") +
  scale_x_continuous(breaks=c(0.05,0.1,0.25,0.5,0.75,1),
                    limits = c(0.05,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(drop=FALSE) +
  scale_fill_manual(values=clr_vec)
p2

#### RSTR data ####
load("results/rstr_tb_fit.RData")

#### ggplot upset ####

p3 <- rstr_result %>% 
  #top performing models
  filter(kinship == "no kinship" & paired == "unpaired" &
           software %in% c("kimma","limma")) %>%
  filter(FDR < 0.05) %>% 
  mutate(group=paste(weights, software)) %>% 
  group_by(gene) %>% 
  summarise(models=list(group)) %>% 
  
  ggplot(aes(x=models)) +
  geom_bar() +
  geom_text(stat='count', aes(label=after_stat(count)), 
            vjust=-1,size=3) +
  scale_x_upset(sets = c("voom weights kimma",
                         "voom weights limma",
                         "no weights kimma",
                         "no weights limma")) +
  theme_classic() +
  labs(x="", y="Genes FDR < 0.05", title="RSTR vs LTBI") +
  lims(y=c(0,7))

p3

#### FDR distrib ####
#Select genes signif in all
# all_methodB <- rstr_result %>% 
#   filter(kinship == "no kinship" & paired == "unpaired" & 
#            software %in% c("kimma","limma") & FDR < 0.2) %>% 
#   dplyr::count(gene) %>% 
#   filter(n==4)

#Select genes signif in at least 1
# method1B <- rstr_result %>% 
#   #top performing models
#   filter(kinship == "no kinship" & paired == "unpaired" &
#            software %in% c("kimma","limma")) %>%
#   filter(FDR < 0.2) %>% 
#   distinct(gene) %>% pull(gene)
# 
# ridge_dataB <- rstr_result %>% 
#   #top performing models
#   filter(kinship == "no kinship" & paired == "unpaired" &
#            software %in% c("kimma","limma")) %>% 
#   #Label
#   mutate(group = paste(weights, software),
#          group = factor(group, c("no weights limma","no weights kimma",
#                                  "voom weights limma","voom weights kimma"))) %>% 
#   #non-overlap genes
#   filter(!(gene %in% all_methodB) & gene %in% method1B
#          & FDR >= 0.2) 
# 
# p4 <- ridge_dataB %>% 
#   ggplot(aes(x = FDR, y = group, fill = software)) +
#   geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.000001) +
#   theme_classic() +
#   theme(legend.position = "none")  +
#   geom_vline(xintercept = c(0.2,0.25), lty="dashed") +
#   labs(y="", x="FDR") +
#   scale_x_continuous(breaks=c(0.2,0.25,0.3),
#                      limits = c(0.2,0.3)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_y_discrete(drop=FALSE) +
#   scale_fill_manual(values=clr_vec)
# p4

#### combine plots #####
# p <- plot_spacer()+p1+p2+plot_spacer()+p3+p4 + 
#   plot_annotation(tag_levels = "A") +
#   plot_layout(widths = c(0.1,1,0.4),
#               design = "
#               123
#               456
#               ")
p <-p1+p2+p3+plot_spacer() +plot_spacer()+plot_spacer()+ 
  plot_annotation(tag_levels = "A") +
  plot_layout(widths = c(1,0.2,0.3),
              design="
              143
              256
              ")
# p

ggsave("figs/real_DEG_noKin.png", p, height=5.5, width = 8)
# ggsave("figs/real_DEG_noKin.pdf", p, height=6, width = 9)

##### Addtl results for text #####
#how many hits between 0.05 and 0.1
length(ridge_data$FDR[ridge_data$FDR < 0.1])
length(ridge_data$FDR[ridge_data$FDR < 0.25])/length(ridge_data$FDR)*100
#how many hits between 0.2 and 0.25
length(ridge_dataB$FDR[ridge_dataB$FDR < 0.25])
length(ridge_dataB$FDR[ridge_dataB$FDR < 0.25])/length(ridge_dataB$FDR)*100

