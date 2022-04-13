vec1 <- c("a","b","c")
vec2 <- c("b","d","e")

library(venn)

vend.dat <- venn(list(vec1,vec2))

######
library(tidyverse)
load("results/rstr_subset_fit.RData")

rstr_subset_metric %>% 
  filter(subset == "unrelated" & weights == "no weights") %>% 
  select(gene, kinship, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) %>% 
  
  ggplot(aes(x=`no kinship`,y=kinship))+
  geom_point(alpha=0.1) +
  geom_abline(intercept = 0, slope = 1, color="red")

rstr_subset_metric %>% 
  filter(subset == "unrelated" & weights == "no weights") %>% 
  select(gene, kinship, AIC) %>% 
  pivot_wider(names_from = kinship, values_from = AIC) %>% 
  
  ggplot(aes(x=`no kinship`,y=kinship))+
  geom_point(alpha=0.1)+
  # geom_abline(intercept = 0, slope = 1, color="red") +
  labs(title="AIC") +
  geom_abline(intercept = 8, slope = 1, color="green") +
  geom_abline(intercept = -8, slope = 1, color="green") 

library(GGally)



p1 <-rstr_subset_metric %>% 
  mutate(name = paste(subset,kinship,weights)) %>% 
  select(gene, name, AIC) %>% 
  pivot_wider(values_from = AIC, names_from = name) %>% 
  ggpairs(columns = 2:9) +
  geom_abline(intercept = 8, slope = 1, color="green") +
  geom_abline(intercept = -8, slope = 1, color="green") 

ggsave(p1, "temp.png", height=20, width=20)







