library(tidyverse)

#Load kinship matrix
load("data/RSTR_RNAseq_kinship.RData")
#Replace diagonal of matrix (self kinship comparisons) with 0 as these are ~1
# and generally the highest kinship value for each individual
kin.matrix0 <- kin.matrix
diag(kin.matrix0) <- 0

#Look for related indiv
kin.summ <- as.data.frame(kin.matrix0) %>% 
  #Move rownames into df
  rownames_to_column() %>% 
  #Move column names into df
  pivot_longer(-rowname, names_to = "ptID") %>% 
  #Calculate max kinship value per person
  group_by(ptID) %>% 
  summarise(maxK = max(value)) %>% 
  #Create groups for individuals with and without related indiv in the data set
  mutate(group = ifelse(maxK >= 0.125, 'related',"unrelated"))

#Get FULLID of indiv to filter by
## Indiv with at least 1 related person in data set
kin.related <- kin.summ %>% 
  filter(group == "related") %>% 
  select(name) %>% unlist(use.names = FALSE)
## Indiv with no related persons in data set
kin.unrelated <- kin.summ %>% 
  filter(group == "unrelated") %>% 
  select(name) %>% unlist(use.names = FALSE)

#Make unrelated matrix
## Note we go back to the orig kinship matrix with real diagonal values
kin.matrix.unrelated <- as.data.frame(kin.matrix) %>% 
  rownames_to_column() %>% 
  filter(rowname %in% kin.unrelated) %>% 
  select(rowname, all_of(kin.unrelated)) %>% 
  column_to_rownames() %>% 
  as.matrix()
  
##check order because a kinship matrix must be square with rows and columns
## in the same order
identical(rownames(kin.matrix.unrelated), colnames(kin.matrix.unrelated))

#Make related matrix
kin.matrix.related <- as.data.frame(kin.matrix) %>% 
  rownames_to_column() %>% 
  filter(rowname %in% kin.related) %>% 
  select(rowname, all_of(kin.related)) %>% 
  column_to_rownames() %>% 
  as.matrix()

##check order
identical(rownames(kin.matrix.related), colnames(kin.matrix.related))


# Check if kinship subsets result in uneven RSTR / LTBI groups
# These look good for linear modeling
load("data/RSTR_kimma_voom.RData")
datV$targets %>% 
  filter(ptID %in% kin.unrelated) %>% 
  count(Sample_Group)

datV$targets %>% 
  filter(ptID %in% kin.related) %>% 
  count(Sample_Group)
