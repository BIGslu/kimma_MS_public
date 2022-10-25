# kimma_MS_public

### kimma: flexible linear mixed effects modeling with kinship covariance for RNA-seq data

Dill-McFarland KA, Mitchell K, Batchu S, Segnitz RM, Benson B, Janczyk T, Cox MS, Mayanja-Kizza HA, Boom WH, Benchek P, Stein C, Hawn TR, Altman MC

### Abstract 

Motivation: The identification of differentially expressed genes (DEGs) from transcriptomic datasets is a major avenue of research across diverse disciplines. However, current bioinformatic tools do not support covariance matrices in DEG modeling. Here, we introduce kimma (Kinship In Mixed Model Analysis), an open-source R package for flexible linear mixed effects modeling including covariates, weights, random effects, covariance matrices, and fit metrics.

Results: In simulated datasets, kimma detects DEGs with similar specificity, sensitivity, and computational time as limma unpaired and dream paired models. Unlike other software, kimma supports covariance matrices as well as fit metrics like AIC. Utilizing genetic kinship covariance, kimma revealed that kinship impacts model fit and DEG detection in a related cohort. Thus, kimma equals or outcompetes current DEG pipelines in sensitivity, computational time, and model complexity. 

Availability and implementation: Kimma is freely available on GitHub https://github.com/BIGslu/kimma with an instructional vignette at https://bigslu.github.io/kimma_vignette/kimma_vignette.html 

Supplemental information: Scripts related to this manuscript can be found at https://github.com/BIGslu/kimma_MS_public 
![image](https://user-images.githubusercontent.com/21342185/197858718-21ad644e-ed45-4bda-b729-4b891201368f.png)
