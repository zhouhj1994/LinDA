# LinDA
Linear model for differential abundance analysis. 

Reference: Huijuan Zhou, Kejun He, Jun Chen, and Xianyang Zhang. (2021). LinDA: Linear Models for Differential Abundance Analysis of Microbiome Compositional Data.

LinDA implements a simple, robust and highly scalable approach to tackle the compositional effects in differential abundance analysis. 
It fits linear regression models on the centered log2-ratio transformed data, identifies a bias term due to the transformation
and compositional effect, and corrects the bias using the mode of the regression coefficients. It could fit mixed-effect models. 

We have integrated the procedure into our CRAN package "MicrobiomeStat" as the "linda" function (https://CRAN.R-project.org/package=MicrobiomeStat).

## Installation
```r
# install.packages(c("modeest", "lmerTest", "foreach", "parallel", "ggplot2", "ggrepel"))
# install.packages("devtools")
devtools::install_github("zhouhj1994/LinDA")
```
## An Example
Apply LinDA to a dataset from the study of the smoking effect on the human upper respiratory tract. (Reference: Charlson et al. (2010). 
Disordered microbial communities in the upper respiratory tract of cigarette smokers.)

```r
#install package "phyloseq" for importing "smokers" dataset
ind <- smokers$meta$AIRWAYSITE == 'Throat'
otu.tab <- as.data.frame(smokers$otu[, ind])
meta <- cbind.data.frame(Smoke = factor(smokers$meta$SMOKER[ind]),
                         Sex = factor(smokers$meta$SEX[ind]),
                         Site = factor(smokers$meta$SIDEOFBODY[ind]),
                         SubjectID = factor(smokers$meta$HOST_SUBJECT_ID[ind]))
ind1 <- which(meta$Site == 'Left')
res.left <- linda(otu.tab[, ind1], meta[ind1, ], formula = '~Smoke+Sex', alpha = 0.1,
                  prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
ind2 <- which(meta$Site == 'Right')
res.right <- linda(otu.tab[, ind2], meta[ind2, ], formula = '~Smoke+Sex', alpha = 0.1,
                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
rownames(res.left$output[[1]])[which(res.left$output[[1]]$reject)]
rownames(res.right$output[[1]])[which(res.right$output[[1]]$reject)]

linda.obj <- linda(otu.tab, meta, formula = '~Smoke+Sex+(1|SubjectID)', alpha = 0.1,
                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
linda.plot(linda.obj, c('Smokey', 'Sexmale'), 
           titles = c('Smoke: n v.s. y', 'Sex: female v.s. male'), alpha = 0.1, lfc.cut = 1,
           legend = TRUE, directory = NULL, width = 11, height = 8)

L <- matrix(c(0, 1, 0, 0, 0, 1), nrow = 2, byrow = TRUE)
#L <- matrix(c(0, 1, 0), nrow = 1, byrow = TRUE)
linda.wald.test(linda.obj, L, 'LMM', alpha = 0.1)
```
