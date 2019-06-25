library(car)
library(bestNormalize)
library(data.table)
library(datasets)
library(devtools)
library(doParallel)
library(foreach)
library(haven)
library(MASS)
library(parallel)
library(survival)
library(zoo)
library(adegenet)
library(glmnet)
library(fmsb)
library(plyr)
library(corrplot)
library(Hmisc)
library(corrplot)

data <- read_sas('./data.sas7bdat')
sig_covar <- c(...)     # this vector should be filled with statisitcally significant variables found by EWAS.


cor_res <- rcorr(as.matrix(data), type = 'spearman')

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

cor_res_flatten <- flattenCorrMatrix(cor_res$r, cor_res$P)
write.csv(cor_res_flatten, './cor_res_flatten.csv')


out_df <- data.frame(matrix(ncol = 3, nrow = 1378))
colnames(out_df) <- c('row', 'column', 'perm_cor')

# permutation
c <- 1
cc <- 1
for (i in sig_covar[1:52]){
  cc <- cc+1
  for (j in sig_covar[cc:53]){
    if (i != j){
      temp <- as.matrix(sub_merge[c(i,j)])
      temp[,2] <- sample(temp[,2])
      out_df[c,]['row'] <- i
      out_df[c,]['column'] <- j
      out_df[c,]['perm_cor'] <- rcorr(temp, type = 'spearman')$r[2]
      c <- c+1
      print(c)
    }
  }
}


out_df <- out_df[
  with(out_df, order(row, column)),
  ]


cor_res_flatten <- cor_res_flatten[
  with(cor_res_flatten, order(row, column)),
  ]


write.csv(out_df, './perm_cor.csv')

cor_merge_1 <- merge(out_df,cor_res_flatten,  by = c('row', 'column'))
cor_res_flatten_2 <- cor_res_flatten
names(cor_res_flatten_2) <- c('column', 'row', 'cor', 'p' )
cor_merge_2 <- merge(out_df,cor_res_flatten_2,  by = c('row', 'column'))
cor_merge_2 <- rbind(cor_merge_2,cor_merge_1)

write.csv(cor_merge_2, './permuted_correlation.csv')


# visualization
thre <- 5.47E-05
corrplot(cor_res$r,  order="hclust", tl.cex = 0.5, tl.col = 'black', type="upper",
         p.mat = cor_res$P, sig.level = thre, insig = "pch", pch.cex = 0.5)
