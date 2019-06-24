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



addToBase <- function(base_formula, adjustingVariables) {
  for (var in adjustingVariables) {
    base_formula <- update.formula( base_formula, as.formula(sprintf('~ . + %s', var)) )
  }
  return(base_formula) }

shuff_cwas <- function(data, depvar, time1, time2,  covars , adjvars){

  out_df <- data.frame(matrix(ncol = length(covars), nrow = 1))
  colnames(out_df) <- covars

  for (i in covars){
    baseform <- NULL
    doForm <- NULL
    baseform <- as.formula(sprintf('Surv(%s, %s, %s) ~ %s', time1, time2, depvar, i))
    doForm <- addToBase(baseform, adjvars)
    mod <- coxph(formula = doForm, data = data)
    out_df[1,i] <- summary(mod)$coefficients[1,][5]
  }
  return(out_df)
}


data <- read_sas("./data.sas7bdat")

covar <- c()     # this vector should be filled with exposure names.
adjustfor <- c()   # this vector should be filled with adjusitng variables.



baseform <- NULL
doForm <- NULL
baseform <- as.formula(sprintf('Surv(%s, %s, %s) ~ %s', 'start_time', 'stop_time', 'chdcase', 'agecon'))
doForm <- addToBase(baseform, adjustfor)
mod <- coxph(formula = doForm, data = data)
num_cases <- 2774            # this is the number of cases in the dataset.


pred <- survival:::predict.coxph(mod, type='expected')
pred <- pred/num_cases
id_idx <- data.frame(data$id)
id_idx$pred <- pred
id_idx$pred2 <- cumsum(id_idx$pred)
id_idx$pred20 <- id_idx$pred2 - id_idx$pred
id_idx$chdcase <- 0

adjustfor <- c()

n_permutation <- 1000
perm_df <- data.frame()


log <- "log.txt"

for (n in 1:n_permutation){

  id_idx$chdcase <- 0
  data$new_chdcase <- 0

  while (sum(id_idx$chdcase) < num_cases){
    rand_p <- runif(1, min=0, max=1)
    id_idx$chdcase[between(rand_p, id_idx$pred20, id_idx$pred2)] <- 1
  }

  data$new_chdcase <- id_idx$chdcase
  data$new_chdcase <- ave(data$new_chdcase, data$id, FUN = function(x) replace(x, which(cumsum(x)>=1)[-1], NA))
  df <- subset(data, ! is.na(data$new_chdcase))
                          
  # parallel computing 
  no_cores <- detectCores() - 2
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  writeLines(c(""), log)
  print(paste('iter',n))

  shuff_pval <- foreach(covar = covar,
                        .packages = c('survival'),
                        .combine = cbind)  %dopar% {
                          sink(log, append = TRUE)
                          cat(paste("\n","Starting iteration",covar,"\n"))
                          sink()
                          shuff_cwas(df, 'new_chdcase', 'start_time', 'stop_time', covar, adjustfor)
                        }

  perm_df <- rbind(shuff_pval, perm_df)
  path <- './permutations.csv'
  write.csv(perm_df, path)

  stopCluster(cl)
  stopImplicitCluster()

}
