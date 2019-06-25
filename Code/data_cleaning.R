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

boxcox_trans <- function(df, temp_var){
  for (i in temp_var){
    lambda <- boxCox(df[[i]]~1, family="yjPower", plotit = FALSE)
    lam_df <- data.frame(lambda$x, lambda$y)
    lambda <- lam_df[with(lam_df, order(-lam_df$lambda.y)),][1,1]
    out <- yjPower(df[[i]],lambda,jacobian.adjusted=TRUE)
  }
  out <- data.frame(out)#, ncol=1)
  colnames(out) <- i
  return(out)
}


df <-  read_sas("./data.sas7bdat")



covar <- c()
adjustfor <- c()

all_var <- c(covar, adjustfor, 'id', 'chdcase', 'start_time', 'stop_time')
df <- df[, (names(df) %in% all_var)]


df <- na.locf(df, na.rm = FALSE)


temp_var <- c(covar,  'calorconv')



no_cores <- detectCores() - 2
cl<-makeCluster(no_cores)
registerDoParallel(cl)
writeLines(c(""), "log.txt")

trans_df <- foreach(temp_var = temp_var,
                    .packages = c('car'),
                    .combine = cbind)  %dopar% {
                      sink("log.txt", append = TRUE)
                      cat(paste("\n","Starting iteration",temp_var,"\n"))
                      sink() #end diversion of output
                      boxcox_trans(df, temp_var)
                    }

stopCluster(cl)
stopImplicitCluster()

add_var <- c('id', 'chdcase', 'start_time', 'stop_time', 'agecon', ...)

for (i in add_var){
  trans_df[[i]] <- df[[i]]
}


z_trans_df <- trans_df
c <- 1
for (i in temp_var){
  z_trans_df[[i]] <- scale(z_trans_df[[i]])[,1]
  print(c)
  c <- c+1
}


write_sas(z_trans_df, './cleaned_data.sas7bdat')
