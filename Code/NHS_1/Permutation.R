# load required libraries

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



# functoin  >> to extend the base formula with variables we are adjusting for
# @arguments >> base_formula: a formula object with the dependent variable as a function of the independent variable
# @arguments >> adjustingVariables: a character vector of column names from the data.frame to adjust for
# @return    >> an updated formula object representing a base formula with adjustment variables
 
addToBase <- function(base_formula, adjustingVariables) 
{
  for (var in adjustingVariables) {
    base_formula <- update.formula( base_formula, as.formula(sprintf('~ . + %s', var)) )
  }
  return(base_formula) 
}


# function   >> to recalculate the cox model for the shuffeled data and collect the null p-value
# @arguments >> data: a dataframe containing the cohort data
# @arguments >> depvar: dependent variable in the survival analysis
# @arguments >> time1: the variable showing the start time of an observation
# @arguments >> time2: the variable showing the end time of an observation
# @arguments >> covars: exposure of interest that we want to find its effect on a disease risk
# @arguments >> adjvars: a vector containing the variables that we want to adjust our analysis for
# @return    >> a dataframe containing the name of the exposure of interest and the null p-value

shuff_ewas <- function(data, depvar, time1, time2,  covar , adjvars)
{

  out_df <- data.frame(matrix(ncol = 1, nrow = 1))
  colnames(out_df) <- covar

  baseform <- NULL
  doForm <- NULL
  baseform <- as.formula(sprintf('Surv(%s, %s, %s) ~ %s', time1, time2, depvar, covar))
  doForm <- addToBase(baseform, adjvars)
  mod <- coxph(formula = doForm, data = data)
  out_df[1,covar] <- summary(mod)$coefficients[1,][5]
  
  return(out_df)
}


# read data
data <- read_sas('cleaned_data_sample.sas7bdat')

# this vector should be filled with exposure names.
covar <- c('fat_intake', 	'protein_intake',	'vitamin_intake',	'fruit_intake',	'meat_intake')  

# this vector should be filled with adjusting variables.
adjustfor <- c('physical_activity',	'calorie_intake')  


# associate the adjusting variables with the time of developing a disease of interest
baseform <- NULL
doForm <- NULL
baseform <- as.formula(sprintf('Surv(%s, %s, %s) ~ %s + %s', 'start_time', 'stop_time',
                               'chdcase', 'physical_activity',	'calorie_intake'))
mod <- coxph(formula = baseform, data = data)

# this is the number of cases in the dataset.
num_cases <- sum(data$chdcase)            

# predict the failure probability for each participant at each time stamp
failure_prob <- survival:::predict.coxph(mod, type='expected')
failure_prob <- failure_prob/num_cases

# below, we generate a dataframe as following:

# | id  |  failure_prob  |  failure_prob_cumsum  |  failure_prob_bound  |   
# -----------------------------------------------------------------------
# |  1  |     0.10       |          0.10         |          0           |    
# |  1  |     0.05       |          0.15         |         0.10 `       |     
# |  1  |     0.15       |          0.30         |         0.15         |  
# |  2  |     0.20       |          0.50         |         0.30         |   
# |  2  |     0.15       |          0.65         |         0.50         |   
# |  2  |     0.25       |          0.90         |         0.65         |     
# |  2  |     0.10       |          1.00         |         0.90         |     


id_idx <- data.frame(data$id)
id_idx$failure_prob <- failure_prob
id_idx$failure_prob_cumsum <- cumsum(id_idx$failure_prob)
id_idx$failure_prob_bound <- id_idx$failure_prob_cumsum - id_idx$failure_prob

# specify the number of permutations
n_permutation <- 1000


perm_df <- data.frame()
log <- "log.txt"

for (n in 1:n_permutation){

  id_idx$chdcase <- 0
  data$new_chdcase <- 0

  while (sum(id_idx$chdcase) < num_cases){
    rand_p <- runif(1, min=0, max=1)
    id_idx$chdcase[between(rand_p, id_idx$failure_prob_bound, id_idx$failure_prob_cumsum)] <- 1
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
                          shuff_ewas(df, 'new_chdcase', 'start_time', 'stop_time', covar, adjustfor)
                        }

  perm_df <- rbind(shuff_pval, perm_df)
  path <- './permutations.csv'
  write.csv(perm_df, path)

  stopCluster(cl)
  stopImplicitCluster()

}
