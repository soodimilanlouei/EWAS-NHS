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


# function  >> to associate an exposure with disease status, while controlling for adjusting variables, using Cox model
# @arguments >> data: a dataframe containing the cohort data
# @arguments >> depvar: dependent variable in the survival analysis
# @arguments >> time1: the variable showing the start time of an observation
# @arguments >> time2: the variable showing the end time of an observation
# @arguments >> covars: exposure of interest that we want to find its effect on a disease risk
# @arguments >> adjvars: a vector containing the variables that we want to adjust our analysis for
# @return    >> a dataframe containing the result of Cox model, proportional hazard assumption, and VIF for an exposure of interest

cont_ewas <- function(data, depvar, time1, time2,  covars , adjvars)
{
  # initiate a dataframe to store the results in
  out_df <- data.frame(matrix(ncol = 9, nrow = length(covars)))
  row.names(out_df) <- covars
  colnames(out_df) <- c('coef','exp_coef', 'se_coef', 'z', 'p_val','LB','UB','PH_pval', 'VIF')

  for (i in covars)
  {
    baseform <- NULL
    doForm <- NULL
    baseform <- as.formula(sprintf('Surv(%s, %s, %s) ~ %s', time1, time2, depvar, i))
    doForm <- addToBase(baseform, adjvars)
    mod <- coxph(formula = doForm, data = data)

    # collect the estimated effect size by Cox model, hazard ratio, standard deviation of the estimated effect,
    # z-statistic, p-value, lower and upper bound estimated for an effect.
    out_df[i,]['coef'] <- summary(mod)$coefficients[1,][1]
    out_df[i,]['exp_coef'] <- summary(mod)$coefficients[1,][2]
    out_df[i,]['se_coef'] <- summary(mod)$coefficients[1,][3]
    out_df[i,]['z'] <- summary(mod)$coefficients[1,][4]
    out_df[i,]['p_val'] <- summary(mod)$coefficients[1,][5]
    out_df[i,]['LB'] <- summary(mod)$conf.int[1,][3]
    out_df[i,]['UB'] <- summary(mod)$conf.int[1,][4]

    # check for the proportional hazard assumption by assessing the correlation between 
    # the Schoenfeld residuals and survival time
    zp <- cox.zph(mod)
    out_df[i,]['PH_pval'] <- zp$table[1,3]

    # check for multi-collinearity
    # VIF calculation
    formula_base <- NULL
    formula_base <- as.formula(sprintf('%s ~ %s', i,paste(adjvars, collapse= "+") ))
    res.lm<-lm(formula = formula_base, data = data)
    out_df[i,]['VIF'] <- VIF(res.lm)
  }
  return(out_df)
}


# read data
data <- read_sas('cleaned_data_sample.sas7bdat')

# this vector should be filled with exposure names.
covar <- c('fat_intake', 	'protein_intake'	'vitamin_intake'	'fruit_intake'	'meat_intake')  

# this vector should be filled with adjusting variables.
adjustfor <- c('physical_activity',	'calorie_intake')  

# parallel computing
no_cores <- detectCores() - 2
cl<-makeCluster(no_cores)
registerDoParallel(cl)
writeLines(c(""), "log.txt")

out <- foreach(covar = covar,
               .packages = c('survival', 'fmsb'),
               .combine = rbind)  %dopar% {
                 sink("log.txt", append = TRUE)
                 cat(paste("\n","Starting iteration",covar,"\n"))
                 sink()
                 cont_ewas(data, 'chdcase', 'start_time', 'stop_time', covar, adjustfor)
               }


out <- out[order(out$p_val),]
write.csv(out, './EWAS_results.csv')

# stop parallel computing
stopCluster(cl)
stopImplicitCluster()
