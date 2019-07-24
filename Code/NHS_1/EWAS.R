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

# this function associates an exposure with disease status, while controlling for adjusting variables, using Cox model.
cont_ewas <- function(data, depvar, time1, time2,  covars , adjvars, alpha, psi){
  c <- 0
  out_df <- data.frame(matrix(ncol = 10, nrow = length(covars)))
  row.names(out_df) <- covars
  colnames(out_df) <- c('coef','exp_coef', 'se_coef', 'z', 'p_val','LB','UB','stat_power','PH_pval', 'VIF')

  for (i in covars){

    c <- c+1

    baseform <- NULL
    doForm <- NULL
    baseform <- as.formula(sprintf('Surv(%s, %s, %s) ~ %s', time1, time2, depvar, i))
    doForm <- addToBase(baseform, adjvars)
    mod <- coxph(formula = doForm, data = data)


    out_df[i,]['coef'] <- summary(mod)$coefficients[1,][1]
    out_df[i,]['exp_coef'] <- summary(mod)$coefficients[1,][2]
    out_df[i,]['se_coef'] <- summary(mod)$coefficients[1,][3]
    out_df[i,]['z'] <- summary(mod)$coefficients[1,][4]
    out_df[i,]['p_val'] <- summary(mod)$coefficients[1,][5]
    out_df[i,]['LB'] <- summary(mod)$conf.int[1,][3]
    out_df[i,]['UB'] <- summary(mod)$conf.int[1,][4]

    # checking for the proportional hazard assumption
    zp <- cox.zph(mod)
    out_df[i,]['PH_pval'] <- zp$table[1,3]

    # power calculation
    n <- length(unique(data$id))
    ua2<-qnorm(1-alpha/2)
    theta <- summary(mod)$coefficients[1,][2]
    formula_base <- as.formula(sprintf('%s ~ %s', i,paste(adjvars, collapse= "+") ))
    sigma2<-var(data[i], na.rm = TRUE)
    res.lm<-lm(formula = formula_base, data = data)
    tt<-summary(res.lm)
    rho2<-tt$r.squared
    part2<-n*(log(theta))^2*sigma2*psi*(1-rho2)
    power<-pnorm(sqrt(part2)-ua2)
    res<-list(power=power, rho2=rho2, sigma2=sigma2, psi=psi)
    out_df[i,]['stat_power'] <- res$power

    # VIF calculation
    formula_base <- NULL
    formula_base <- as.formula(sprintf('%s ~ %s', i,paste(adjvars, collapse= "+") ))
    res.lm<-lm(formula = formula_base, data = data)
    out_df[i,]['VIF'] <- VIF(res.lm)
  }
  return(out_df)
}


# read data
data_path <- './data.sas7bdat'
data <- read_sas(data_path)

covar <- c(...)     # this vector should be filled with exposure names.
adjustfor <- c(...)  # this vector should be filled with adjusting variables, such as age, gender, race, etc..

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
