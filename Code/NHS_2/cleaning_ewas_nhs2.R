setwd("/udd/nhsom/EWAS/Validation/")

library(car)
#library(bestNormalize)
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
#library(adegenet)
library(glmnet)
library(fmsb)


categorize_smk <- function(new_data, old_data, method){

  new_data$smkdrcon <- old_data$smkdrcon

  if (method == 'harvard'){
    new_data$smkdrcon <- ifelse(new_data$smkdrcon == 0 | new_data$smkdrcon == 15 , 0,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon == 1 , 1,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon > 1 & new_data$smkdrcon < 9, 2,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon > 8 & new_data$smkdrcon < 11, 3,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon == 11 , 4,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon > 11 , 5,new_data$smkdrcon)
  }
  if (method == 1){
    new_data$smkdrcon <- ifelse(new_data$smkdrcon == 0 | new_data$smkdrcon == 15 , 0,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon == 1 , 1,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon > 1 & new_data$smkdrcon < 5, 2,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon > 4 & new_data$smkdrcon < 12, 3,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon ==12 , 4,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon >12 , 5,new_data$smkdrcon)
  }
  if (method == 2){
    new_data$smkdrcon <- ifelse(new_data$smkdrcon == 0 | new_data$smkdrcon == 15 , 0,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon == 1 , 1,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon > 1 & new_data$smkdrcon < 5 | new_data$smkdrcon == 8, 2,new_data$smkdrcon)

    new_data$smkdrcon <- ifelse(new_data$smkdrcon > 4 & new_data$smkdrcon < 9, 3,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon > 8 & new_data$smkdrcon < 12, 4,new_data$smkdrcon)
    new_data$smkdrcon <- ifelse(new_data$smkdrcon >11 , 5,new_data$smkdrcon)
  }
  return(new_data)
}


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



cont_ewas <- function(data, depvar, time1, time2,  covars , adjvars){
  c <- 0


  out_df <- data.frame(matrix(ncol = 10, nrow = length(covars)))#, row.names = b)


  row.names(out_df) <- covars
  colnames(out_df) <- c('coef','exp_coef', 'se_coef', 'z', 'p_val','LB','UB','stat_power','PH_pval', 'VIF')
  #beta_time <- list()
  # for power calculation
  alpha <- 0.05
  n <- length(unique(data$id))
  psi<- 0.001667384 #sum(data$chdcase==1)/(2*n-sum(data$chdcase==1))/2
  ua2<-qnorm(1-alpha/2)

  for (i in covars){


    c <- c+1

    baseform <- NULL
    doForm <- NULL
    baseform <- as.formula(sprintf('Surv(%s, %s, %s) ~ %s', time1, time2, depvar, i))
    doForm <- addToBase(baseform, adjvars)

    mod <- coxph(formula = doForm, data = data)
    # print(mod)


    out_df[i,]['coef'] <- summary(mod)$coefficients[1,][1]
    out_df[i,]['exp_coef'] <- summary(mod)$coefficients[1,][2]
    out_df[i,]['se_coef'] <- summary(mod)$coefficients[1,][3]
    out_df[i,]['z'] <- summary(mod)$coefficients[1,][4]
    out_df[i,]['p_val'] <- summary(mod)$coefficients[1,][5]
    out_df[i,]['LB'] <- summary(mod)$conf.int[1,][3]
    out_df[i,]['UB'] <- summary(mod)$conf.int[1,][4]


    zp <- cox.zph(mod)
    #beta_time[[i]] <- zp[1]
    out_df[i,]['PH_pval'] <- zp$table[1,3]

    #power
    ## overal incidence prop
    theta <- summary(mod)$coefficients[1,][2]
    formula_base <- as.formula(sprintf('%s ~ %s', i,paste(adjvars, collapse= "+") ))
    #formula_base <- as.formula(sprintf('%s ~ %s', i,paste(covars[! covars %in% c(i)], collapse= "+") ))
    sigma2<-var(data[i])
    res.lm<-lm(formula = formula_base, data = data)
    tt<-summary(res.lm)
    rho2<-tt$r.squared
    part2<-n*(log(theta))^2*sigma2*psi*(1-rho2)
    power<-pnorm(sqrt(part2)-ua2)
    res<-list(power=power, rho2=rho2, sigma2=sigma2, psi=psi)
    out_df[i,]['stat_power'] <- res$power
    #print(c)

    formula_base <- NULL
    formula_base <- as.formula(sprintf('%s ~ %s', i,paste(adjvars, collapse= "+") ))
    res.lm<-lm(formula = formula_base, data = data)
    out_df[i,]['VIF'] <- VIF(res.lm)
  }
  #output <- list("out_df" = out_df, "beta_time" = beta_time)
  return(out_df)
}



df <- read_sas('./DATA/nhs2_data.sas7bdat')    # 990363
df <- subset(df, start_time < stop_time)    # no change




df <- categorize_smk(df, df, 'harvard')
df$hbpcon <- ifelse(is.na(df$hbpcon), 0, df$hbpcon)
df$cholcon <- ifelse(is.na(df$cholcon), 0, df$cholcon)
df$aspicon <- ifelse(df$aspicon == 9, 0, df$aspicon)
df$veyncon <- ifelse(is.na(df$veyncon), 0, df$veyncon)
df$mvitcon <- ifelse(is.na(df$mvitcon) , 0, df$mvitcon)
df$mvitcon <- ifelse(df$mvitcon==2 , 1, 0)


indx <- apply(df, 2, function(x) any(is.nan(x)| is.na(x) | is.infinite(x)))
names(df)[indx]

df <- na.locf(df)

indx <- apply(df, 2, function(x) any(is.nan(x)| is.na(x) | is.infinite(x)))
names(df)[indx]



raw_adjustfor <- c('white', 'mifh', 'hbpfh', 'agecon', 'bmicon', 'smkdrcon', 'hbpcon',
                   'cholcon', 'aspicon', 'mvitcon', 'nhorcon', 'actcon', 'veyncon', 'calorconvconv')

covar <- c("wwinevconv",	"dressvconv",	"yogvconv",	"rcarvconv",	"liqvconv",	"rwinevconv",
           "donutvconv",	"hotdogvconv",	"whbrvconv",	"ajvconv",	"pnutvconv",
           "allprocvconv",	"hambuvconv",	"cervconv",	"sugbevvconv",	"cookvconv",
           "raisvconv",	"bmanvconv",	"branaconvconv",	"t161vconv",
           "addfat_liqvconv",	"addfat_solvconv",	"f161vconv",	"afatvconv",	"f180vconv",
           "mill_grnvconv",	"phytvconv",	"navconv",	"hydprovconv",	"satvconv",
           "uisorvconv",	"sum_whgcarbvconv",	"cerafvconv",	"t182vconv",	"uttocovconv",
           "f160vconv",	"fol_nvconv",	"ubT3vconv",	"mufa_plvconv",	"germnvconv",
           "uapigvconv",	"ubtocovconv",	"brannvconv",	"sevconv",	"f140vconv",
           "mn_nvconv",	"uaT3vconv",	"cholvconv",	"vcfolicvconv",	"t181vconv",
           "mnvconv",	"b6svconv",	"hemevconv",	"es2vconv",	"trnssvconv", "mufa_avconv", "alcovconv")

all_var <- c(raw_adjustfor, covar, 'id', 'chdcase', 'start_time', 'stop_time')

df <- df[, (names(df) %in% all_var)]

temp_var <- c(covar, 'calorconvconv')


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



add_var <- c('id', 'chdcase', 'start_time', 'stop_time', 'white', 'mifh', 'hbpfh', 'agecon', 'bmicon', 'smkdrcon', 'hbpcon',
             'cholcon', 'aspicon', 'mvitcon', 'nhorcon', 'actcon', 'veyncon')
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


write_sas(z_trans_df,'./DATA/cleaned_nhs2_data.sas7bdat')



adjustfor <- c('factor(white)', 'factor(mifh)', 'factor(hbpfh)', 'agecon', 'bmicon',
                   'factor(smkdrcon)', 'hbpcon',
                   'factor(cholcon)', 'factor(aspicon)', 'factor(mvitcon)',
                   'factor(nhorcon)', 'actcon', 'factor(veyncon)', 'calorconvconv')


no_cores <- detectCores() - 2
cl<-makeCluster(no_cores)
registerDoParallel(cl)
writeLines(c(""), "log.txt")

out <- foreach(covar = covar,
               .packages = c('survival', 'fmsb'),
               .combine = rbind)  %dopar% {
                 sink("log.txt", append = TRUE)
                 cat(paste("\n","Starting iteration",covar,"\n"))
                 sink() #end diversion of output
                 cont_ewas(z_trans_df, 'chdcase', 'start_time', 'stop_time', covar, adjustfor)
               }


out <- out[order(out$p_val),]

write.csv(out, './ewas_results_nhs2.csv')

stopCluster(cl)
stopImplicitCluster()
