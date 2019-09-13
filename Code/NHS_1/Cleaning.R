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


# function  >> to transform a continuous variable using Box-Cox transformation
# arguments >> data: the dataframe that we wish to transform the variable from. 
# arguments >> continuous_vars: a character indicating the continuous variable's name that we wish to transfrom
# return    >> a dataframe with one column containing the transformed values

boxcox_trans <- function(data, continuous_var)
{
  lambda <- boxCox(data[[continuous_var]]~1, family="yjPower", plotit = FALSE)
  lam_df <- data.frame(lambda$x, lambda$y)
  lambda <- lam_df[with(lam_df, order(-lam_df$lambda.y)),][1,1]
  out <- yjPower(data[[continuous_var]],lambda,jacobian.adjusted=TRUE)
  out <- data.frame(out)
  colnames(out) <- continuous_var
  return(out)
}

# read the data
data <-  read_sas("./data_sample.sas7bdat")

# this vector should be filled with exposure names.
covar <- c('fat_intake', 	'protein_intake'	'vitamin_intake'	'fruit_intake'	'meat_intake')  

# this vector should be filled with adjusting variables.
adjustfor <- c('physical_activity',	'calorie_intake')  


# put all required variable in a vector
all_var <- c(covar, adjustfor, 'id', 'chdcase', 'start_time', 'stop_time')
# id is a unique identifier for each participant
# chdcase is a binary variable showing whether a participant developed CHD or not
# start_time and stop_time indicate the beginning and end of an observation for a participant

data <- data[, (names(data) %in% all_var)]

# in our particular data, the missing values denoted by NA only appeared in the last observation 
# of each participant. The below function applies a last observation carried forward (LOCF) approach,
# a common method when working with time series data, to impute these missing values.
data <- na.locf(data, na.rm = FALSE)

# put continuous variables in a vector to apply Box-Cox transformation on.
continuous_var <- c('fat_intake', 	'protein_intake'	'vitamin_intake'	'fruit_intake'	'meat_intake', 
                     'physical_activity',	'calorie_intake')

# paraller computing for Box-Cox transformation
no_cores <- detectCores() - 2
cl<-makeCluster(no_cores)
registerDoParallel(cl)
writeLines(c(""), "log.txt")

trans_data <- foreach(continuous_var = continuous_var,
                    .packages = c('car'),
                    .combine = cbind)  %dopar% {
                      sink("log.txt", append = TRUE)
                      cat(paste("\n","Starting iteration",continuous_var,"\n"))
                      sink() #end diversion of output
                      boxcox_trans(data, continuous_var)
                    }

stopCluster(cl)
stopImplicitCluster()

add_var <- c('id', 'chdcase', 'start_time', 'stop_time')

for (i in add_var)
{
  trans_data[[i]] <- data[[i]]
}

# z-transformation
z_trans_data <- trans_data

for (i in continuous_var)
{
  z_trans_data[[i]] <- scale(z_trans_data[[i]])[,1]
}


write_sas(z_trans_data, './cleaned_data_sample.sas7bdat')
