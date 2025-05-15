
##1. lmm: function to perform the LMM of warming effects on 
# each variable (each column)
library(lme4)
library(car)
lmm <- function(j, dat = divindex, treat = treatused) {
  # dat is the dataframe to store the variables
  # treat is the treatment info, including the year, block, warm/control treatment
  # j is the column id
  message("Now j=",j," in ",ncol(dat),". ",date())
  if (length(unique(dat[,j]))<4){
    result<-rep(NA,8)
  } else{
    div<-data.frame(divtest=dat[,j],treat)
    fm1<-lmer(divtest~warm+(1|year)+(1|block),data=div)
    
    presult<-car::Anova(fm1,type=2)
    coefs<-coef(summary(fm1))[ , "Estimate"]  ##four coefs
    names(coefs)<-paste0(names(coefs),".mean")
    
    SEvalues<-coef(summary(fm1))[ , "Std. Error"] ##standard errors
    names(SEvalues)<-paste0(names(SEvalues),".se")
    
    tvalues<-coef(summary(fm1))[ , "t value"] ##t values
    names(tvalues)<-paste0(names(tvalues),".t")
    
    chisqP<-c(presult[,1],presult[,3])
    names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
    
    result<-c(coefs,tvalues,SEvalues,chisqP)
  }
  if (length(result) < 8) {
    result<-rep(NA,8)
  }
  result 
}


##2. lmm2: function to perform the LMM of warming, precipitation, and clipping effects on 
# each variable (each column)
# suitable for the data of year 2019, which contain all 48 plots

lmm2 <- function(j, dat = divindex, treat = treatused) {
  # paramters as in 'lmm'
  # make sure the samples names matched in dat and treat
  message("Now j=",j," in ",ncol(dat),". ",date())
  if (length(unique(dat[,j]))<3){
    result<-rep(NA,38)
  } else {
    div<-data.frame(divtest=dat[,j],treat)
    
    fm1<-lmer(divtest~warm*precip*clip+(1|block),data=div)
    
    presult<-car::Anova(fm1,type=2)
    coefs<-coef(summary(fm1))[ , "Estimate"]  ##four coefs
    names(coefs)<-paste0(names(coefs),".mean")
    
    SEvalues<-coef(summary(fm1))[ , "Std. Error"] ##standard errors
    names(SEvalues)<-paste0(names(SEvalues),".se")
    
    tvalues<-coef(summary(fm1))[ , "t value"] ##t values
    names(tvalues)<-paste0(names(tvalues),".t")
    
    chisqP<-c(presult[,1],presult[,3])
    names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
    result<-c(coefs,tvalues,SEvalues,chisqP)}
  if (length(result) < 38) {
    result<-rep(NA,38)
  }
  result 
}


## 3. warm_eff: function to perform LMMs based on pure warming and control samples
warm_eff <- function(dt, env) {
  # dt is the data frame of variables
  # env is the metadata of all samples
  dat_all <- dt
  
  envCW<-env[(env$precip==0 & env$clip==0) | env$year==2009,] # select pure warming and control samples
  
  cov_sumCW<-dat_all[match(row.names(envCW),row.names(dat_all)),]
  message("Number of samples not matched: ", sum(row.names(envCW) != row.names(cov_sumCW)))

  treatused<-envCW
  treatused$year<-treatused$year-2009

  divindex<-scale(cov_sumCW)  # scale each variable to standard normal distribution
  
  divs1 <- sapply(1:ncol(divindex), function(x) tryCatch(lmm(j = x, dat =divindex, treat = treatused), error=function(e) rep(NA,8)))
  colnames(divs1)<-colnames(divindex)
  return(divs1)  # return the LMM effect size, p, t, etc.
}

## 3.2 warm_eff_unscaled: function to perform LMMs based on pure warming and control samples
# originnal rather than scaled variables

warm_eff_unscaled <- function(dt, env) {
  # dt is the data frame of variables
  # env is the metadata of all samples
  dat_all <- dt
  
  envCW<-env[(env$precip==0 & env$clip==0) | env$year==2009,] # select pure warming and control samples
  
  cov_sumCW<-dat_all[match(row.names(envCW),row.names(dat_all)),]
  message("Number of samples not matched: ", sum(row.names(envCW) != row.names(cov_sumCW)))
  
  treatused<-envCW
  treatused$year<-treatused$year-2009
  
  divindex<- cov_sumCW  
  
  divs1 <- sapply(1:ncol(divindex), function(x) tryCatch(lmm(j = x, dat =divindex, treat = treatused), error=function(e) rep(NA,8)))
  colnames(divs1)<-colnames(divindex)
  return(divs1)  # return the LMM effect size, p, t, etc.
}


## 4. treat_eff: function to perform LMMs based on 2019 samples (48 plots)
treat_eff <- function(dt, env) {
  # dt is the data frame of variables
  # env is the metadata of all samples
  dat_all <- dt
  env2019<-env[env$year==2019,]
  cov_sum2019<-dat_all[match(row.names(env2019),row.names(dat_all)),]
  message("Number of samples not matched: ", sum(row.names(env2019) != row.names(cov_sum2019)))
  
  treatused<-env2019
  divindex<-scale(cov_sum2019)
  
  divs1 <- sapply(1:ncol(divindex), function(x) tryCatch(lmm2(j = x, dat =divindex, treat = treatused), error=function(e) rep(NA, 38)))
  colnames(divs1)<-colnames(divindex)
  return(divs1)  # return the LMM effect size, p, t, etc.
}
## 5. warm_eff_flt: function as 3, but exclude potential outliers (samples with HQ base > 50GB)

warm_eff_flt <- function(dt, env) {
  # dt is the data frame of variables
  # env is the metadata of all samples
  dat_all <- dt
  
  envCW<-env[(env$precip==0 & env$clip==0) | env$year==2009,] # select pure warming and control samples
  
  cov_sumCW<-dat_all[match(row.names(envCW),row.names(dat_all)),]
  # exclude samples with large HQ size
  envCW2 <- envCW[envCW$HQ.base.number..Gb. < 50, ]
  cov_sumCW2 <- cov_sumCW[envCW$HQ.base.number..Gb. < 50, ]
  
  message("Number of samples not matched: ", sum(row.names(envCW2) != row.names(cov_sumCW2)))
  
  treatused<-envCW2
  treatused$year<-treatused$year-2009
  
  divindex<-scale(cov_sumCW2)  # scale each variable to standard normal distribution
  
  divs1 <- sapply(1:ncol(divindex), function(x) tryCatch(lmm(j = x, dat =divindex, treat = treatused), error=function(e) rep(NA,8)))
  colnames(divs1)<-colnames(divindex)
  return(divs1)  # return the LMM effect size, p, t, etc.
}

## 6. treat_eff_flt: function as 4, but exclude potential outliers (samples with HQ base > 50GB)

treat_eff_flt  <- function(dt, env) {
  # dt is the data frame of variables
  # env is the metadata of all samples
  dat_all <- dt
  env2019<-env[env$year==2019,]
  cov_sum2019<-dat_all[match(row.names(env2019),row.names(dat_all)),]
  # exclude samples with large HQ size
  env2 <- env2019[env2019$HQ.base.number..Gb. < 50, ]
  cov_sum2 <- cov_sum2019[env2019$HQ.base.number..Gb. < 50, ]
  
  
  message("Number of samples not matched: ", sum(row.names(env2) != row.names(cov_sum2)))
  
  treatused<-env2
  divindex<-scale(cov_sum2)
  
  divs1 <- sapply(1:ncol(divindex), function(x) tryCatch(lmm2(j = x, dat =divindex, treat = treatused), error=function(e) rep(NA, 38)))
  colnames(divs1)<-colnames(divindex)
  return(divs1)  # return the LMM effect size, p, t, etc.
}

## 7. treat_eff_all: function perform LMMs based on all 2016-2020 samples  (552 samples)
lmm3 <- function(j, dat = divindex, treat = treatused) {
  # paramters as in 'lmm'
  # make sure the samples names matched in dat and treat
  message("Now j=",j," in ",ncol(dat),". ",date())
  if (length(unique(dat[,j]))<3){
    result<-rep(NA,38)
  } else {
    div<-data.frame(divtest=dat[,j],treat)
    
    fm1<-lmer(divtest~warm*precip*clip+(1|year)+(1|block),data=div)
    presult<-car::Anova(fm1,type=2)
    coefs<-coef(summary(fm1))[ , "Estimate"]  ##four coefs
    names(coefs)<-paste0(names(coefs),".mean")
    
    SEvalues<-coef(summary(fm1))[ , "Std. Error"] ##standard errors
    names(SEvalues)<-paste0(names(SEvalues),".se")
    
    tvalues<-coef(summary(fm1))[ , "t value"] ##t values
    names(tvalues)<-paste0(names(tvalues),".t")
    
    chisqP<-c(presult[,1],presult[,3])
    names(chisqP)<-c(paste0(row.names(presult),".chisq"),paste0(row.names(presult),".P"))
    result<-c(coefs,tvalues,SEvalues,chisqP)}
  if (length(result) < 38) {
    result<-rep(NA,38)
  }
  result 
}


treat_eff_all <- function(dt, env) {
  # dt is the data frame of variables
  # env is the metadata of all 552 samples
  dat_all <- dt
  cov_sum<-dat_all[match(row.names(env),row.names(dat_all)),]
  message("Number of samples not matched: ", sum(row.names(env) != row.names(cov_sum)))
  
  treatused<-env
  treatused$year<-treatused$year-2009  
  divindex<-scale(cov_sum)
  
  divs1 <- sapply(1:ncol(divindex), function(x) tryCatch(lmm3(j = x, dat =divindex, treat = treatused), error=function(e) rep(NA, 38)))
  colnames(divs1)<-colnames(divindex)
  return(divs1)  # return the LMM effect size, p, t, etc.
}