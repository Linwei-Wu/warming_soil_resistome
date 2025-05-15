# 0. load data
setwd("/Users/linwei/Dropbox/warming_ARG/code_source_data/data")
# read-based coverage depth table
reads <- read.csv("Abundance_different cats.csv",row.names = 1)
# geochip summary table 
geochip <- read.csv("/Users/linwei/Dropbox/warming_ARG/0_raw_data/geochip/geochip_args_different cats.csv",row.names = 1)
# meta data
env <- read.csv("summary of environmental data2009-2020.csv",row.names = 1)

# 1. LMMs on the variables
source("../LMM_functions.R")
# 1.1 read-based allele coverage
# normalize coverages by HQ size 
reads <- reads[match(row.names(env), row.names(reads)),]
sum(row.names(env) == row.names(reads))
reads_pergb <- reads/env$HQ.base.number..Gb.
warm_eff_reads <- warm_eff(dt = reads_pergb, env = env)
write.csv(warm_eff_reads,"reads_LMM_warm.csv")
treat_eff_reads <- treat_eff(dt = reads_pergb, env = env)
write.csv(treat_eff_reads,"reads_LMM_treat2019.csv")
#warm_eff_reads_flt <- warm_eff_flt(dt = reads_pergb, env = env)
#write.csv(warm_eff_reads_flt,"reads_LMM_warm_rmoutlr.csv")
#treat_eff_reads_flt   <- treat_eff_flt(dt = reads_pergb, env = env)
#write.csv(treat_eff_reads_flt,"reads_LMM_treat2019_rmoutlr.csv")


# 1.2 Geochip
warm_eff_geochip <- warm_eff(dt = geochip, env = env)
write.csv(warm_eff_geochip,"geochip_LMM_warm.csv")
treat_eff_geochip <- treat_eff(dt = geochip, env = env)
write.csv(treat_eff_geochip,"geochip_LMM_treat2019.csv")
# test the treatment effects on all 552 samples of geochip data
geochip2 <- read.csv("geochip_args_different cats_all552.csv",row.names = 1)
env2 <- read.csv("env_all552samples.csv",row.names = 2)
treat_eff_geochip_all552 <- treat_eff_all(dt = geochip2, env = env2)
write.csv(treat_eff_geochip_all552,"geochip_LMM_treat_all552.csv")

