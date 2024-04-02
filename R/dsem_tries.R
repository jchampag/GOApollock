###############################################################
#script : tries to perform some dsem analyses / GOA pollock assessment

# Juliette Champagnat
# based on Cole's script dsem_esp.R
# started 4/1/24

###############################################################

rm(list=ls())
## devtools::install_github("James-Thorson-NOAA/dsem")
## install.packages("dsem")
library(TMB)
library(dsem)
library(ggplot2)
library(dplyr)
## devtools::install_github("afsc-assessments/GOApollock", ref='dev')
library(GOApollock)
library(readxl)
theme_set(theme_bw())

##### DATA READING + FORMATTING #######

# ESP data
dat0 <- read_xlsx('data/2023/goa_pollock_2023_dsem.xlsx', sheet=1, skip=3, na='NA')

# additional information about their lag
lag <- read_xlsx('data/2023/goa_pollock_2023_dsem.xlsx', sheet=1, range="B2:Z2",
                 na='NA', col_names=FALSE) %>% as.numeric
inc <- which(!is.na(lag))

## only tine series used in the BAS analysis
dat <- dat0[,inc]

## plot ESP data
g1 <- pivot_longer(dat0, -Year) %>% filter(!is.na(value)) %>%
  ggplot(aes(Year, value)) + geom_point() + facet_wrap('name', scales='free_y', ncol=4)
g1 # ggsave("esp_time_series_all.png", g1, width=15, height=8)
g2 <- pivot_longer(dat, -Year) %>% filter(!is.na(value)) %>%
  ggplot(aes(Year, value)) + geom_point() + facet_wrap('name', scales='free_y', ncol=4)
g2 # ggsave("esp_time_series_BAS.png", g2, width=15, height=6)


#### FIT POLLOCK USUAL ASSESSMENT MODEL #####

input_assess_usual <- prepare_pk_input(path='source', modfile='goa_pk_tmb',
                          datfile='../data/2023/pk23_10.txt', version='sigmaR')

input_assess_usual$map$sigmaR <- factor(1) # estimate
input_assess_usual$random <- 'dev_log_recruit' # recdevs

# dyn.load(dynlib('source/goa_pk_tmb'))
fit_assess_usual <- fit_pk(input=input_assess_usual, getsd=TRUE, newtonsteps=1, do.fit=1)
# fit_assess_usual$opt$objective

# plot recruits
fit_assess_usual$sd %>% filter(name=='log_recruit') %>%
  ggplot (aes(year, est, ymin=lwr, ymax=upr)) +
  geom_pointrange() + labs(y='log Recruits (billions)') +
  theme_bw()


#### FIT EXTERNAL DESM ####

#reduce ESP dataset to 3 variables
dat_names <-- dat %>% names();dat_names
dat_ext_dsem <- dat %>%
  select(2:4)# %>% ts()

#create the dsem formulation
sem = "
  # link, lag, param_name, start_value
  Annual_Heatwave_GOA_Model -> Annual_Heatwave_GOA_Model, 1, AR1, 1
  Spring_Chlorophylla_Peak_WCGOA_Satellite -> Spring_Chlorophylla_Peak_WCGOA_Satellite, 1, AR2, 1
  Annual_Auklet_Reproductive_Success_Chowiet_Survey -> Annual_Auklet_Reproductive_Success_Chowiet_Survey, 1, AR3, 1
#  Summer_Pollock_CPUE_YOY_Nearshore_Kodiak_Survey -> Summer_Pollock_CPUE_YOY_Nearshore_Kodiak_Survey, 1, AR4, 1
#  Summer_Pollock_Euphausiid_Diet_Juvenile_GOA_Survey -> Summer_Pollock_Euphausiid_Diet_Juvenile_GOA_Survey, 1, AR5, 1
#  Fall_Pollock_Condition_Adult_GOA_Fishery -> Fall_Pollock_Condition_Adult_GOA_Fishery, 1, AR6, 1
#  Summer_Pollock_Area_Occupied_WCGOA_Model -> Summer_Pollock_Area_Occupied_WCGOA_Model, 1, AR7, 1
#  Annual_Arrowtooth_Biomass_GOA_Model -> Annual_Arrowtooth_Biomass_GOA_Model, 1, AR8, 1
#  Annual_Pacific_Ocean_Perch_Biomass_GOA_Model -> Annual_Pacific_Ocean_Perch_Biomass_GOA_Model, 1, AR9, 1
#  Annual_Sablefish_Biomass_GOA_Model -> Annual_Sablefish_Biomass_GOA_Model, 1, AR10, 1
"

## prepare using package
family <- rep('fixed', ncol(dat_ext_dsem))
## family <- c('normal', 'normal', 'normal', 'fixed', 'fixed', 'fixed')
control <- dsem_control(use_REML=FALSE, run_model=TRUE,
                        getsd=TRUE, trace=100, newton_loops=0)
fit_ext_dsem <- dsem(sem=sem, tsdata=ts(dat_ext_dsem), family=family,
            estimate_delta0=FALSE, control=control)
dsem:::summary.dsem(fit_ext_dsem)
ParHat = fit$obj$env$parList()

## Quick plots, adapted from vignette cold pool example
oldpar <- par(no.readonly = TRUE)
par( mfcol=c(3,1), mar=c(2,2,2,0), mgp=c(2,0.5,0), tck=-0.02 )
for(i in 1:ncol(dat_ext_dsem)){
  tmp = dat_ext_dsem[,i,drop=FALSE]
  tmp = cbind( tmp, "PSEM"=ParHat$x_tj[,i] )
  SD = as.list(fit_ext_dsem$sdrep,what="Std.")$x_tj[,i]
  tmp = cbind( tmp, outer(tmp[,2],c(1,1)) +
                 outer(ifelse(is.na(SD),0,SD),c(-1,1)) )
  #
  plot( x=rownames(dat_ext_dsem), y=tmp[,1], ylim=range(tmp,na.rm=TRUE),
        type="p", main=colnames(dat_ext_dsem)[i], pch=20, cex=2 )
  lines( x=rownames(dat_ext_dsem), y=tmp[,2], type="l", lwd=2,
         col="blue", lty="solid" )
  polygon( x=c(rownames(dat_ext_dsem),rev(rownames(dat_ext_dsem))),
           y=c(tmp[,3],rev(tmp[,4])), col=rgb(0,0,1,0.2), border=NA )
}


### TRY FITING DSEM INSIDE ASSESSMENT MODEL ####

## pad extra years on beginning of time series
dat_int_dsem <- bind_rows(data.frame(Year=1970:1976), dat[,-1]) %>%
  select(2:4) %>%
  mutate(recdevs=NA) %>% relocate(recdevs)
# quick check that we are matching assessment duration #
dat_int_dsem %>% dim() # nyrs <- 54;stopifnot(length(dat_int_dsem$Year)==nyrs)

# plotting tentative
pairs(dat_int_dsem[,-1], upper.panel=NULL)

## prepare data using dsem package
sem = "
  # link, lag, param_name, start_value
  recdevs <-> recdevs, 0, sigmaR, 1
  Annual_Heatwave_GOA_Model -> Annual_Heatwave_GOA_Model, 1, AR1, 1
  Annual_Heatwave_GOA_Model -> recdevs, 0, beta, 1
  Spring_Chlorophylla_Peak_WCGOA_Satellite -> Spring_Chlorophylla_Peak_WCGOA_Satellite, 1, AR2, 1
  Annual_Auklet_Reproductive_Success_Chowiet_Survey -> Annual_Auklet_Reproductive_Success_Chowiet_Survey, 1, AR3, 1
#  Summer_Pollock_CPUE_YOY_Nearshore_Kodiak_Survey -> Summer_Pollock_CPUE_YOY_Nearshore_Kodiak_Survey, 1, AR4, 1
#  Summer_Pollock_Euphausiid_Diet_Juvenile_GOA_Survey -> Summer_Pollock_Euphausiid_Diet_Juvenile_GOA_Survey, 1, AR5, 1
#  Fall_Pollock_Condition_Adult_GOA_Fishery -> Fall_Pollock_Condition_Adult_GOA_Fishery, 1, AR6, 1
#  Summer_Pollock_Area_Occupied_WCGOA_Model -> Summer_Pollock_Area_Occupied_WCGOA_Model, 1, AR7, 1
#  Annual_Arrowtooth_Biomass_GOA_Model -> Annual_Arrowtooth_Biomass_GOA_Model, 1, AR8, 1
#  Annual_Pacific_Ocean_Perch_Biomass_GOA_Model -> Annual_Pacific_Ocean_Perch_Biomass_GOA_Model, 1, AR9, 1
#  Annual_Sablefish_Biomass_GOA_Model -> Annual_Sablefish_Biomass_GOA_Model, 1, AR10, 1
"
family <- rep('fixed', ncol(dat_int_dsem))
control <- dsem_control(use_REML=FALSE, run_model=FALSE)
fit_null_dsem = dsem(sem=sem, tsdata=ts(dat_int_dsem), family=family, control=control)
dsem:::summary.dsem(fit_null_dsem)

## try to fit it inside the model, recdevs are the first column of dat_int_dsem and x_tj
input_assess <- prepare_pk_input(path="source", modfile='goa_pk_dsem',
                          datfile='../data/2023/pk23_10.txt', version='dsem')


## stack on the dsem inputs to the assessment ones
input_int_dsem <- input_assess
input_int_dsem$dat <- c(input_assess$dat, fit_null_dsem$tmb_inputs$dat)
summary(input_int_dsem$dat)
# input_int_dsem$dat$y_tj[,1] <- NA # so not counted in NLL? already done i think
input_int_dsem$pars <- c(input_assess$pars, fit_null_dsem$tmb_inputs$parameters)
input_int_dsem$map <- c(input_assess$map, fit_null_dsem$tmb_inputs$map)
input_int_dsem$map$mu_j <- factor(c(NA,2:ncol(dat_int_dsem))) ## ??? ### map off recdev mean since already a parameter
input_int_dsem$pars$mu_j[1] <- 0
input_int_dsem$map$sigmaR <- NULL ## took out of model
input_int_dsem$random <- c(input_assess$random, fit_null_dsem$tmb_inputs$random) #?? ici on n'a que x_tj de random , etonne que ds l'assessment y en ait pas
input_int_dsem %>% summary()


#fit
# before fitting compile should have been done and ddl uploaded
# compile('source/goa_pk_dsem.cpp')
# dyn.load(dynlib('source/goa_pk_dsem'))
fit_assess_dsem <- fit_pk(input=input_int_dsem, getsd=1, newtonsteps=1, do.fit=1)
# xx <- sdreport(fit2$obj)
fit_assess_dsem$opt$objective

#plot recruitment
fit_assess_dsem$sd %>% filter(name=='log_recruit') %>%
  ggplot (aes(year, est, ymin=lwr, ymax=upr)) +
  geom_pointrange() + labs(y='log Recruits (billions)') +
  theme_bw()

filter(fit_assess_dsem$sd, name %in% c('mean_log_recruit', 'beta_z', 'mu_j'))


#### COMPARE RESCRUITMENT ESTIMATES ####

# estimated time series of recruitment
fit_assess_dsem$sd %>% mutate(version='dsem') %>%
  bind_rows(fit_usual_assess$sd) %>% filter(name=='log_recruit') %>%
  mutate(year=ifelse(version=='dsem', year-.2, year+.2)) %>% #pull(version) %>% unique()
  ggplot(aes(year, est, ymin=lwr, ymax=upr, color=version)) +
  geom_pointrange() + labs(y='log Recruits (billions)') +
  theme_bw()
fit_assess_dsem$opt$objective;fit_usual_assess$opt$objective


