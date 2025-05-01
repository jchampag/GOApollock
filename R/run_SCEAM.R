
# rm(list=ls())
## devtools::install_github("James-Thorson-NOAA/dsem")
## install.packages("dsem")
library(TMB)
library(dsem)
library(ggplot2)
library(dplyr)
## devtools::install_github("afsc-assessments/GOApollock", ref='dev')
# devtools::install_github("jchampag/GOApollock", ref='dev')
library(GOApollock)
library(readxl)
# theme_set(theme_bw())
library(phylopath)
library(ggpubr)
library(reshape)
library(gridExtra)
source("GitHub/GOApollock/R/AssessDsem_fns.R")


#----------------------------#
##### Get ESP data ready #####
#----------------------------#

## Load ESP data
ESPdata <- read.csv('GitHub/GOApollock/Data/2023/ESP_2024_4SCEAM.csv')

## set number of years for projection
ny_proj = 5

## get them ready
ESPdata_MS  <- ESPdata %>% 
  # adding NA before (match assessment begining) 
  bind_rows(data.frame(Year=c(1970:1976))) %>% arrange(Year) %>% 
  #  sdding NA after (forecast)
  add_row(Sum_OffYOY_Shelikof=rep(NA,ny_proj)) %>% 
  # selecting time series of interest
  select(Spr_SST,Wind_NS,Sum_LCopepod_Shelikof,Sum_Euph_Kodiak,
         Spr_Larvae_Shelikof,Sum_OffYOY_Shelikof,Sum_NearYOY_Kodiak,
         Sum_OffYOY_Cond_Shelikof,
         Sum_Juv_EuphDiet,Fal_Adult_Cond_Fishery) %>% 
  #scale and log data
  mutate(across(.cols=c(Spr_Larvae_Shelikof,Sum_OffYOY_Shelikof,Sum_NearYOY_Kodiak),.fns=log)) %>% 
  mutate(across(.fns =scale,.cols=-c(Spr_Larvae_Shelikof,Sum_OffYOY_Shelikof,Sum_NearYOY_Kodiak))) 
    
View(ESPdata_MS)  

#---------------------------------#
#### Run a SCEAM from scratch #####
#---------------------------------#  


sem_MapMod = "
  # link, lag, param_name, start_value

  #internal data relation
  Spr_SST -> Spr_SST,1,AR_SST
  Wind_NS -> Wind_NS,1,AR_Wind
  Sum_LCopepod_Shelikof -> Sum_LCopepod_Shelikof,1, AR_Cope
  Sum_Euph_Kodiak -> Sum_Euph_Kodiak,1,AR_Euph
  Spr_Larvae_Shelikof - > Spr_Larvae_Shelikof,1,AR_larvae
  Sum_OffYOY_Cond_Shelikof -> Sum_OffYOY_Cond_Shelikof,1,AR_CondOffYOY
  Sum_OffYOY_Shelikof -> Sum_OffYOY_Shelikof,1,AR_OffYOY
  Sum_NearYOY_Kodiak -> Sum_NearYOY_Kodiak,1,AR_NearYOY
  Sum_Juv_EuphDiet -> Sum_Juv_EuphDiet,1,AR_JuvEuphDiet
  Fal_Adult_Cond_Fishery -> Fal_Adult_Cond_Fishery,1,AR_CondAd
  recdevs <-> recdevs, 0, sigmaR, 1

  #causal relation
  Sum_Euph_Kodiak -> Fal_Adult_Cond_Fishery, 0, Euph_to_CondAd 
  Fal_Adult_Cond_Fishery -> Spr_Larvae_Shelikof,1,CondAd_to_Larvae
  #Spr_SST -> Spr_Larvae_Shelikof,0,SST_to_Larvae
  Wind_NS -> Spr_Larvae_Shelikof,0,Wind_to_Larvae
  Spr_Larvae_Shelikof -> Sum_OffYOY_Shelikof,0,Larvae_to_OffYOY
  Sum_OffYOY_Shelikof <-> Sum_NearYOY_Kodiak,0,OffYOY_to_NearYOY
  #Spr_SST -> Sum_OffYOY_Shelikof,0, SST_to_OffYOY
  Sum_OffYOY_Shelikof -> recdevs, 1, offYOY_to_R
  Spr_SST -> recdevs, 1,SST_to_R
  Spr_SST -> Sum_OffYOY_Cond_Shelikof,0,SST_to_CondYOY
  Sum_OffYOY_Cond_Shelikof -> recdevs,1,CondYOY_to_R
  "

#dataset
ESPdata_MSr <- ESPdata_MS %>% mutate(recdevs=NA) %>% relocate(recdevs)

# first non fit
family <- rep('normal',ncol(ESPdata_MSr))
DSEMcontrol <- dsem_control(use_REML=F, run_model=F,quiet=TRUE, getJointPrecision = TRUE,newton_loops=0)
fit_dsem = dsem(sem=sem_MapMod, tsdata=ts(ESPdata_MSr), family=family, control=DSEMcontrol )

#create new map & pars
new_map = fit_dsem$tmb_inputs$map
new_map$lnsigma_j <- factor(rep(NA, length=length(new_map$lnsigma_j)))
new_map <- fit_dsem$tmb_inputs$parameters
new_map$lnsigma_j <- rep(log(0.1), length=length(new_map$lnsigma_j))

input_assess <- prepare_pk_input(path="GitHub/GOApollock/data/2023/", modfile='goa_pk_dsem',
                                 datfile='pk23_10.txt', version='MapMod')#

input_assess$path <- 'GitHub/GOApollock/source'

## stack on the dsem inputs to the assessment ones
input_SCEAM <- input_assess

#dat
input_SCEAM$dat <- c(input_SCEAM$dat, fit_dsem$tmb_inputs$dat)
input_SCEAM$dat$y_tj[,1] <- NA
input_SCEAM$dat$Ftarget <- rep(input_SCEAM$dat$Ftarget[1],ny_proj) # set duration of projections
input_SCEAM$dat$indxsurv_log_sd4 <- input_SCEAM$dat$indxsurv_log_sd4*5 #decrease weight of some indices
input_SCEAM$dat$indxsurv_log_sd5 <- input_SCEAM$dat$indxsurv_log_sd5*5 #decrease weight of some indices

#pars
input_SCEAM$pars <- c(input_assess$pars, fit_dsem$tmb_inputs$parameters)
input_SCEAM$pars$lnsigma_j <- rep(log(0.1), length=length(fit_dsem$tmb_inputs$parameters$lnsigma_j)) #fixed value

#map
input_SCEAM$map <- c(input_assess$map, fit_dsem$tmb_inputs$map)
input_SCEAM$map$mu_j <- factor(c(NA,2:ncol(ESPdata_MSr))) # map off recdev mean since already a parameter
input_SCEAM$map$sigmaR <- NULL ## took out of model
input_SCEAM$map$lnsigma_j <- factor(rep(NA, length=length(fit_dsem$tmb_inputs$map$lnsigma_j))) #no estimated

#random
input_SCEAM$random <- c(input_assess$random, fit_dsem$tmb_inputs$random)

#INTERNAL MODEL RUN

# dyn.unload(dynlib('source/goa_pk_dsem'))
# compile('../source/goa_pk_dsem.cpp')
# dyn.load(dynlib('../source/goa_pk_dsem'))


fit_SCEAM <- fit_pk(input=input_SCEAM, getsd=TRUE, newtonsteps=1,
        do.fit=TRUE,save.sdrep = TRUE,filename = NULL)

fit_SCEAM$sem_full <-  fit_dsem$sem_full
class(fit_SCEAM) <- c(class(fit_SCEAM),'assessdsem')
fit_SCEAM


#--------------------------------------#
#### Run a SCEAM w function wraper #####
#--------------------------------------#  

sem_MapMod = "
  # link, lag, param_name, start_value

  #internal data relation
  Spr_SST -> Spr_SST,1,AR_SST
  Wind_NS -> Wind_NS,1,AR_Wind
  Sum_LCopepod_Shelikof -> Sum_LCopepod_Shelikof,1, AR_Cope
  Sum_Euph_Kodiak -> Sum_Euph_Kodiak,1,AR_Euph
  Spr_Larvae_Shelikof - > Spr_Larvae_Shelikof,1,AR_larvae
  Sum_OffYOY_Cond_Shelikof -> Sum_OffYOY_Cond_Shelikof,1,AR_CondOffYOY
  Sum_OffYOY_Shelikof -> Sum_OffYOY_Shelikof,1,AR_OffYOY
  Sum_NearYOY_Kodiak -> Sum_NearYOY_Kodiak,1,AR_NearYOY
  Sum_Juv_EuphDiet -> Sum_Juv_EuphDiet,1,AR_JuvEuphDiet
  Fal_Adult_Cond_Fishery -> Fal_Adult_Cond_Fishery,1,AR_CondAd
  recdevs <-> recdevs, 0, sigmaR, 1

  #causal relation
  Sum_Euph_Kodiak -> Fal_Adult_Cond_Fishery, 0, Euph_to_CondAd 
  Fal_Adult_Cond_Fishery -> Spr_Larvae_Shelikof,1,CondAd_to_Larvae
  #Spr_SST -> Spr_Larvae_Shelikof,0,SST_to_Larvae
  Wind_NS -> Spr_Larvae_Shelikof,0,Wind_to_Larvae
  Spr_Larvae_Shelikof -> Sum_OffYOY_Shelikof,0,Larvae_to_OffYOY
  Sum_OffYOY_Shelikof <-> Sum_NearYOY_Kodiak,0,OffYOY_to_NearYOY
  #Spr_SST -> Sum_OffYOY_Shelikof,0, SST_to_OffYOY
  Sum_OffYOY_Shelikof -> recdevs, 1, offYOY_to_R
  Spr_SST -> recdevs, 1,SST_to_R
  Spr_SST -> Sum_OffYOY_Cond_Shelikof,0,SST_to_CondYOY
  Sum_OffYOY_Cond_Shelikof -> recdevs,1,CondYOY_to_R
  "

fit_SCEAM2 <- AssessDsem_fit(fit_type = "AssessDsem",
                                            fit_name = '_MS_MapMod' ,
                                            ESPdata=ESPdata_MS,
                                            sem=sem_MapMod,
                                            sem_family =  rep('normal', ncol(ESPdata_MS)+1),
                                            family_sd = rep(0.1, ncol(ESPdata_MS)+1),
                                            ny_proj=ny_proj,
                                            Assess_recdevs= NULL,
                                            save_fit = FALSE,
                                            save.sdrep=FALSE,
                                            src_path = 'GitHub/GOApollock/source',
                                            data_path = 'GitHub/GOApollock/data/2023/')


#---------------------------------------#
#### Run predicitive skill-testing ######
#---------------------------------------#

predskill_MapMod <- retro_proj_analysis_AssessDsem(sem=sem_MapMod,
                                                        fit=fit_SCEAM2,
                                                        peels=0:10,env_data = 'real',ny_proj = 5)


retrAD_real_MS_MapMod %>% length()
skillpred_real_MS_MapMod <- unlist(retrAD_real_MS_MapMod,recursive = FALSE)
length(skillpred_real_MS_MapMod)
save(skillpred_real_MS_MapMod,file='skilltest_real_MS_MapMod.Rdata')


#------------------------#
#### Run self-testing ####
#------------------------#


