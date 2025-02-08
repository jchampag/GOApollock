##########################
# Script: to run some simulation (self-test) with //

# Juliette started 4/3/24
###########################
setwd('C:/Users/jchampag/Documents/Github/GOApollock')

library(TMB)
## devtools::install_github("James-Thorson-NOAA/dsem", ref='dev')
## devtools::install_github("Cole-Monnahan-NOAA/dsem", ref='simulate_NAs')
library(dsem)
library(ggplot2)
library(dplyr)
## devtools::install_github("afsc-assessments/GOApollock", ref='dev')
# devtools::install_github("jchampag/GOApollock", ref='dsem')
library(GOApollock)
library(readxl)
# theme_set(theme_bw())
# library(phylopath)
# library(ggpubr)
# library(reshape)
library(gridExtra)
library(purrr)
source("R/AssessDsem_fns.R")

# compile('source/goa_pk_dsem.cpp')
# dyn.load(dynlib('source/goa_pk_dsem'))
#--------------------------------------------------------------------------------#
###### WITH PARALLEL AND PARLAPPLY - not working ####
#--------------------------------------------------------------------------------#

# library(parallel)
# cl <- makeCluster(60)
#
# clusterEvalQ(cl, {library(tidyverse);library(TMB);library(dsem);library(GOApollock)})
#
# clusterExport(cl,
#               c("simulate_AssessDsem", #functions
#                 'sem_DAG1','fit_assess_dsem_DAG1_sd0.1','ESPdata_DAG1_rec'#objects
#               ))
#
# results <- parLapply(cl,seq(1,nrow(simu_par_all)),fun=MSYanalysis_wraper,par=simu_par_all,F_seq = seq(0,1.5,0.01),
#                      years=c(1:1000),
#                      data=list_params)
# results <- parLapply(cl,seq(1,nrow(simu_par_all)),fun=simulate_AssessDsem,sem=sem_DAG1,fit_assess_dsem =fit_assess_dsem_DAG1_sd0.1,
#                      family= family <- rep('normal',ncol(ESPdata_DAG1_rec)),nsims = 1)
# stopCluster(cl)
#
# df_results <- data.table::rbindlist(results)
#
# fit_assess_dsem = fit_assess_dsem_DAG1_sd0.1
# fit_assess_dsem_DAG1_sd0.1$path <-'source'
# fit_assess_dsem$path <-'source'
# fit_assess_dsem_DAG1_sd0.1$input$path <-'source'
# # input_int_dsem4simu$path <-'source'
# sem=sem_DAG1
# #
# trysimu2 <- simulate_AssessDsem(sem=sem_DAG1,fit_assess_dsem =fit_assess_dsem_DAG1_sd0.1,
#                     family= family <- rep('normal',ncol(ESPdata_DAG1_rec)),nsims = 1)
#
# get_std(fits=re_fit_assess_dsem) %>%
#   filter(name %in%c('Espawnbio','mean_log_recruit','recruit','log_recruit','beta_z'))


#--------------------------------------------------------------------------------#
###### WRAPED IN A FUNCTION ####
#--------------------------------------------------------------------------------#

# library(snowfall)
# repretros <- do_skill_test(reps=1:2,sem=sem_DAG1,fit=fit,
#                            peels=0:1,ny_proj=3,parallel=TRUE)
#
#
#
# 'do_skill_test' <- function(reps=1:20,sem,fit,peels=0:10,ny_proj=3,env_data='idealistic',parallel=TRUE){
#   if(class(fit)[1]!='pkfit')
#     stop("fit argument is not a fitted model")
#   if(parallel){
#     message("Preparing parallel session..")
#     if(!require(snowfall))
#       stop("snowfall package required for parallel execution")
#     sfInit(parallel=TRUE, cpus=parallel::detectCores()-1)
#     sfLibrary(GOApollock);sfLibrary(dsem)
#     sfSource(file='R/AssessDsem_fns.R')
#     sfSource(file='R/retro_fns.R')
#     # sfExport(list('env_data'=env_data),local = TRUE)
#     rep_retros <- sfLapply(reps, function(i) retro_proj_analysis_AssessDsem(reps=i,
#                                                                             sem=sem,fit=fit,peels=peels,ny_proj = ny_proj))
#     sfStop()
#   } else {
#     rep_retros <- lapply(reps, function(i) retro_proj_analysis_AssessDsem(reps=i,sem=sem,fit=fit,peels=peels,ny_proj = ny_proj,env_data = env_data))
#   }
#   return(rep_retros)
# }
#
#
#
# library(parallel)
# cl <- parallel::makeCluster(2)#parallel::detectCores()-5)
#
# parallel::clusterEvalQ(cl, {library(tidyverse);library(GOApollock);library(dsem)})
#
# parallel::clusterExport(cl,
#                         c(#functions
#                           'retro_proj_analysis_AssessDsem','fit_retro_proj_AssessDsem','simulate_AssessDsem',
#                           'peel_pars','peel_map','peel_data',
#                           #object
#                           'fit','reps','sem','peels','ny_proj'))
#
# rep_retros <- parallel::parLapply(cl,1:2, fun= retro_proj_analysis_AssessDsem,
#                                   sem=sem,fit=fit,peels=0:1,ny_proj = 3)
# # results <- parLapply(cl,seq(1,nrow(simu_par_all)),fun=MSYanalysis_wraper,par=simu_par_all,F_seq = seq(0,1.5,0.01),
# #                      years=c(1:1000),
# #                      data=list_params)
# stopCluster(cl)
#
# df_results <- data.table::rbindlist(results)

#--------------------------------------------------------#
######## WITH FOREACH ########
#--------------------------------------------------------#



# install.packages("doSNOW")
# install.packages("doParallel")
# install.packages("doMPI")

library(doParallel)
library(foreach)

load('runs/fit_MS_MapMod.Rdata')
fit <- fit_assess_dsem_MS_MapMod
fit$path <- fit$input$path <-  "source"
parallel::detectCores()

ncores= 30

start_time <- Sys.time()
registerDoParallel(cl <- makeCluster(ncores))
results_list <- foreach(i = 1:ncores,.errorhandling = 'pass') %dopar% {

  library(dsem)
  library(GOApollock)
  source("R/AssessDsem_fns.R")
  simss <- simulate_AssessDsem(sem=fit$sem,fit_assess_dsem =fit,
                               family= rep('normal',ncol(fit$input$dat$y_tj)),
                               nsims = 20,fit_sim=TRUE,simverbose = FALSE,seed=30+i,
                               dsem_sim_control=list(resimulate_gmrf=TRUE,variance='none',ignoreNAs=FALSE))

  simss

}

results_list_noNA <- foreach(i = 1:ncores,.errorhandling = 'pass') %dopar% {
  
  library(dsem)
  library(GOApollock)
  source("R/AssessDsem_fns.R")
  simss <- simulate_AssessDsem(sem=fit$sem,fit_assess_dsem =fit,
                               family= rep('normal',ncol(fit$input$dat$y_tj)),
                               nsims = 20,fit_sim=TRUE,simverbose = FALSE,seed=30+i,
                               dsem_sim_control=list(resimulate_gmrf=TRUE,variance='none',ignoreNAs=TRUE))
  
  simss
  
}

stopCluster(cl)
end_time <- Sys.time()
print(end_time - start_time)


simu

save(simu,file='runs/selftest_MapMod_noNAs.Rdata')
# results_list_NA <- results_list
#-----------------------------------#
##### POST SIMULATON ANALYSIS #######
#-----------------------------------#
load('runs/selftest_MapMod.Rdata')
onelist <- unlist(results_list_NA,recursive = F)
twolist <- unlist(results_list,recursive = F)
onetwolist <- list(onelist,twolist)
simu <- unlist(results_list_noNA,recursive = F) #results_list_noNA
length(simu)
which(names(simu)!='')
simu <- discard_at(simu,which(names(simu)!=''))#c(56,57)
length(simu)

### Check cv

purrr::map_dfr(simu,function(x){data.frame(sim=x$version,max_grad=x$opt$max_gradient)})
purrr::map_dfr(simu,function(x){data.frame(sim=x$version,max_grad=x$opt$max_gradient)}) %>% filter(max_grad>0.1) %>% nrow()
bad_grad <- purrr::map_dfr(simu,function(x){list(sim=x$version,max_grad=x$opt$max_gradient)}) %>% filter(max_grad>0.1) %>% pull(sim)
get_std(simu) %>% filter(is.na(est))
bad_se <- get_std(simu) %>% filter(is.na(se)) %>% pull(version) %>% as.character() %>% unique()

length(bad_grad);length(bad_se)
match(bad_se,bad_grad);match(bad_grad,bad_se)

# sim_names <- purrr::map_vec(simu,function(x){x$version})
# sim_names[!sim_names %in%bad_se] %>% length()
N = length(simu)-length(bad_se);N

## RE for Parameters (FE)

parplot <- purrr::map_dfr(simu,function(x){as.data.frame(fit$sem_full) %>% 
    mutate(est=x$parList$beta_z,true=fit$parList$beta_z,sim=x$version) %>% #type='selftest',
    select(c(name,lag,first,second,est,true,sim)) %>% 
    bind_rows(data.frame(name='mean_log_recruit',#lag=NA,first=NA,second=
                         est=x$parList$mean_log_recruit,true=fit$parList$mean_log_recruit,
                         sim=x$version))}) %>%
  filter(est!=0) %>%
  filter(!sim%in%bad_se) %>% 
  mutate(RE=(est-true)/true) %>% 
  ggplot(aes(x=name,y=RE,group=name))+#fill=par_type
  # geom_point(aes(x=name2,y=RE,group=name2,col=version))+
  geom_boxplot()+
  geom_hline(yintercept=0)+
  theme(legend.position = 'bottom',
        axis.text.x =element_text(angle=45))+ #coord_flip()+
  labs(y='Relative error',x="Parameter names")
    #title=paste0('N=',N,' resimulatedGMRF=',simu[[1]]$DSEM_simu_settings$resimulate_gmrf,
    #                ' variance=',simu[[1]]$DSEM_simu_settings$variance,' fillNA=',simu[[1]]$DSEM_simu_settings$ignoreNAs))

### RE for quantities

qplot <- purrr::map_dfr(simu,function(x){data.frame(what='Recruitment',year=x$rep$years,est=x$rep$recruit,true=x$simu_int_dsem$recruit) %>% 
    bind_rows(data.frame(what='SSB',year=x$rep$years,est=x$rep$Espawnbio,true=x$simu_int_dsem$Espawnbio)) %>% 
    bind_rows(data.frame(what='Catch',year=x$rep$years,est=x$rep$Ecattot,true=x$simu_int_dsem$Ecattot)) %>% 
    mutate(sim=x$version)}) %>%
  # filter(est!=0) %>% 
  filter(!sim%in%bad_se) %>% 
  mutate(RE=(est-true)/true) %>% 
  ggplot(aes(x=year,y=RE,group=year))+facet_wrap(~what,ncol=1,scale='free_y')+
  # geom_point(aes(x=name2,y=RE,group=name2,col=version))+
  geom_boxplot()+#outlier.shape = NA
  # scale_y_continuous(limits=c(0,30))+#quantile(RE, c(0.1, 0.9)))+
  # gg.layers::geom_boxplot2(width.errorbar = 0.5)+
  geom_hline(yintercept=0)+
  # theme(legend.position = 'bottom',
  #       axis.text.x =element_text(angle=45))+
  labs(y='Relative error',x="Year")
# title=paste0('N=',N,' resimulatedGMRF=',simu[[1]]$DSEM_simu_settings$resimulate_gmrf,
#                     ' variance=',simu[[1]]$DSEM_simu_settings$variance,' fillNA=',simu[[1]]$DSEM_simu_settings$ignoreNAs))


## Compose a plot w parameter and quatites RE
ggpubr::ggarrange(plotlist = list(parplot,qplot), nrow=2,ncol = 1)

## RE for Rec

purrr::map_dfr(simu,function(x){data.frame(what='Recruitment',year=x$rep$years,est=x$rep$recruit,true=x$simu_int_dsem$recruit)%>% #fit$rep$recruit,
    mutate(sim=x$version)}) %>%
  # filter(est!=0) %>% 
  filter(!sim%in%bad_se) %>% 
  mutate(RE=(est-true)/true) %>% 
  ggplot(aes(x=year,y=RE,group=year))+
  geom_boxplot()+ #outlier.shape = NA
  # scale_y_c#ontinuous(limits=c(0,40))+
  geom_hline(yintercept=0)+
  theme(legend.position = 'bottom',
        axis.text.x =element_text(angle=45))+
  labs(title=paste0('Recruitment - N=',N,' resimulatedGMRF=',simu[[1]]$DSEM_simu_settings$resimulate_gmrf,
                    ' variance=',simu[[1]]$DSEM_simu_settings$variance,' fillNA=',simu[[1]]$DSEM_simu_settings$ignoreNAs))


## RE for Catch

purrr::map_dfr(simu,function(x){data.frame(what='Catch',year=x$rep$years,est=x$rep$Ecattot,true=x$simu_int_dsem$Ecattot)%>% 
    mutate(sim=x$version)}) %>%
  # filter(est!=0) %>% 
  filter(!sim%in%bad_se) %>% 
  mutate(RE=(est-true)/true) %>% 
  ggplot(aes(x=year,y=RE,group=year))+
  geom_boxplot()+ #outlier.shape = NA
  # scale_y_continuous(limits=c(0,20))+
  geom_hline(yintercept=0)+
  theme(legend.position = 'bottom',
        axis.text.x =element_text(angle=45))+
  labs(title=paste0('Catch - N=',N,' resimulatedGMRF=',simu[[1]]$DSEM_simu_settings$resimulate_gmrf,
                    ' variance=',simu[[1]]$DSEM_simu_settings$variance,' fillNA=',simu[[1]]$DSEM_simu_settings$ignoreNAs))


## RE for SSB

purrr::map_dfr(simu,function(x){data.frame(what='SSB',year=x$rep$years,est=x$rep$Espawnbio,true=x$simu_int_dsem$Espawnbio)%>% 
    mutate(sim=x$version)}) %>%
  # filter(est!=0) %>% 
  filter(!sim%in%bad_se) %>% 
  mutate(RE=(est-true)/true) %>% 
  ggplot(aes(x=year,y=RE,group=year))+
  geom_boxplot()+ #outlier.shape = NA
  # scale_y_continuous(limits=c(0,5))+
  geom_hline(yintercept=0)+
  theme(legend.position = 'bottom',
        axis.text.x =element_text(angle=45))+
  labs(title=paste0('SSB - N=',N,' resimulatedGMRF=',simu[[1]]$DSEM_simu_settings$resimulate_gmrf,
                    ' variance=',simu[[1]]$DSEM_simu_settings$variance,' fillNA=',simu[[1]]$DSEM_simu_settings$ignoreNAs))
