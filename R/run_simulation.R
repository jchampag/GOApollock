##########################
# Script: to run some simulation (self-test) with //

# Juliette started 4/3/24
###########################

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
source("AssessDsem_fns.R")

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

start_time <- Sys.time()

install.packages("doSNOW")
install.packages("doParallel")
install.packages("doMPI")

library(doParallel)
library(foreach)
registerDoParallel(cl <- makeCluster(2))
results_list <- foreach(i = 1:2) %dopar% {

  simss <- simulate_AssessDsem(sem=sem_MS_MapMod,fit_assess_dsem =fit_assess_dsem_MS_MapMod,
                               family= family <- rep('normal',ncol(ESPdata_DAG1_rec)),
                               nsims = 1,fit_sim=FALSE,simverbose = FALSE,seed=i,
                               dsem_sim_control=list(resimulate_gmrf=TRUE,variance='none',ignoreNAs=FALSE))

  simss

}
stopCluster(cl)


end_time <- Sys.time()
print(end_time - start_time)


#-----------------------------------#
##### POST SIMULATON ANALYSIS #######
#-----------------------------------#


### Check cv

purrr::map_dfr(simu10,function(x){list(sim=x$version,max_grad=x$opt$max_gradient)})
get_std(simu10) %>% filter(is.na(est),is.na(se))

sim30 <- list(simu20,simu10)
simu30 <- unlist(sim30,recursive = FALSE)
simu10[[2]]$version

for(i in 1:10){
  simu10[[i]]$version <- 20+i
  simu10[[i]]$sd$version <- (20+i)
}



### Plot results

#### Parameters (FE)


sem_full<- fit_assess_dsem_DAG8_simple_2SSTdarrow_sd0.1$sem_full %>% as.data.frame()
get_std(fits=simu30) %>%
  filter(name %in%c('mean_log_recruit','beta_z')) %>%
  mutate(type='re-estimated',
         name2=rep(c('mean_log_recruit',sem_full$name),length(simu30))) %>% select(name2,year,version,type,est) %>%

  #add true (initially estimated) values
  bind_rows(data.frame(name2='mean_log_recruit',year=1970,
                       est=fit_assess_dsem_DAG8_simple_2SSTdarrow_sd0.1$parList$mean_log_recruit,
                       type='true')) %>%
  bind_rows(data.frame(name2=sem_full$name,
                       year=1970:(1970+length(fit_assess_dsem_DAG8_simple_2SSTdarrow_sd0.1$parList$beta_z)-1),
                       est=fit_assess_dsem_DAG8_simple_2SSTdarrow_sd0.1$parList$beta_z,
                       type='true')) %>%

  #compute RE
  group_by(year,name2) %>%
  mutate(RE=(est-est[which(type=='true')])/est[which(type=='true')]) %>%

  filter(type=='re-estimated') %>%
  # filter(!version%in%bad_grad) %>% #pull(version) %>% unique()
  #improve dsem par name
  mutate(par_type=ifelse(name2=='mean_log_recruit','assessement','dsem')) %>%
  filter(name2 %in% c('AR_JuvEuphDiet')) %>%

  ggplot(aes(x=name2,y=RE,group=name2,fill=par_type))+
  # geom_point(aes(x=name2,y=RE,group=name2,col=version))+
  geom_boxplot()+
  geom_hline(yintercept=0)+
  theme(legend.position = 'bottom',
        axis.text.x =element_text(angle=45))
# facet_wrap(~par_type,ncol=1)#,scale='free_y')

get_std(fits=simu20) %>%
  filter(name %in%c('mean_log_recruit','beta_z')) %>%
  mutate(type='re-estimated',
         name2=rep(c('mean_log_recruit',sem_full$name),length(simu20))) %>% select(name2,year,version,type,est) %>%
  filter(name2=='AR_JuvEuphDiet') %>% summary()

