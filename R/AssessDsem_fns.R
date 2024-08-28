##########################
# Script: some functions modified from dsem/GOApollock for my use + new functions

# Juliette started 4/3/24
###########################

#Improvements to think of : none for now

#' Fit a model (assessment, dsem or assessment+dsem) to data
#' @param fit_type character chain indicating which model to fit
#'   "AssessOnly", "ExtDsem" or "AssessDsem"
#' @param fit_name name to be passed in \code{version} aurgument from \code{\link{fit_pk}}() and used if \code{save_fit=TRUE}
#' @param ESPdata Dataset with environmental date if a DSEM model is fit
#' @param sem sem statement from \code{\link{dsem}}() if a DSEM model is fit
#' @param sem_family family type for each variables of the sem model from \code{\link{dsem}}()
#' @param ny_proj number of year to project
#' @param Assess_recdevs if \code{fit_type = "ExtDSEM"} the needs a recruitment deviation series to fit to
#' @param save_fit TRUE/FALSE
#' @return A list object of class depending of  \code{fit_type}
#'  * 'pkfit' if \code{fit_type = "AssessOnly"}
#'  * 'dsem' if \code{fit_type = "ExtDSEM"}
#'  * 'pkfit' AND 'assessdsem' if \code{fit_type = "AssessDsem"}
#' @details This function is beta still and subject to change
#'   without warning.
#' @export

'AssessDsem_fit' <- function(fit_type = "AssessOnly", #or "ExtDsem" or "AssessDsem"
                             fit_name = '_1' ,
                             ESPdata,
                             sem=sem,
                             sem_family =  rep('fixed', ncol(ESPdata)),
                             ny_proj=15,
                             newtonsteps=1,
                             Assess_recdevs= NULL, #if ExtDsem it it requiered to have it
                             save_fit = TRUE,
                             Rdsem=TRUE
){
  #some warnings
  if(is.null(Assess_recdevs)& fit_type=="ExtDsem"){print('cannot fit without recdev estimates')}


  if(fit_type != 'ExtDsem'){ #if AssessOnly or AssessDsem need to do Assess stuff

    input_assess <- prepare_pk_input(path='../source',
                                     modfile=ifelse(fit_type=='AssessOnly','goa_pk_tmb','goa_pk_dsem'),
                                     datfile='../data/2023/pk23_10.txt', version=paste0(fit_type,fit_name))

    input_assess$dat$Ftarget <- rep(input_assess$dat$Ftarget[1],ny_proj) #/!

    # decreasing influence of spawner survey
    input_assess$dat$indxsurv_log_sd4 <- input_assess$dat$indxsurv_log_sd4*5
    input_assess$dat$indxsurv_log_sd5 <- input_assess$dat$indxsurv_log_sd5*5

    if(fit_type == 'AssessOnly'){
      input_assess$map$sigmaR <- factor(1) # estimate
      input_assess$random <- 'dev_log_recruit' #

      #fit
      fit_assess <- fit_pk(input=input_assess,
                           getsd=TRUE, newtonsteps=newtonsteps,
                           save.sdrep = TRUE,
                           filename = if(save_fit){paste0(fit_type,fit_name,'.RDS')}else{NULL})

      out <- fit_assess
    }
  }

  if(fit_type != 'AssessOnly'){ # if AssessOnly, do not need all the DSEM stuff

    #prepare dsem fit
    ESPdata <- ESPdata %>%
      mutate(recdevs=if(fit_type=="ExtDsem"){c(Assess_recdevs,rep(NA,ny_proj))}else{rep(NA,nrow(ESPdata))}) %>%
      relocate(recdevs)

    # should it be done externally because then I'm loosing the possibility to change the control object...?
    control_dsem <- dsem_control(use_REML=FALSE, run_model=ifelse(fit_type == 'ExtDsem',TRUE,FALSE),
                                 getsd=TRUE, trace=100, newton_loops=1)

    fit_dsem <- dsem(sem=sem, tsdata=ts(ESPdata), family=sem_family,
                     estimate_delta0=FALSE, control=control_dsem)
    fit_dsem$version <- paste0(fit_type,fit_name)

    out <- fit_dsem
  }

  if(fit_type == "AssessDsem"){

    #prepare data
    input_assessDsem <- input_assess
    input_assessDsem$dat <- c(input_assess$dat, fit_dsem$tmb_inputs$dat)
    input_assessDsem$dat$y_tj[,1] <- NA # so not counted in NLL? already done i think, is it usefull? #kind of a double check Iguess
    input_assessDsem$pars <- c(input_assess$pars, fit_dsem$tmb_inputs$parameters)
    input_assessDsem$map <- c(input_assess$map, fit_dsem$tmb_inputs$map)
    input_assessDsem$random <- c(input_assess$random, fit_dsem$tmb_inputs$random)

    #specific mapping
    input_assessDsem$map$mu_j <- factor(c(NA,2:ncol(ESPdata))) #
    input_assessDsem$pars$mu_j[1] <- 0 #already done by dsem, usefull?
    input_assessDsem$map$sigmaR <- NULL ## took out of model

    if(!Rdsem){
      ids <-  as.data.frame(fit_dsem$sem_full)
      id_Rlink <- which(ids$first!=ids$second & ids$second=='recdevs')
      input_assessDsem$map$beta_z <- 1:length(input_assessDsem$pars$beta_z )
      input_assessDsem$map$beta_z[id_Rlink] <- NA
      input_assessDsem$map$beta_z <- as.factor(input_assessDsem$map$beta_z)
      input_assessDsem$pars$beta_z[id_Rlink] <- 0
    }

    fit_assessDsem <- fit_pk(input=input_assessDsem,
                             getsd=TRUE, newtonsteps=newtonsteps, do.fit=1,
                             save.sdrep = TRUE,
                             filename = if(save_fit){paste0(fit_type,fit_name,'.RDS')}else{NULL})

    fit_assessDsem$sem_full <- fit_dsem$sem_full
    class(fit_assessDsem) <- c(class(fit_assessDsem),'assessdsem')
    out <- fit_assessDsem
  }

  return(out)

}

#--------------------------------------------------------------------------#

#Improvements to think of : none for now + made some adjustement for matrix output for plotting purpose, does it need to be a special case?

#' Return objects able to display estimated internal and external relation between variables
#' based on \code{\link{as_fitted_DAG}}()
#' @param fit fit, should be of class 'dsem' or 'assessdsem'
#' @param fit_null_dsem if the fit is not directly from \code{\link{dsem}}()
#' @param lag which lag to output
#' @param direction whether to include one-sided arrows direction=1, or both one- and two-sided arrows direction=c(1,2)
#' @param what whether to output estimates what="Estimate", standard errors what="Std_Error" or p-values what="Std_Error"
#' @param which_link with link to display, 'all' or causal only
#' @param out_type type of output to return (see return) else it will be
#' @return An object
#'  * a matrix ready to be use for a DAG plot if \code{out_type = "matrix"}
#'  * a table with all link value elsewise
#' @details This function is beta still and subject to change
#'   without warning.
#' @export

'as_fitted_DAG_assessdsem' <- function (fit,
                                     #fit_null_dsem=NULL,#/!\ this part could be simplified, to evaluate later
                                     lag = c(0,1), direction = c(1,2),
                                     what = "Estimate",which.link='causal',
                                     out.type="matrix")
{

  if(any(class(fit)=='dsem')){ coefs = summary(fit);  vars = colnames(fit$tmb_inputs$data$y_tj)
  }else{#cannot use summary(fit) for a 'assessdsem' output, this part replace what it does
    if(all((class(fit)!='assessdsem'))){print('This function can only work with a dsem like output')}

    model = fit$sem_full#fit_null_dsem$sem_full #/!\ this part could be simplified, to evaluate later
    ParHat <- list('beta_z'= fit$sd %>% filter(name== 'beta_z') %>% pull(est),
                   'mu_j'= fit$sd %>% filter(name== 'mu_j') %>% pull(est))#,

    coefs = data.frame( model, "Estimate"=c(NA,ParHat$beta_z)[ as.numeric(model[,'parameter'])+1 ] ) # parameter=0 outputs NA
    coefs$Estimate = ifelse( is.na(coefs$Estimate), as.numeric(model[,4]), coefs$Estimate )
    coefs = data.frame( coefs, "Std_Error"=fit$sd %>% filter(name== 'beta_z') %>% pull(se))
    coefs = data.frame( coefs, "z_value"=coefs[,'Estimate']/coefs[,'Std_Error'] )
    coefs = data.frame( coefs, "p_value"=pnorm(-abs(coefs[,'z_value'])) * 2 )
    vars = colnames(fit$input$dat$y_tj)

  }

  coefs = coefs[which(coefs[, 2] %in% lag), ]#this allows selecting more than 1 lag to represent
  coefs = coefs[which(coefs[, "direction"] %in% direction), ]

  if(out.type!='matrix'){
    if(which.link!= 'all'){coefs = coefs[which(coefs[, 6] !=coefs[, 7] ),-c(4,5) ];return(coefs)}
    return(coefs[,-c(4,5)])
  }else{
    coefs = coefs %>% filter(direction==1,first!=second)
    vars = c(coefs %>% pull(first),coefs %>% pull(second)) %>% unique()

    out = list(coef = array(0, dim = rep(length(vars), 2), dimnames = list(vars,
                                                                             vars)))
    out$coef[as.matrix(coefs[, c("first", "second")])] = coefs[,
                                                                 what]
    class(out) = "fitted_DAG"
    return(out)
  }
}
#
# as_fitted_DAG_assessdsem(fit=fits[[1]],#fit_assess_dsem_proj,#fit_ext_dsem_proj_dropidx,#fits[[1]],
#                          #fit_null_dsem = NULL,#fit_ext_dsem_proj,#NULL,#
#                                         #fit_null_dsem=NULL,#/!\ this part could be simplified, to evaluate later
#                                         lag = c(0,1), direction = c(1,2),
#                                         what = "p_value",which.link='all',
#                                         out.type="matrix")
#
#--------------------------------------------------------------------------#

#Improvements to think of : none for now

#' Return a plot comparing the log recruitment estimations for different type of fits
#' @param fits a list of fits to compare
#' @param ny_proj nb of year projected
#' @param mean_log_rec value to add to recruitment deviations from dsem fit
#' @param extDsem_name name to be passed out for the plot (only for ext dsem fits, pkfit uses version)
#' @param points only lines or also diaply estimate points
#' @return A plot
#' @export


'rec_plot' <- function(fits,
                       ny_proj=15,
                       mean_log_rec=0,
                       extDsem_name=NULL,
                       points= FALSE)


{
  tmp <- NULL
  log_mean_rec_pointer = 0

  for(i in seq_along(fits)){
    if(any(class(fits[[i]])== 'pkfit')){
      tmp <- tmp %>%
        bind_rows(fits[[i]]$sd %>%  filter(name %in% c('log_recruit_proj','log_recruit')) %>%
                                             mutate(year=ifelse(name=='log_recruit',year,year+54)))

    }else{

      if(class(fits[[i]])!= 'dsem'){print('This methods cannot be applied to this class of fit')}
      log_mean_rec_pointer = log_mean_rec_pointer+1
      tmp <- tmp %>%
        bind_rows(data_frame(name="", est=as.list(fits[[i]]$sdrep,what="Estimate")$x_tj[,1]+mean_log_rec[log_mean_rec_pointer],
                             se=as.list(fits[[i]]$sdrep,what="Std.")$x_tj[,1],
                             year=1970:(2023+ny_proj)) %>%
                    mutate(lwr=est-qnorm(0.975)*se,
                           upr=est+qnorm(0.975)*se,
                           version=fits[[i]]$version))
    #                          ifelse(is.null(extDsem_name),paste0("extDsem_",log_mean_rec_pointer),
    #                                       paste0("extDsem_",extDsem_name[log_mean_rec_pointer]))))
    }
  }

  plt <- tmp %>%
    ggplot(aes(x=year,y=est,ymin=upr,ymax=lwr))+
    geom_line(aes(col=version))+
    geom_ribbon(aes(fill=version),alpha=0.2)+
    labs(y="Estimated log recruitment")+
    theme(legend.position='bottom',legend.key.size = unit(1, 'cm'),legend.text = element_text(size=12),
          axis.title = element_text(size=12),axis.text = element_text(size=10))+
    scale_x_continuous(n.breaks = 10)+
    geom_vline(xintercept=2023)

  if(points){return(plt+geom_point(aes(col=version)))}else{return(plt)}
}

# fits=list(fit_assess_usual_proj,
#                              fit_ext_dsem_proj,
#                              fit_ext_dsem_proj_dropidx)



# rec_plot(fits=list(fit_assess_usual_proj,
#                    fit_ext_dsem_proj,
#                    fit_ext_dsem_proj_dropidx),
#          ny_proj=15,
#          mean_log_rec=c(0,1))#,
#          #extDsem_name = c('nodropidx','dropidx')     )

#--------------------------------------------------------------------#


#Improvements to think of : none for now

#' Return a plot comparing the estimations of sigmaR for different type of fits
#' @param fits a list of fits to compare
#' @return A plot
#' @export

'sigmaR_plot' <- function(fits){

  #fits_class <- purrr::map(fits,.f=class)
  tmp <- NULL

  for(i in seq_along(fits)){

    if(any(class(fits[[i]])== 'assessdsem')){
      #have to grab the nb of beta_z corresponding to sigmaR
      id <-  which(as_fitted_DAG_assessdsem(fits[[i]],out.type = 'summary',which.link = 'all')$name=='sigmaR')
      tmp <- tmp %>%
        bind_rows(fits[[i]]$sd %>%  filter(name=='beta_z') %>% filter(year==(1969+id)))

    }else{
      # is this one necessary? if yes have to think to get same str than fit$sd + give a version name
      # if(any(class(fits[[i]])== 'dsem')){
      #   tmp <- tmp %>%
      #     bind_rows(fits[[i]]$sd %>%  filter(name=='sigmaR'))
      # }else{
        tmp <- tmp %>%
          bind_rows(fits[[i]]$sd %>%  filter(name=='sigmaR'))
    }
  }

  plt <- tmp %>%
    ggplot(aes(x=version,ymin=lwr,y=est,ymax=upr))+geom_pointrange()+
    geom_hline(yintercept=0)+
    coord_flip()+labs(y="Sigma R estimates")+
    theme(axis.title = element_text(size=14),axis.text = element_text(size=12))


  return(plt)
}

# sigmaR_plot(fits=list(fit_assess_dsem_proj,fit_assess_dsem_proj_dropidx,
#                       fit_assess_usual_proj,fit_assess_usual_proj_dropidx))



#--------------------------------------------------------------#

#Improvements to think of :

#' Return a plot comparing the parameters estimated by dsem
#' @param fits a list of fits to compare
#' @param what which type of parameters are in the plot all/causal/variance/AR1
#' @return A plot
#' @export

'compute_var_red' <- function(new_fit,ref_fit){
  ref <- ref_fit$sd %>% filter(name=='sigmaR')%>% pull(est)
  new <- new_fit$sd %>% filter(name=='beta_z') %>% filter(year==1970) %>% pull(est)

  red = ((new^2)-(ref^2))/(ref^2)
  return(red)
}

#----------------------------------------------------------------------#


#Improvements to think of : add subtitle causal-link/p-value , have an automatic label related to the fit version
# improve the DAG plot in general for ex:
#   * changing color/size of arrow if the pvalue is low
#   * having the pvalue plotted on the same graph
#   * add a marker for different temporal lag

#' Return a plot comparing the estimated causal links
#' @param fits a list of fits to compare
#' @return A plot
#' @export


'causal_plot' <- function(fits,
                          fit_label = c('fit1','fit2')){

  plt_link <- plt_pval <- list()
  # tour =1
  for(i in seq_along(fits)){

    plt_link[[i]] <- plot( (as_fitted_DAG_assessdsem(fit=fits[[i]])),
                           edge.width=1, type="width",
                           text_size=4, show.legend=FALSE,
                           arrow = grid::arrow(type='closed', 18, grid::unit(10,'points')) ) +
      scale_x_continuous(expand = c(0.4, 0.1))


    # tour=tour+2
  }
  for(i in seq_along(fits)){
    plt_link[[length(fits)+i]]  <-  plot(  (as_fitted_DAG_assessdsem(fit=fits[[i]],
                                                             what = 'p_value')), edge.width=1, type="width",
                                   text_size=4, show.legend=FALSE, colors=c('black', 'black'),
                                   arrow = grid::arrow(type='closed', 18, grid::unit(10,'points')) ) +
      scale_x_continuous(expand = c(0.4, 0.1))
  }

  # plt <-
  #   ggarrange(plt_link,plt_pval,
  #             labels = fit_label,
  #             ncol = length(fits), nrow = 2)
  # plt <- grid.arrange(grobs =list(plt_link,plt_pval))
  # gridExtra::grid.arrange(grobs = list(plt_link,plt_pval))
  # gridExtra::grid.arrange(plt_link)
  #
  # ggarrange(for(i in seq_along(fits)){plt_link[[i]]},
  #                       labels = '',#fit_label,
  #                       ncol = length(fits), nrow = 2)
  plt <- ggpubr::ggarrange(plotlist = plt_link, nrow=2,ncol = length(fits),labels=fit_label)
  # cowplot::plot_grid(plt_link)

  return(plt)

}


#--------------------------------------------------------------#

#Improvements to think of :

#' Return a plot comparing the parameters estimated by dsem
#' @param fits a list of fits to compare
#' @param what which type of parameters are in the plot all/causal/variance/AR1
#' @return A plot
#' @export

# dsem_par_plot(fits=list(fit_ext_dsem_proj_dropidx,fit_assess_dsem_proj_dropidx,fit_ext_dsem_proj),what='all')

'dsem_par_plot' <- function(fits,what='all'){

  tmp <- NULL
  for(i in seq_along(fits)){

    tmp <- tmp %>%
      bind_rows(as_fitted_DAG_assessdsem(fit = fits[[i]],
                          out.type = 'other',
                          which.link = 'all') %>%
      mutate(model=fits[[i]]$version))
  }

  tmp <- tmp %>%
    mutate(lwr=Estimate-qnorm(0.975)*Std_Error,upr=Estimate+qnorm(0.975)*Std_Error) %>%
    mutate(param=ifelse(first!=second,'causal_link',
                        ifelse(direction==1,'AR1correlation','variance')))

  minval <- min(tmp %>% filter(param!='variance') %>% pull(lwr))
  maxval <- max(tmp %>% filter(param!='variance') %>% pull(upr))
  if(abs(minval)<abs(maxval)){interval <- round(abs(maxval)+0.2,digits = 1)}else{interval <- round(abs(minval)+0.2,digits = 1)}#ceiling(abs(maxval))

  maxval_var <- max(tmp %>% filter(param=='variance') %>% pull(upr))

  plt1 <- tmp  %>%
    filter(param=='causal_link') %>%
    ggplot(aes(x=path,y=Estimate,ymin=lwr,ymax=upr,col=model))+
    geom_pointrange(position=position_dodge(width=1)) +
    facet_wrap(~param)+
    geom_hline(yintercept = 0)+
    ylim(-interval,interval)+
    coord_flip()+theme(legend.position='none')

  plt2 <- tmp %>%
    filter(param=='AR1correlation') %>%
    ggplot(aes(x=path,y=Estimate,ymin=lwr,ymax=upr,col=model))+
    geom_pointrange(position=position_dodge(width=1)) +
    facet_wrap(~param)+
    geom_hline(yintercept = 0)+
    ylim(-interval,interval)+
    coord_flip()+theme(legend.position='bottom')

  plt3 <- tmp %>%
    filter(param=='variance') %>%
    ggplot(aes(x=path,y=Estimate,ymin=lwr,ymax=upr,col=model))+
    geom_pointrange(position=position_dodge(width=1)) +
    facet_wrap(~param)+
    geom_hline(yintercept = 0)+
    ylim(0,maxval_var)+
    coord_flip()+theme(legend.position='none')

  plt <- ggarrange(plt1,plt2,plt3,ncol=3)#grid.arrange
  if(what=='all'){return(plt)}
  else if(what== 'causal'){return(plt1+theme(legend.position='bottom'))}
  else if(what== 'AR1'){return(plt2)}
  else if(what== 'variance'){return(plt3+theme(legend.position='bottom'))}

}

#--------------------------------------------------------------------#


#Improvements to think of :

#' Return a plot showing the estimated imput time series
#' @param fits a list of fits to compare
#'
#' @return A plot
#' @export

# fits=list(fit_ext_dsem_proj_dropidx,fit_assess_dsem_proj_dropidx)
# TS_years = 1970:(2023+15)
#
# 'TS_plot' <- function(fits,TS_years){
#
#   tmp <- NULL
#
#   for(i in seq_along(fits)){
#
#     #get & format x_tj
#     par_tmp = fits[[i]]$obj$env$parList()
#     X_long <- par_tmp$x_tj %>% as.data.frame() %>%
#       mutate(Year= TS_years,version=fits[[i]]$version) %>%
#       pivot_longer(cols=-c(Year,version),values_to = 'estimates')
#
#     #get & format y_tj
#     dat_tmp = if(any(class(fits[[i]])=='dsem')){fits[[i]]$tmb_inputs$data$y_tj
#     }else{fits[[i]]$input$dat$y_tj}
#
#     Y_long <- dat_tmp %>% as.data.frame() %>%
#       mutate(Year= TS_years,version=fits[[i]]$version) %>%
#       pivot_longer(cols=-c(Year,version),values_to = 'data')
#
#     XY_long <- X_long %>% full_join(Y_long,join_by(Year, version, name)) %>%
#       pivot_longer(cols=c(data,estimates),names_to='value_type')
#
#     tmp <- tmp %>%
#       bind_rows(XY_long)
#
#   }
#
#   tmp %>%
#     ggplot(aes(x=Year,))
#
# }

#--------------------------------------------------------------------#


#Improvements to think of : get rid of the family arg -> can be deduced from fit

#' Perform a simulation analysis
#' @param
#'
#' @return
#' @export
#'
#'
#'

# fit_assess_dsem = fit_assess_dsem_DAG1_sd0.1
# fit_assess_dsem_DAG1_sd0.1$path <-'source'
# fit_assess_dsem$path <-'source'
# fit_assess_dsem_DAG1_sd0.1$input$path <-'source'
# # input_int_dsem4simu$path <-'source'
# sem=sem_DAG1
# #
# trysimu2 <- simulate_AssessDsem(sem=sem_DAG1,fit_assess_dsem =fit_assess_dsem_DAG1_sd0.1,
#                     family= family <- rep('normal',ncol(ESPdata_DAG1_rec)),nsims = 1)

# get_std(fits=re_fit_assess_dsem) %>%
#   filter(name %in%c('Espawnbio','mean_log_recruit','recruit','log_recruit','beta_z'))


'simulate_AssessDsem' <- function(sem, #ESPdata -> not needed?
                                  fit_assess_dsem,
                                  nsims=2,
                                  family,
                                  dsem_sim_control=list(resimulate_gmrf=TRUE,variance='none',ignoreNAs=TRUE),
                                  simverbose=TRUE,
                                  fit_sim=TRUE,seed=5){

  start_time <- Sys.time()
  set.seed(seed)
  re_fit_assess_dsem <- list()
  ### Simulation step 1: fit + simulate external dsem ---------

  # built objects we need for a dsem fit
  ext_dsem_data <- fit_assess_dsem$input$dat$y_tj
  # ext_dsem_data[1:nrow(ext_dsem_data),1] <- fit_assess_dsem$sd %>% filter(name =='x_tj') %>%

  ext_dsem_par = fit_assess_dsem$parList[c("beta_z","lnsigma_j","mu_j","delta0_j","x_tj")]
  ext_dsem_map = fit_assess_dsem$input$map[c("x_tj","lnsigma_j","mu_j")]
  if(!is.null(fit_assess_dsem$input$map[c("beta_z")])){ext_dsem_map$beta_z <-fit_assess_dsem$input$map$beta_z }
  ext_dsem_map$mu_j <- factor(c(1:ncol(ext_dsem_data)))

  # ext_dsem_family = fit_assess_dsem$input$dat$familycode_j
  ext_dsem_control=dsem_control(use_REML=F, run_model=FALSE,quiet=TRUE,
                                getJointPrecision = TRUE,
                                parameters=ext_dsem_par,
                                map=ext_dsem_map)
  # fit
  fit_ext_dsem4sim = dsem( sem=sem, tsdata=ext_dsem_data, family=family, control=ext_dsem_control)
  class(fit_ext_dsem4sim) <- 'dsem'

  #simulate with ext dsem
  sim_ext_dsem <- simulate(object = fit_ext_dsem4sim,nsim=nsims,seed=seed,
                           variance = dsem_sim_control$variance,
                           resimulate_gmrf = dsem_sim_control$resimulate_gmrf,
                           ignoreNAs = dsem_sim_control$ignoreNAs,
                           full = TRUE)

  for(i in 1:nsims){
    ### Simulation step 2: build a intDsem obj + simulate -----------

    #create a new input object
    input_int_dsem4simu <- fit_assess_dsem$input

    # stick X_tj from ext_dsem simulation
    input_int_dsem4simu$pars$x_tj <- sim_ext_dsem[[i]]$z_tj

    #fix something to avoid NAs
    input_int_dsem4simu$dat$multN_fsh[1] <- 1

    #build obj
    # input_int_dsem4simu$path <-'source'
    fit_int_dsem4simu <- fit_pk(input=input_int_dsem4simu,newtonsteps=1,
                                do.fit=FALSE,filename = NULL,
                                getsd=FALSE,save.sdrep = FALSE)

    #simulate with it
    sim_int_dsem <- fit_int_dsem4simu$simulate()

    ### Simulation step 3: new intDsem obj + fit -------

    #create another input object
    input_int_dsem4fit <- fit_assess_dsem$input
    input_int_dsem4fit$version <- i

    #stack the simulated data in this input object
    for (x in 1:length( input_int_dsem4fit$dat)){
      for(z in 1:length(sim_int_dsem)){
        if(names(input_int_dsem4fit$dat)[x]==names(sim_int_dsem)[z]){
          input_int_dsem4fit$dat[x] <-  sim_int_dsem[z]}
      }
    }
    # add the simulated y_tj from ext dsem simulation
    input_int_dsem4fit$dat$y_tj <- sim_ext_dsem[[i]]$y_tj
    input_int_dsem4fit$dat$y_tj[,1]<- NA #set NAs for the recdev col

    #fit
    # input_int_dsem4fit$path <-'source'
    re_fit_assess_dsem[[i]] <- suppressMessages(fit_pk(input=input_int_dsem4fit,newtonsteps=2,
                                                       do.fit=fit_sim,filename = NULL,
                                                       control=list(trace=0),verbose=FALSE,
                                                       save.sdrep = FALSE,getsd=TRUE))#can probably be turn to false

    #add sim data
    re_fit_assess_dsem[[i]]$simu_int_dsem <-  sim_int_dsem
    re_fit_assess_dsem[[i]]$simu_ext_dsem <- sim_ext_dsem[[i]]
    # re_fit_assess_dsem[[i]]$simu <- c(sim_ext_dsem,sim_int_dsem)
    #if the model has not be fitted, add manually some slot usefull
    if(!fit_sim){
    re_fit_assess_dsem[[i]]$input <- input_int_dsem4fit;
    re_fit_assess_dsem[[i]]$modfile<-input_int_dsem4fit$modfile;
    re_fit_assess_dsem[[i]]$path<-input_int_dsem4fit$path
    re_fit_assess_dsem[[i]]$version <- fit_assess_dsem$version}
  }

  end_time <- Sys.time()
  if(simverbose){print(end_time - start_time)}
  return(re_fit_assess_dsem)
}


#--------------------------------------------------------------------#

#Improvements to think of : add parallelization capabilities? (interest?)

#' Fit a sequence of extDsem models to peeled data set +
#' project based on these for skill prediction analysis
#'
#' @export

# retr <- retro_proj_analysis_extDsem(fit=fit_dsem4int_DAG1_sd0.1,peels=0:3,ny_proj = 3)
# ret <- fit_retro_proj_extDsem(fit=fit_dsem4int_DAG1_sd0.1,peel=5,ny_proj = 3)

'retro_proj_analysis_extDsem' <- function(fit,peels=0:5,ny_proj=3,new_control=NULL){
  retros <- lapply(peels, function(i) fit_retro_proj_extDsem(fit, peel=i,ny_proj=ny_proj,new_control=new_control))
  return(retros)
}

#--------------------------------------------------------------------#

#Improvements to think of : detail cases for which all family are not normal + if you want to change the controls

#' Internal wraper: Fit 1 extDsem models with peeled data set + project based on it
#'
#' @export


'fit_retro_proj_extDsem' <- function(fit,peel,ny_proj,new_control=NULL){

  dat <- fit$tmb_inputs$data$y_tj
  ny <- nrow(dat)
  # trim y_tj + add number of year for projection
  dat2 <- dat[1:(ny-peel),] %>% as.data.frame() %>% add_row(recdevs=rep(NA,ny_proj))

  if(is.null(new_control)){
    control1 <- dsem_control(nlminb_loops = fit$internal$control$nlminb_loops,
                             newton_loops=fit$internal$control$newton_loops,
                             trace=0,
                             eval.max = fit$internal$control$eval.max,
                             iter.max = fit$internal$control$iter.max,
                             getsd=fit$internal$control$getsd,
                             quiet=TRUE, run_model=FALSE,
                             gmrf_parameterization=fit$internal$control$gmrf_parameterization,
                             use_REML=fit$internal$control$use_REML,
                             profile=fit$internal$control$profile,
                             getJointPrecision=fit$internal$control$getJointPrecision,
                             extra_convergence_checks=fit$internal$control$extra_convergence_checks)

    fit_null <- dsem( sem=fit$internal$sem, tsdata=ts(dat2), family=fit$internal$family, control=control1)
    new_map = fit_null$tmb_inputs$map
    new_pars <- fit_null$tmb_inputs$parameters

    if(any(fit$internal$family=='normal')){ #needs to fix the sd parameter estimation
      new_map$lnsigma_j <- factor(rep(NA, length=length(new_map$lnsigma_j)))
      new_pars$lnsigma_j <- rep(log(0.1), length=length(new_pars$lnsigma_j))}

    control2 <- dsem_control(nlminb_loops = fit$internal$control$nlminb_loops,
                            newton_loops=fit$internal$control$newton_loops,
                            trace=0,
                            eval.max = fit$internal$control$eval.max,
                            iter.max = fit$internal$control$iter.max,
                            getsd=fit$internal$control$getsd,
                            quiet=TRUE, run_model=TRUE,
                            gmrf_parameterization=fit$internal$control$gmrf_parameterization,
                            use_REML=fit$internal$control$use_REML,
                            profile=fit$internal$control$profile,
                            parameters=new_pars,
                            map=new_map,
                            getJointPrecision=fit$internal$control$getJointPrecision,
                            extra_convergence_checks=fit$internal$control$extra_convergence_checks)


  }else{
      control2 = new_control}
  fit2 <- dsem( sem=fit$internal$sem, tsdata=ts(dat2), family=fit$internal$family, control=control2)
  return(fit2)
}


# repretros <- do_skill_test(reps=1:20,sem=sem_DAG1,fit=fit,peels=0:10,ny_proj=3,env_data='idealistic',parallel=TRUE)

# 'do_skill_test' <- function(reps=1:20,sem,fit,peels=0:10,ny_proj=3,env_data='idealistic',parallel=TRUE){
#   if(class(fit)[1]!='pkfit')
#     stop("fit argument is not a fitted model")
#   if(parallel){
#     message("Preparing parallel session..")
#     if(!require(snowfall))
#       stop("snowfall package required for parallel execution")
#     sfInit(parallel=TRUE, cpus=parallel::detectCores()-1)
#     sfLibrary(GOApollock);sfLibrary(dsem)
#     sfExport(c('retro_proj_analysis_AssessDsem','fit_retro_proj_AssessDsem','simulate_AssessDsem'))
#     rep_retros <- sfLapply(reps, function(i) retro_proj_analysis_AssessDsem(reps=i,sem=sem,fit=fit,peels=peels,ny_proj = ny_proj,env_data = env_data))
#     sfStop()
#   } else {
#     rep_retros <- lapply(reps, function(i) retro_proj_analysis_AssessDsem(reps=i,sem=sem,fit=fit,peels=peels,ny_proj = ny_proj,env_data = env_data))
#   }
#   return(rep_retros)
# }



'SIMretro_proj_analysis_AssessDsem' <- function(reps=1,sem,fit,peels=0:2,ny_proj=3,env_data='idealistic'){
  #start_time <- Sys.time()

  if(env_data != 'real'){ # need to run a simulation (1 for all the peels)
    seed = reps
    sim <- suppressMessages(simulate_AssessDsem(sem=sem,fit_assess_dsem =fit,family= family <- rep('normal',ncol(ESPdata_DAG1_rec)),
                                                nsims = 1,fit_sim=FALSE,simverbose = FALSE,seed=seed))}

  retros <- sim

  # add simulated x_tj in the retros
  # sim_ext_dsem[[1]]$  sim[[1]]$sd %>% filter(name%in%c('recruit','recruit_proj'))
  for(p in peels){
    retros[[p+1]]$reps = reps
    retros[[p+1]]$Rtrue <-c(sim[[1]]$simu_int_dsem$recruit,sim[[1]]$simu_int_dsem$recruit_proj)}#[(fit$input$dat$endyr-fit$input$dat$styr+1-p)]}


  end_time <- Sys.time()
  print(end_time - start_time)

  return(retros)
}


#--------------------------------------------------------------------#

#Improvements to think of : add parallelization capabilities? (interest?)

#' Fit a sequence of AssessDsem models to peeled data set +
#' project based on these for skill prediction analysis
#'
#' @export

# fit = fit_assess_dsem_DAG1_sd0.1
# fit$path <-'source'
# fit$input$path <-'source'
# retros2 <- retro_proj_analysis_AssessDsem(sem=sem_DAG1noRlink,
#                                fit=fit,peels=0:10,ny_proj = 3,env_data = 'idealistic')
# retros <- retro_proj_analysis_AssessDsem(sem=sem_DAG1,
#                                           fit=fit,peels=0:1,ny_proj = 3,env_data = 'idealistic')

'retro_proj_analysis_AssessDsem' <- function(reps=1,sem,fit,peels=0:2,ny_proj=3,env_data='idealistic'){
  start_time <- Sys.time()

  # Fit with the exact same model if no simulation  ------------
  if(env_data == 'real'){
    retros <- lapply(peels, function(i) fit_retro_proj_AssessDsem(sem=sem,fit_assess_dsem=fit,
                                                                peel=i,ny_proj=ny_proj,
                                                                env_data=env_data,reps=reps))

  # Simulate data if needed  (1 for all the peels) ---------
  }else{
    seed = reps
    #trick for a stronger sigmaR
    # ids <-  as.data.frame(fit$sem_full)
    # id_Rsigma <- which(ids$name=='sigmaR')
    # fit$parList$beta_z[id_Rsigma] <- .1
    sim <- suppressMessages(simulate_AssessDsem(sem=sem,fit_assess_dsem =fit,family= family <- rep('normal',ncol(ESPdata_DAG1_rec)),
                        nsims = 1,fit_sim=FALSE,simverbose = FALSE,seed=seed,
                        dsem_sim_control=list(resimulate_gmrf=TRUE,variance='none',
                                              ignoreNAs = if(env_data=='realistic'){FALSE}else{TRUE})))


    ## Fit with the simulated data ------
    retros2 <- lapply(peels, function(i) fit_retro_proj_AssessDsem(sem=sem,fit_assess_dsem=sim[[1]],
                                                                  peel=i,ny_proj=ny_proj,reps=reps))
    ## Add true recruitment values in every peel -----
    # /!\ 2 fits nested in each attribute of the list
    for(p in peels){
      retros2[[p+1]][[1]]$reps <- retros2[[p+1]][[2]]$reps <- reps
      retros2[[p+1]][[1]]$Rtrue <-c(sim[[1]]$simu_int_dsem$recruit,sim[[1]]$simu_int_dsem$recruit_proj)#[(fit$input$dat$endyr-fit$input$dat$styr+1-p)]}
      retros2[[p+1]][[2]]$Rtrue <- retros2[[p+1]][[1]]$Rtrue
      retros2[[p+1]][[1]]$env_data <- retros2[[p+1]][[2]]$env_data <- env_data}

    ## Unlist the nested list
    retros <- unlist(retros2,recursive = FALSE)

  }


  end_time <- Sys.time()
  print(end_time - start_time)

  return(retros)
}


#-------------------------------------------------------------------------------#


#Improvements to think of : detail cases for which all family are not normal + if you want to change the controls + more elegant way to not pass dsem dat/par/map in cole's functions

#' Internal wraper: Fit 1 extDsem models with peeled data set + project based on it
#'
#' @export


'fit_retro_proj_AssessDsem' <- function(sem,fit_assess_dsem,peel,ny_proj=3,reps=1){

  stopifnot(peel>=0)

  ## peel assessment data -----
  control_assess <- list(eval.max=10000, iter.max=10000, trace=0)

  dat2_assess <- peel_data(dat=fit_assess_dsem$input$dat[-c(105:109)], peel)#peel_data(dat=fit_assess_dsem$obj$env$data, peel)#
  # attributes(dat2_assess) <- attributes(fit_assess_dsem$obj$env$data)
  # attributes(dat2_assess)$check.passed <- NULL
  pars2_assess <- peel_pars(pars=fit_assess_dsem$input$pars[-c(40:44)], dat=dat2_assess, peel)#peel_pars(pars=fit_assess_dsem$obj$env$parList(), dat=dat2_assess, peel)#
  map2_assess <- peel_map(map=fit_assess_dsem$input$map[-c(18:20)], pars2_assess)#peel_map(map=fit_assess_dsem$obj$env$map, pars2_assess)#
  if(dat2_assess$endyr<2013){
    warning("No survey 6 data so mapping off log_q6")
    map2_assess$log_q6 <- factor(NA)
  }
  input2_assess <- list(version=paste0(fit_assess_dsem$version,'_rep',reps,'_peel',peel), path=fit_assess_dsem$path,
                 modfile=fit_assess_dsem$modfile,
                 dat=dat2_assess, pars=pars2_assess, map=map2_assess, random=fit_assess_dsem$input$random)

  ## peel dsem data -----
  dat_dsem <- fit_assess_dsem$input$dat$y_tj
  #here we want to account the nb of proj year already included in y_tj
  ny <- nrow(dat_dsem)-length(fit_assess_dsem$input$dat$Ftarget)

  # trim y_tj + add number of year for projection
  dat2_dsem <- dat_dsem[1:(ny-peel),] %>% as.data.frame() %>% add_row(recdevs=rep(NA,ny_proj))

  # create new control for the null fit
  control1_dsem <- dsem_control( trace=0,quiet=TRUE, run_model=FALSE, use_REML=FALSE)
  #make a null fit so all the parameters and map are trimmed to new y_tj dimensions
  fit_dsem_null <- dsem(sem=sem, tsdata=ts(dat2_dsem), family=names(fit_assess_dsem$input$dat$familycode_j), control=control1_dsem)

  ## gather everything -----
  input2_assessDsem <- input2_assess

  #dat
  # input2_assessDsem$dat <- c(input2_assess$dat, ) #this doesnot work, will have to replace dsem input manually
  input2_assessDsem$dat$options <- fit_dsem_null$tmb_inputs$dat$options
  input2_assessDsem$dat$RAM <- fit_dsem_null$tmb_inputs$dat$RAM
  input2_assessDsem$dat$RAMstart <- fit_dsem_null$tmb_inputs$dat$RAMstart
  input2_assessDsem$dat$familycode_j <- fit_dsem_null$tmb_inputs$dat$familycode_j
  input2_assessDsem$dat$y_tj <- fit_dsem_null$tmb_inputs$dat$y_tj
  input2_assessDsem$dat$y_tj[,1] <- NA

  input2_assessDsem$dat$Ftarget <- rep(0.2630716,ny_proj) # set duration of projections
  input2_assessDsem$dat$indxsurv_log_sd4 <- input2_assess$dat$indxsurv_log_sd4*5 #decrease weight of some indices
  input2_assessDsem$dat$indxsurv_log_sd5 <- input2_assess$dat$indxsurv_log_sd5*5 #decrease weight of some indices

  #pars
  # input2_assessDsem$pars <- c(input2_assess$pars, fit_dsem_null$tmb_inputs$parameters) #idem
  input2_assessDsem$pars$beta_z <- fit_assess_dsem$input$pars$beta_z#fit_dsem_null$tmb_inputs$parameters$beta_z
  input2_assessDsem$pars$x_tj <- fit_dsem_null$tmb_inputs$parameters$x_tj
  input2_assessDsem$pars$mu_j <- fit_dsem_null$tmb_inputs$parameters$mu_j
  input2_assessDsem$pars$delta0_j <- fit_dsem_null$tmb_inputs$parameters$delta0_j
  input2_assessDsem$pars$lnsigma_j <- rep(log(0.1), length=length(fit_dsem_null$tmb_inputs$parameters$lnsigma_j)) #fixed value

  #map
  # input2_assessDsem$map <- c(input2_assess$map, fit_dsem_null$tmb_inputs$map) #same
  input2_assessDsem$map$x_tj <- fit_dsem_null$tmb_inputs$map$x_tj
  input2_assessDsem$map$mu_j <- factor(c(NA,2:ncol(dat2_dsem))) # map off recdev mean since already a parameter
  input2_assessDsem$map$lnsigma_j <-factor(rep(NA, length=length(fit_dsem_null$tmb_inputs$map$lnsigma_j))) #noT estimated
  if(!is.null(fit_assess_dsem$input$map[c("beta_z")])){input2_assessDsem$map$beta_z <-fit_assess_dsem$input$map$beta_z }
  input2_assessDsem$map$sigmaR <- NULL ## took out of model

  #random #should not have to change this, because it was already an internal dsem fit
  # input2_assessDsem$random <-  fit_dsem_null$tmb_inputs$random#c(input2_assess$random, fit_dsem_null$tmb_inputs$random)

  ## Fit peel ------
  message("Starting optimization for peel=",peel)
  fit2 <- suppressWarnings(fit_pk(input=input2_assessDsem, getsd=TRUE, control=control_assess))
  # fit3 <- (fit_pk(input=input2_assessDsem, getsd=TRUE, control=control_assess))#,newtonsteps = 0))

  # Fit with the 'complementary' one - no Rlink if cond has, Rlink if cond doesn't have it ----------
  input2_assessDsem3 <- input2_assessDsem
  input2_assessDsem3$version <- paste0('oppRlink',input2_assessDsem3$version )

  #find where the Rlink is
  ids <-  as.data.frame(fit_dsem_null$sem_full)
  id_Rlink <- which(ids$first!=ids$second & ids$second=='recdevs')

  ### If it was a no Rlink fit ------
  if(any(is.na(input2_assessDsem$map$beta_z))){
    # re set the Rlink to be estimated
    input2_assessDsem3$map$beta_z <- factor(1:length(input2_assessDsem3$map$beta_z))
    input2_assessDsem3$pars$beta_z[id_Rlink] <- 0.01

  }else{ ### If it was a Rlink fit -------
  # Set the par to 0 and map it out
    input2_assessDsem3$map$beta_z <- 1:length(input2_assessDsem$pars$beta_z )
    input2_assessDsem3$map$beta_z[id_Rlink] <- NA
    input2_assessDsem3$map$beta_z <- as.factor(input2_assessDsem3$map$beta_z)
    input2_assessDsem3$pars$beta_z[id_Rlink] <- 0   ## checl if that do the job!
  }

  fit3 <- suppressWarnings(fit_pk(input=input2_assessDsem3, getsd=TRUE, control=control_assess))

  fits_w_wo_Rlink <- list(fit2,fit3)
  return(fits_w_wo_Rlink)
}
