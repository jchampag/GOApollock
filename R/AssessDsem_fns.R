##########################
# Script: some functions modified from dsem/GOApollock for my use + new functions

# Juliette started 4/3/24
###########################

#Improvements to think of : none for now

#' Fit a model (assessment, dsem or assessment+dsem) to data
#' @param fit_type character chain indicating which model to fit
#'   "AssessOnly", "ExtDSEM" or "AssessDsem"
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

'AssessDsem_fit' <- function(fit_type = "AssessOnly", #or "ExtDSEM" or "AssessDsem"
                             fit_name = '_1' ,
                             ESPdata,
                             sem=sem,
                             sem_family =  rep('fixed', ncol(ESPdata)),
                             ny_proj=15,
                             Assess_recdevs= NULL, #if ExtDsem it it requiered to have it
                             save_fit = TRUE
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
                           getsd=TRUE, newtonsteps=1,
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
                                 getsd=TRUE, trace=100, newton_loops=0)

    fit_dsem <- dsem(sem=sem, tsdata=ts(ESPdata), family=sem_family,
                     estimate_delta0=FALSE, control=control_dsem)

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
    input_assessDsem$map$mu_j <- factor(c(NA,2:ncol(ESPdata))) #pourquoi pas 1:ncol()? on s'en fiche pe?
    input_assessDsem$pars$mu_j[1] <- 0 #already done by dsem, usefull?
    input_assessDsem$map$sigmaR <- NULL ## took out of model

    fit_assessDsem <- fit_pk(input=input_assessDsem,
                             getsd=TRUE, newtonsteps=1, do.fit=1,
                             save.sdrep = TRUE,
                             filename = if(save_fit){paste0(fit_type,fit_name,'.RDS')}else{NULL})

    fit_assessDsem <- c(fit_assessDsem,sem_full=fit_dsem$sem_full)
    class(fit_assessDsem) <- c(class(fit_assessDsem),'assessdsem')
    out <- fit_assessDsem
  }

  return(out)

}

#--------------------------------------------------------------------------#

#Improvements to think of : none for now

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
                                     fit_null_dsem=NULL,#/!\ this part could be simplified, to evaluate later
                                     lag = c(0,1), direction = c(1,2),
                                     what = "Estimate",which.link='causal',
                                     out.type="matrix")
{

  if(any(class(fit)=='dsem')){ coefs = summary(fit);  vars = colnames(fit$tmb_inputs$data$y_tj)
  }else{#cannot use summary(fit) for a 'assessdsem' output, this part replace what it does
    if(all((class(fit)!='assessdsem'))){print('This function can only work with a dsem like output')}

    model = fit_null_dsem$sem_full #/!\ this part could be simplified, to evaluate later
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

    out = list(coef = array(0, dim = rep(length(vars), 2), dimnames = list(vars,
                                                                             vars)))
    out$coef[as.matrix(coefs[, c("first", "second")])] = coefs[,
                                                                 what]
    class(out) = "fitted_DAG"
    return(out)
  }
}
#
as_fitted_DAG_assessdsem(fit=fits[[2]],#fit_assess_dsem_proj,#fit_ext_dsem_proj_dropidx,#fits[[1]],
                         fit_null_dsem = NULL,#fit_ext_dsem_proj,
                                        #fit_null_dsem=NULL,#/!\ this part could be simplified, to evaluate later
                                        lag = c(0,1), direction = c(1,2),
                                        what = "p_value",which.link='all',
                                        out.type="matrix")

#--------------------------------------------------------------------------#

#Improvements to think of : none for now

#' Return a plot comparing the log recruitment estimations for different type of fits
#' @param fits a list of fits to compare
#' @param ny_proj nb of year projected
#' @param mean_log_rec value to add to recruitment deviations from dsem fit
#' @param extDsem_name name to be passed out for the plot (only for ext dsem fits, pkfit uses version)
#' @return A plot
#' @export


'rec_plot' <- function(fits,
                       ny_proj=15,
                       mean_log_rec=0,
                       extDsem_name=NULL)


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
                           version=ifelse(is.null(extDsem_name),paste0("extDsem_",log_mean_rec_pointer),
                                          paste0("extDsem_",extDsem_name[log_mean_rec_pointer]))))
    }
  }

  plt <- tmp %>%
    ggplot(aes(x=year,y=est,ymin=upr,ymax=lwr))+
    geom_line(aes(col=version))+
    geom_ribbon(aes(fill=version),alpha=0.2)+
    labs(y="Estimated log recruitment")+
    theme(legend.position='bottom')


  return(plt)

}

# fits=list(fit_assess_usual_proj,
#                              fit_ext_dsem_proj,
#                              fit_ext_dsem_proj_dropidx)



# rec_plot(fits=list(fit_assess_usual_proj,
#                    fit_ext_dsem_proj_dropidx,
#                    fit_ext_dsem_proj_dropidx),
#          ny_proj=15,
#          mean_log_rec=c(0,1),
#          extDsem_name = c('nodropidx','dropidx')     )

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
      tmp <- tmp %>%
        bind_rows(fits[[i]]$sd %>%  filter(name=='beta_z') %>% filter(year==1970))

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
    coord_flip()+labs(x="Sigma R estimates")

  return(plt)
}

# sigmaR_plot(fits=list(fit_assess_dsem_proj,fit_assess_dsem_proj_dropidx,
#                       fit_assess_usual_proj,fit_assess_usual_proj_dropidx))


#----------------------------------------------------------------------#


#Improvements to think of : add subtitle causal-link/p-value
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

    plt_link[[i]] <- plot( (as_fitted_DAG_assessdsem(fit_null_dsem = fit_null_dsem_proj,fit=fits[[i]])),
                           edge.width=1, type="width",
                           text_size=4, show.legend=FALSE,
                           arrow = grid::arrow(type='closed', 18, grid::unit(10,'points')) ) +
      scale_x_continuous(expand = c(0.4, 0.1))


    # tour=tour+2
  }
  for(i in seq_along(fits)){
    plt_link[[length(fits)+i]]  <-  plot(  (as_fitted_DAG_assessdsem(fit_null_dsem = fit_null_dsem_proj,fit=fits[[i]],
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
#' @return A plot
#' @export
#'
'dsem_par_plot' <- function(fits){



}
