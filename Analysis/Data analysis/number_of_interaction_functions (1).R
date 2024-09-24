## functions to select the best E

best_E_fct <- function(data, target, interactor, libsize, step=10, maxE=6){
  libsize <- paste(libsize, libsize, step)
  max1 <- c()
  max2 <- c()
  for(i in 1:maxE){
    ccm_out <- CCM(dataFrame = data, target = target, column = interactor,
                   libSizes = libsize, Tp = -1, E = i, sample = 10)
    # max1[i-1] <- ccm_out[1,2]
    max2[i] <- ccm_out[1,3]
  }
  df <- data.frame(target=target,
                   interactor=interactor,
                   Best_E = which.max(max2))
  return(df)
}

CCM_convergence_test <- function(data, target, interactor, libsize, step=3, E, lmin=NULL, lmax=NULL, maxlag=4){
  seeds <- sample(10000:99999, 100, replace = F)
  if(is.null(lmin)){
    lmin <- floor(nrow(data)*0.2)
  }
  lmin <- max(lmin,E+2)
  
  if(is.null(lmax)){
    lmax <- nrow(data) - (maxlag+E-2) 
  }
  
  rhos_bool <- sapply(seeds, function(seed){
    ccm_out <- CCM(dataFrame = data, target = target, column = interactor,
                   libSizes = paste(lmin, lmax, lmax-lmin), Tp = 0, E = E,
                   sample = 1, seed = seed)
    rho_lmin <- ccm_out[1,3]
    rho_lmax <- ccm_out[2,3]
    return(rho_lmin >= rho_lmax)
  })
  
  pvalue <- sum(rhos_bool)/length(rhos_bool)
}


CCM_wrapper <- function(data, target, interactor, libsize, step=3, E, lmin=NULL, lmax=NULL, maxlag=4){
  libsizes <- paste(E+2, libsize, step)
  ccm_out <- CCM(dataFrame = data, target = target, column = interactor,
                 libSizes = libsizes, Tp = 0, E = E, sample = 100)
  rho <- ccm_out[nrow(ccm_out),3]
  plotdd <- data.frame(x=seq(E+2,libsize,by = step), y=ccm_out[,3], target=target, interactor=interactor)
  if(rho<=0) {
    return(list(df=data.frame(rho=rho, skill_sign="negative", slope=NA, convergence.pvalue=NA), plot=plotdd))
  }
  slope <- lm(ccm_out[,3] ~ seq(E+2,libsize,by = step))$coefficients[2]
  slope <- ifelse(slope>0,"positive","negative")
  convergence.pvalue <- CCM_convergence_test(data, target, interactor, libsize, step, E, lmin, lmax, maxlag)
  return(list(df=data.frame(rho=rho, skill_sign="positive", slope=slope, convergence.pvalue=convergence.pvalue), plot=plotdd))
}

CCM_significance <- function(data, target, interactor, libsize, step=3, E, surrogateMethod="seasonal", T_period=12, rho){
  ts <- unlist(data[,interactor])
  if(surrogateMethod=="seasonal") {
    surr_interactor = SurrogateData(ts, method = surrogateMethod, T_period = T_period,
                                    num_surr = 1000, alpha = 3)
  } else {
    surr_interactor = SurrogateData(ts, method = "random_shuffle",
                                    num_surr = 1000, alpha = 3)      
  }
  
  rho_surr <- data.frame(interactor_rho = numeric(1000))
  
  interactor_data = as.data.frame(cbind(data[,1], ts, surr_interactor))
  names(interactor_data) = c("time", target, paste("T", as.character(seq(1, 1000)),	sep = ""))
  
  for (i in 1:1000) {
    targetCol = paste("T", i, sep = "") 
    ccm_out = CCM(dataFrame = interactor_data, E = E, Tp = 0, columns = target,
                  target = targetCol, libSizes = paste(libsize, libsize, step), sample = 1)
    col = paste(target, ":", targetCol, sep = "")
    rho_surr$interactor_rho[i] = ccm_out[1, col]
  }
  significance <- 1 - ecdf(rho_surr$interactor_rho)(rho)
  return(significance)
}


# e.g. https://ha0ye.github.io/rEDM/articles/rEDM-time-delay-ccm.html

CCM_lag_test <- function(data, target, interactor, libsize, step=3, E, maxlag=4){
  lags <- -maxlag:maxlag
  seeds <- sample(10000:99999, 100, replace = F)
  libsize <- nrow(data) - (maxlag+E-2)
  lag_rhos_bool <- lapply(1:length(seeds), function(j){
    seed <- seeds[j]
    skills <- c()
    i <- 0
    for(l in lags){
      i <- i + 1
      ccm_out = CCM(dataFrame = data, target = target, column = interactor, random = T, seed = seed,
                    libSizes = paste(libsize, libsize, step), Tp = l, E = E, sample = 1)
      skills[i] <- ccm_out[1,3]
    }
    negLag_maxRho <- max(skills[0:maxlag+1])
    posLag_maxRho <- max(skills[(maxlag+2):(2*maxlag+1)])
    dd <- data.frame(lag=lags, seed=seed, rho = skills, target=target, interactor=interactor)
    return(list(bool = posLag_maxRho >= negLag_maxRho, df=dd))
  })
  
  plotdd <- do.call("rbind", map(lag_rhos_bool,2)) %>%
    group_by(lag, target, interactor) %>%
    summarize(mean_CCM = mean(rho),
              sd_CCM = sd(rho))
  
  lag_rhos_bool <- unlist(map(lag_rhos_bool,1))
  pvalue <- sum(lag_rhos_bool)/length(lag_rhos_bool)
  optimal.lag <- unlist(plotdd[which.max(plotdd$mean_CCM),"lag"])
  
  return(list(df = data.frame(lag.pvalue=pvalue, optimal.lag = optimal.lag),
              plot = plotdd))
}


numberInteractions <- function(data, libsize=nrow(data), step=3, maxE=6, targets, interactors, trait_variables,
                               surrogateMethod="seasonal", T_period=12, maxlag=4, lag.test=T,
                               lmin=NULL, lmax=NULL) {
  libsize <- min(nrow(data) - (maxlag+maxE-2), libsize)
  var_pairs <- expand.grid(target=targets, interactor=interactors)
  var_pairs$target <- as.character(var_pairs$target)
  var_pairs$interactor <- as.character(var_pairs$interactor)
  var_pairs <- var_pairs %>%
    dplyr::filter(!(target %in% trait_variables & interactor %in% trait_variables))
  
  CCM_list <- mclapply(1:nrow(var_pairs), function(row){
    # determine best E
    Best_E_df <- best_E_fct(data, var_pairs[row,1],var_pairs[row,2],
                            step = step, maxE = maxE, libsize = libsize)
    
    libsize.specific <- nrow(data) - (maxlag+Best_E_df$Best_E-2) 
    # Calculate CCM skill and assess convergence and sign
    ccm_out <- CCM_wrapper(data, Best_E_df$target, Best_E_df$interactor, libsize=libsize.specific,
                           step=step, E=Best_E_df$Best_E, lmin, lmax, maxlag=maxlag)
    
    rho_df <- cbind(Best_E_df, ccm_out$df)
    convergence.plot <- ccm_out$plot
    
    # Check whether CCM skill is significantly bigger than surrogate TS 
    rho_df$surrogate.pvalue <- CCM_significance(data, rho_df$target, rho_df$interactor,
                                                libsize = libsize.specific, step=step, E = rho_df$Best_E, 
                                                surrogateMethod=surrogateMethod, T_period=T_period, rho = rho_df$rho)
    
    ## lag test
    if(lag.test){
      lag_test <- CCM_lag_test(data = data, target = rho_df$target, 
                               interactor = rho_df$interactor, maxlag=maxlag,
                               libsize = libsize.specific, step=step, E = rho_df$Best_E)
      
      rho_df <- cbind(rho_df, lag_test$df)
      lag.plot <- lag_test$plot
    } else {
      lag.plot <- NULL
    }
    
    return(list(df=rho_df, convergence.plot=convergence.plot, lag.plot=lag.plot))
  }, mc.cores = detectCores()-1)
  
  pairwise_ccm_df <- do.call("rbind",map(CCM_list,1))
  convergence.plots <- map(CCM_list, 2)
  lag.plots <- map(CCM_list, 3)
  names(lag.plots) <- names(convergence.plots) <- paste(var_pairs$target,"xmap", var_pairs$interactor)
  
  return(list(ccm_interactions=pairwise_ccm_df, convergence.plots=convergence.plots, lag.plots=lag.plots))
}





plotCCM <- function(lag.dfs, convergence.dfs, ccm.dfs){
  plots.all <- list()
  l <- 1
  for(i in seq(1,length(lag.dfs),10)){
    k <- 1
    plots <- list()
    
    tar <- lag.dfs[[i]][1,"target"]
    int <- lag.dfs[[i]][1,"interactor"]
    
    p1 <- ggplot(convergence.dfs[[i]], aes(x,y)) + geom_line() + geom_point() +
      labs(title=paste0(tar," predicting \n",int,":\ndetermine whether \n",int," drives \n",tar),
           x="Library size", y="CCM skill")
    
    p2 <- ggplot(lag.dfs[[i]], aes(lag,mean_CCM)) + geom_line() + geom_point() +
      geom_pointrange(aes(ymin=mean_CCM-sd_CCM,ymax=mean_CCM+sd_CCM)) +
      labs(title=paste0(tar," predicting \n",int,":\nlag test"))
    
    plots[[k]] <- p1+theme(plot.title = element_text(size = 7))
    k <- k + 1
    plots[[k]] <- p2+theme(plot.title = element_text(size = 7))
    k <- k + 1
    
    for(j in (i+1):(i+9)){
      if(j > nrow(ccm.dfs)) break
      
      tar <- lag.dfs[[j]][1,"target"]
      int <- lag.dfs[[j]][1,"interactor"]
      
      p1 <- ggplot(convergence.dfs[[j]], aes(x,y)) + geom_line() + geom_point() +
        labs(title=paste0(tar," predicting \n",int,":\ndetermine whether \n",int," drives \n",tar),
             x="Library size", y="CCM skill")
      
      p2 <- ggplot(lag.dfs[[j]], aes(lag,mean_CCM)) + geom_line() + geom_point() +
        geom_pointrange(aes(ymin=mean_CCM-sd_CCM,ymax=mean_CCM+sd_CCM)) +
        labs(title=paste0(tar," predicting \n",int,":\nlag test"))
      
      plots[[k]] <- p1+theme(plot.title = element_text(size = 7))
      k <- k + 1
      plots[[k]] <- p2+theme(plot.title = element_text(size = 7))
      k <- k + 1
    }
    p <- Reduce("+", plots) + plot_layout(ncol = 4) + plot_annotation(title = paste("plots ", i, " to ",i+9))
    plots.all[[l]] <- p 
    l <- l + 1
  }
  return(plots.all)
}
