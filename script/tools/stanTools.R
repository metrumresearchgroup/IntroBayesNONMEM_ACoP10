## 5/20/2008 WRG: Fixed bug in save.model. Added n.chains argument.
## 7/2/2010 WRG: Added na.rm = TRUE to summary(x1) statement
## 4/21/2014 WRG: renamed using camel case; added bugsChain

mcmc.history = function(sims.array){
  
  n1 = dim(sims.array)[1]
  n2 = dim(sims.array)[2]
  n3 = dim(sims.array)[3]
  x = data.frame(value=as.vector(sims.array),
                 chain=rep(rep(1:n2,ea=n1),n3),
                 parameter=rep(dimnames(sims.array)[[3]],ea=n1*n2))
  print(xyplot(value~rep(1:n1,n2*n3)|parameter,x,groups=chain,
               panel=panel.superpose,type="l",col=c(1,2,3,4),
               layout=c(1,6),scales=list(cex=1,y=list(relation="free")),
               xlab=list(label="sample",cex=1.2),ylab=list(label="value",cex=1.2),
               par.strip.text=list(cex=1),strip = function(...) strip.default(..., style = 1)))
  NULL
}


mcmcHistory = function(simsTable, nParPerPage = 6){
  posterior <- reshape::melt(simsTable, id.vars = c("chain", "iteration"))
  posterior <- posterior %>% dplyr::rename(parameter = variable)
  posterior$parameter <- as.character(posterior$parameter)
  
  parameters <- sort(unique(posterior$parameter))
  nParameters <- length(parameters)
  nPages <- ceiling(nParameters / nParPerPage)
  parameters <- data.frame(parameter = parameters,
                           page = sort(rep(1:nPages, length = nParameters)),
                           stringsAsFactors = FALSE)
  posterior <- metrumrg::stableMerge(posterior, parameters)
  posterior$chain <- as.factor(posterior$chain)
  
  for(i in 1:nPages){
    xplot <- subset(posterior, page == i)
    p1 <- ggplot2::ggplot(xplot, ggplot2::aes(x = iteration, y = value))
    print(p1 + ggplot2::aes(color = chain) + ggplot2::geom_line() + 
            ggplot2::labs(x = "iteration", y = "value") +
            ggplot2::theme(text = ggplot2::element_text(size = 12), axis.text = ggplot2::element_text(size = 12),
                           legend.position = "none", strip.text = ggplot2::element_text(size = 8)) +
            ggplot2::facet_wrap(~ parameter, ncol = 1, scales = "free_y"))
  }
  NULL
}


mcmc.density = function(sims.array,n=50){
  
  # n = number of points at which the density is calculated
  
  n1 = dim(sims.array)[1]
  n2 = dim(sims.array)[2]
  n3 = dim(sims.array)[3]
  x = data.frame(value=as.vector(sims.array),
                 chain=rep(rep(1:n2,ea=n1),n3),
                 parameter=rep(dimnames(sims.array)[[3]],ea=n1*n2))
  x = data.frame(parameter=rep(unique(x$parameter),ea=n),
                 value=as.vector(sapply(unique(x$parameter),function(x,n,sim.list){
                   density(sim.list$value[sim.list$parameter==x],
                           n=n,na.rm=T)$x},n=n,sim.list=x)),
                 frequency=as.vector(sapply(unique(x$parameter),function(x,n,sim.list){density(
                   sim.list$value[sim.list$parameter==x],
                   n=n,na.rm=T)$y},n=n,sim.list=x)))
  
  print(xyplot(frequency~value|parameter,x,scales="free",type="l",col=1,
               layout=c(0,min(16,length(unique(x$parameter)))),
               par.strip.text=list(cex=1),strip = function(...) strip.default(..., style = 1)))
  NULL
}

mcmcDensity = function(simsTable, byChain = FALSE, nParPerPage = 16, prior = NULL){
  posterior <- reshape::melt(simsTable, id.vars = c("chain", "iteration"))
  posterior <- posterior %>% dplyr::rename(parameter = variable)
  posterior$parameter <- as.character(posterior$parameter)
  
  parameters <- sort(unique(posterior$parameter))
  nParameters <- length(parameters)
  nPages <- ceiling(nParameters / nParPerPage)
  parameters <- data.frame(parameter = parameters,
                           page = sort(rep(1:nPages, length = nParameters)),
                           stringsAsFactors = FALSE)
  posterior <- metrumrg::stableMerge(posterior, parameters)
  posterior$chain <- as.factor(posterior$chain)
  
  if(!is.null(prior)) prior <- metrumrg::stableMerge(prior, parameters)
  
  for(i in 1:nPages){
    xplot <- subset(posterior, page == i)
    p1 <- ggplot(xplot, aes(x = value))
    if(byChain) p1 <- p1 + aes(color = chain)
    p1 <- p1 + geom_density() + 
      labs(x = "value", y = "density") +
      theme(text = element_text(size = 12), axis.text = element_text(size = 12),
            legend.position = "none", strip.text = element_text(size = 8)) +
      facet_wrap(~ parameter, ncol = 4, nrow = 4, scales = "free")
    if(!is.null(prior))
      p1 <- p1 + geom_line(data = subset(prior, page == i), aes(x = value, y = density),
                           color = "red")
    print(p1)
  }
  NULL
}


summary.mcmc.list <- function (object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), 
                               ...) 
{
  x <- mcmc.list(object)
  statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
  varstats <- matrix(nrow = nvar(x), ncol = length(statnames), 
                     dimnames = list(varnames(x), statnames))
  xtsvar <- matrix(nrow = nchain(x), ncol = nvar(x))
  if (is.matrix(x[[1]])) {
    for (i in 1:nchain(x)) for (j in 1:nvar(x)) xtsvar[i, 
                                                       j] <- coda:::safespec0(x[[i]][, j])
    xlong <- do.call("rbind", x)
  }
  else {
    for (i in 1:nchain(x)) xtsvar[i, ] <- coda:::safespec0(x[[i]])
    xlong <- as.matrix(x)
  }
  xmean <- apply(xlong, 2, mean, na.rm = TRUE)
  xvar <- apply(xlong, 2, var, na.rm = TRUE)
  xtsvar <- apply(xtsvar, 2, mean, na.rm = TRUE)
  varquant <- t(apply(xlong, 2, quantile, quantiles, na.rm = TRUE))
  varstats[, 1] <- xmean
  varstats[, 2] <- sqrt(xvar)
  varstats[, 3] <- sqrt(xvar/(niter(x) * nchain(x)))
  varstats[, 4] <- sqrt(xtsvar/(niter(x) * nchain(x)))
  varquant <- drop(varquant)
  varstats <- drop(varstats)
  out <- list(statistics = varstats, quantiles = varquant, 
              start = start(x), end = end(x), thin = thin(x), nchain = nchain(x))
  class(out) <- "summary.mcmc"
  return(out)
}

parameter.plot.table = function(parameter.array){
  # create history, density and Gelman-Rubin-Brooks plots
  # return value is a table of summary stats for the parameters
  
  require(rstan)
  
  ## create history, density and Gelman-Rubin-Brooks plots
  mcmc.history(parameter.array)
  mcmc.density(parameter.array,n=50)
  if(length(dim(parameter.array))<3){
    n.chains = 1
  }else{
    n.chains = dim(parameter.array)[2]
  }
  x1 = mcmc.list(lapply(1:n.chains,function(i) mcmc(parameter.array[,i,]))) # format required by CODA
  try(gelman.plot(x1, ask = FALSE),TRUE)
  
  ## summary stats on parameters
  ptable <- rstan:::monitor(parameter.array, warmup = 0, print = FALSE)
  
  ptable
}

parameterTable <- function(simsTable){
  
  nPerChain <- table(simsTable$chain)
  if(any(nPerChain != nPerChain[1])) stop("ERROR in parameterTable: Unequal chain lengths")
  
  `.` <-plyr::`.`
  simsArray <- plyr::daply(simsTable, .(chain),
                           function(x, iterMax){
                             as.matrix(x[,-(1:2)])
                           })
  rm(.)
  if(length(dim(simsArray)) == 2){
    oldNames <- dimnames(simsArray)
    dim(simsArray) <- c(dim(simsArray)[1], 1, dim(simsArray)[2])
    dimnames(simsArray) <- list(oldNames[[1]], NULL, oldNames[[2]])
  }else{
    simsArray <- aperm(simsArray, c(2, 1, 3))
  }
  rstan:::monitor(simsArray, warmup = 0, print = FALSE)
}



colVars <- function(a) {
  vars <- a[1,]
  for (n in 1:ncol(a))
    vars[n] <- var(a[,n])
  return(vars)
}

waic <- function(log_lik){
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
    matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
                    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  rbind(total, se)
}
