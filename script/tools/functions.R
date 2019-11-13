 ## Miscellaneous functions

trapez = function(x,y){
  nx = length(x)
  0.5*sum((y[-1]+y[-nx])*(x[-1]-x[-nx]))
}

pkpar = function(x,y){
  nx = length(x)
  auc = sum((x[-1]-x[-nx])*(y[-1]+y[-nx]),na.rm=T)/2
  nmax = order(y,na.last=F)[nx]
  cmax = y[nmax]
  tmax = x[nmax]
  c(auc=auc,cmax=cmax,tmax=tmax)
}

frac = function(x) sum(x)/length(x)

swap = function(x, oldval, newval)
{
  if(length(oldval) != length(newval))
    stop("length(oldval)!=length(newval)")
  x1 = match(x, oldval, 0)
  x[as.logical(x1)] = newval[x1]
  x
}

mvdnorm = function(y,mu,sigma2){
  # multivariate norma density
  d = length(mu)
  det.sigma2 = det(sigma2)
  inv.sigma2 = solve(sigma2)
  exp(-t(y-mu)%*%inv.sigma2%*%(y-mu)/2)/sqrt(((2*pi)^d)*det.sigma2)
}


# generalized logit and inverse logit functions
logit.inv = function(x, lower = 0, upper = 1)
{
  x1 <- exp(x)
  ifelse(is.infinite(x1),upper,lower + ((upper - lower) * x1)/(1 + x1))
}

logit = function(x, lower = 0, upper = 1)
{
  x <- (x - lower)/(upper - lower)
  log(x/(1 - x))
}

factor2char.data.frame = function(x){
  for(i in 1:ncol(x)){
    if(class(x[[i]])=="factor") x[[i]] = as.character(x[[i]])
  }
  x
}

strip.lead.blanks = function(x){
  w = regexpr(" +",x)
  as.character(ifelse(w==1,substring(x,attr(w,"match.length")+1),x))
}

symmat = function(a){
  # generate symmetric matrix from its lower triangle given as
  # (a[1,1],a[2,1],a[2,2],a[3,1],a[3,2],a[3,3],...)
  n = (-1+sqrt(1+8*length(a)))/2
  x = matrix(rep(0,n*n),ncol=n)
  x[t(lower.tri(x,T))] = a
  y = x + t(t(lower.tri(x))*x)
  y
}

qqhist <- function(x, label = NULL, main=NULL)
{
  # plot histogram and quantile-normal quantile plot for vector x
  
  plot1 <- histogram(~x,
                     ylab=list(cex=1.2),	xlab=list(label=label,cex=1.2),
                     scales=list(cex=1.2))
  plot2 <- qqmath(~x,prepanel = prepanel.qqmathline,
                  panel=function(x){
                    panel.qqmath(x,cex=1.2,pch=16)
                    panel.qqmathline(x,distribution=qnorm,col=3,lwd=3)
                  },
                  ylab=list(label=label,cex=1.2),
                  xlab=list(label="quantiles of standard normal",cex=1.2),
                  scales=list(cex=1.2),ylim=range(x))
  list(plot1=plot1,plot2=plot2)
}

qnorm.trunc = function(p,mean=0,sd=1,lower=-Inf,upper=Inf)
  qnorm(p*pnorm(upper,mean,sd)+(1-p)*pnorm(lower,mean,sd),mean,sd)

rnorm.trunc = function(n,mean=0,sd=1,lower=-Inf,upper=Inf)
  qnorm.trunc(runif(n),mean,sd,lower,upper)

# Prepanel and panel plot functions for error bars taken from R help response by
# Deepayan Sarkar
prepanel.ci <- function(y, ly, uy, subscripts, ...)
{
  y <- as.numeric(y)
  ly <- as.numeric(ly[subscripts])
  uy <- as.numeric(uy[subscripts])
  list(ylim = range(y, uy, ly, finite = TRUE))
}

panel.ci <- function(x, y, ly, uy, subscripts, pch = 16, col.line =
                       'black', ...)
{
  x <- as.numeric(x)
  y <- as.numeric(y)
  ly <- as.numeric(ly[subscripts])
  uy <- as.numeric(uy[subscripts])
  panel.arrows(x, ly, x, uy, col = col.line,
               length = 0.25, unit = "native",
               angle = 90, code = 3)
  panel.xyplot(x, y, pch = pch, col.line = col.line, ...)
}

mbmaPpcPlot <- function(data, x = "x", treatment = "treatment", type = "type", xlab = "x"){
  
  names(data)[names(data) == x] <- "x"
  names(data)[names(data) == treatment] <- "treatment"
  names(data)[names(data) == type] <- "type"
  xobs <- data[data$type == "observed",]
  
  ## Reorder assuming that the predicted values are in the same order as the observed values.
  xOrder <- order(xobs$treatment, xobs$x)
  xOrder <- rep(xOrder, 4) + rep(0:3, ea = length(xOrder)) * length(xOrder)
  data <- data[xOrder,]
  xobs <- data[data$type == "observed",]
  
  ## Calculate number of treatments and number of arms of each treatment
  treats <- unique(xobs$treatment)
  nTreat <- length(treats)
  nPerTreat <- sapply(treats, function(treat, data)
    sum(data$treatment == treat),
    data = xobs)
  
  ## Calculate number of equally spaced positions on the y axis
  ny <- sum(nPerTreat) + nTreat - 1
  dy <- 1 / (ny - 1)
  y <- (1:ny) * dy
  
  ## Drop y positions between treatments to create gaps
  iSkip <- cumsum(nPerTreat[-nTreat]) + 1:(nTreat - 1)
  y <- y[-iSkip]
  
  ## Calculate label locations
  iStart <- c(1, 1 + cumsum(nPerTreat[-nTreat]))
  iEnd <- cumsum(nPerTreat)
  yLab <- (y[iStart] + y[iEnd]) / 2
  
  data$y <- rep(y, 4)
  xyplot(y ~ x, data, groups = type, panel =
           function(x, y, subscripts, groups, ...){
             iObs <- (groups[subscripts] == "observed")
             iMed <- (groups[subscripts] == "median")
             i5 <- (groups[subscripts] == "5%ile")
             i95 <- (groups[subscripts] == "95%ile")
             panel.xyplot(x[iObs], y[iObs], pch = 19, col = "black", ...)
             panel.segments(x[iMed], y[iMed] - dy/2, x[iMed], y[iMed] + dy/2,
                            col = "red", lwd = 2, ...)
             panel.segments(x[i5], y[iMed], x[i95], y[iMed],
                            col = "blue", lwd = 2, ...)
           },
         scales = list(cex = 1.2, y = list(at = yLab, labels = treats)),
         xlab = xlab, ylab = "")
}

makeNONMEMScript_BAYES <- function(scriptFile, scriptStubFile, mcmcFile, inits, priors,
                                   nIter, nBurnIn, seed, print, modelName, OLKJDFvalue, AUTO){
  ## Make a NONMEM script by taking a NMTRAN script stub containing everything
  ## except the $THETA, $OMEGA, $SIGMA and $ESTIMATE statements.
  scriptStub <- scan(scriptStubFile, what="", sep="\n", blank.lines.skip = F)
  
  # OLKJDF added to $EST
  # When 0, the usual inverse Wishart prior is used for Omegas. 
  # When OLKJDF>0, then the LKJ density is used as the prior, with OLKJDF degrees of freedom for all OMEGA blocks
  estimate <- c("$EST METHOD=BAYES INTERACTION LAPLACIAN",
                "SIGL=6 NSIG=3",      ## might need to remove this, not tested with NUTS
                #paste("AUTO=", as.integer(AUTO)),
                #paste("OLKJDF=", as.integer(OLKJDFvalue)),           ## OLKJDF =0(default)
                paste("NBURN =", as.integer(nBurnIn), 
                      "NITER =", as.integer(nIter), 
                      "SEED =", seed),
                paste("PRINT =", print,
                      "NOABORT FILE =", mcmcFile))
  
  table1 <- paste(c("$TABLE ID EVID TIME DV IPRED CWRES CWRESI NPDE WT",
                    "NOPRINT ONEHEADER FILE=./.tab"
  ),
  sep = "\n")
  
  table2 <- paste(c("$TABLE ID WT CL V2 Q KA V3 ETA1 ETA2 PER FORM AUC",
                    "NOPRINT ONEHEADER FILE=./par.tab"
  ),
  sep = "\n")
  
  script <- c(scriptStub, inits, priors, estimate, table1, table2)
  if(file.exists(scriptFile)) file.remove(scriptFile)
  write(script, scriptFile)
  return(NULL)
}

makeNONMEMScript <- function(scriptFile, scriptStubFile, mcmcFile,
                             inits, priors, tables,
                             nPost, nBurn, nThin = 1, seed, print = 100,parafprint=10000,
                             OLKJDF = 1, OVARF = -1, AUTO = 1, NUTS_DELTA = 0.8,
                             method = "BAYES"){
    ## Make a NONMEM script by taking a NMTRAN script stub containing everything
    ## except the $THETA, $OMEGA, $SIGMA and $ESTIMATE statements.
    scriptStub <- scan(scriptStubFile, what="", sep="\n", blank.lines.skip = F)
    
    ## OLKJDF added to $EST
    ## When 0, the usual inverse Wishart prior is used for Omegas. 
    ## When OLKJDF>0, then the LKJ density is used as the prior, with OLKJDF degrees of freedom for all OMEGA blocks
    
    if(method == "BAYES"){
        method <- c("METHOD = BAYES INTERACTION LAPLACIAN")
##        method <- c("METHOD = BAYES INTERACTION LAPLACIAN",
##                    "SIGL=6 PSCALE_MIN=1.0E-10 PSCALE_MAX=1.0E+10",
##                    "GRD=TN(1-19) ISAMPLE_M1=5 ISAMPLE_M2=5 ISAMPLE_M3=5 PSAMPLE_M1=5")
    }else if(method == "NUTS"){
        method <- c("METHOD=NUTS INTERACTION")
##        method <- c("METHOD=NUTS INTERACTION",
##                    "SIGL=6")
    }else
        stop("method = ", method, " not supported, Must be BAYES or NUTS.")
    
    estimate <- c("$ESTIMATION",
                  method,
                  paste("AUTO=", as.integer(AUTO)),
                  "CTYPE = 0",
                  paste("OLKJDF=", as.integer(OLKJDF)),
                  paste("OVARF=", as.integer(OVARF)),
                  paste("NUTS_DELTA=", NUTS_DELTA),
                  paste("NBURN =", as.integer(nBurn * nThin), 
                        "NITER =", as.integer(nPost * nThin),
##                        "THIN =", as.integer(nThin),
                        "SEED =", seed),
                  paste('parafprint =', parafprint),
                  paste("PRINT =", print,
                        "NOABORT FILE =", mcmcFile))
        
    script <- c(scriptStub, inits, priors, estimate, tables)
    if(file.exists(scriptFile)) file.remove(scriptFile)
    write(script, scriptFile)
    return(NULL)
}

makeNONMEMScriptOpt <- function(scriptFile, scriptStubFile, mcmcFile,
                             inits, priors, tables,
                             nPost, nBurn, nThin = 1, seed, print = 100,
                             OLKJDF = 1, OVARF = -1, AUTO = 1,
                             method = 1){
    ## Make a NONMEM script by taking a NMTRAN script stub containing everything
    ## except the $THETA, $OMEGA, $SIGMA and $ESTIMATE statements.

    ## Intended for use with NONMEM optimization methods
    
    scriptStub <- scan(scriptStubFile, what="", sep="\n", blank.lines.skip = F)

    method = paste("METHOD =", method)
    
    estimate <- c("$ESTIMATION",
                  method,
                  paste("AUTO=", as.integer(AUTO)),
                  paste("OLKJDF=", as.integer(OLKJDF)),
                  paste("OVARF=", as.integer(OVARF)),
                  paste("PRINT =", print, "NOABORT"),
                  paste("LAPLACIAN INTERACTION"),
                  paste("SIGL = 12"),
                  paste("SIG = 4"))
        
    script <- c(scriptStub, inits, priors, estimate, tables)
    if(file.exists(scriptFile)) file.remove(scriptFile)
    write(script, scriptFile)
    return(NULL)
}

runChain <- function(chain, modelName, modelDir, priors, tables,
                     nPost, nBurn, nThin = 1, seed, print = 100,parafprint=10000,
                     OLKJDF = 1, OVARF = -1, AUTO = 1, NUTS_DELTA = 0.8,
                     grid = TRUE, method = "BAYES", pe = "orte 32",
                     mode = "para"){
    inits <- geninit()
    iOmega <- grep("OMEGA", inits)[1]
    inits2 <- double(length(inits) + 1)
    inits2[1:(iOmega - 1)] <- inits[1:(iOmega - 1)]
    inits2[iOmega] <- paste("(", nThin, " FIX)", sep = "")
    inits2[(iOmega + 1):length(inits2)] <- inits[iOmega:length(inits)]
    
    scriptFileRoot <- paste(modelName, ".", chain, sep = "")
    thisModelDir <- file.path(modelDir, modelName)
    makeNONMEMScript(scriptFile = file.path(thisModelDir,
                                            paste(scriptFileRoot, ".ctl", sep = "")),
                     scriptStubFile = file.path(modelDir,
                                                paste(modelName, "stub.ctl", sep = "")),
                     mcmcFile = "/dev/null",
                     ## mcmcFile = paste(scriptFileRoot, "MCMC.txt", sep = ""),
                     inits = inits2, priors = priors, tables = tables,
                     nPost = nPost, nBurn = nBurn, nThin = nThin,
                     seed = seed[chain], print = print, parafprint=parafprint,
                     OLKJDF = OLKJDF, OVARF = OVARF, AUTO = AUTO,
                     NUTS_DELTA = NUTS_DELTA,
                     method = method)
    if(mode == "para"){
      file.copy(file.path(modelDir, paste(modelName, ".pnm", sep = "")),
              file.path(thisModelDir, paste(scriptFileRoot, ".pnm", sep = "")),
              overwrite = TRUE)
    NONR(run = scriptFileRoot,
         command = "/opt/NONMEM/nm74gf/nmqual/autolog.pl",
         project = thisModelDir,
         grid = grid,
         wait = FALSE,
         diag = FALSE,
         fdata = FALSE,
         purge = TRUE,
         checkrunno = TRUE,
         pe = pe,
         mode = mode
         )
    }else{
      NONR(run = scriptFileRoot,
           command = "/opt/NONMEM/nm74gf/nmqual/autolog.pl",
           project = thisModelDir,
           grid = grid,
           wait = FALSE,
           diag = FALSE,
           fdata = FALSE,
           purge = TRUE,
           checkrunno = TRUE
      )
      
    }
}

makeNONMEMScriptOpt <- function(scriptFile, scriptStubFile, mcmcFile,
                             inits, priors, tables,
                             nPost, nBurn, nThin = 1, seed, print = 100,
                             OLKJDF = 1, OVARF = -1, AUTO = 1,
                             method = 1, Lapl = FALSE){
    ## Make a NONMEM script by taking a NMTRAN script stub containing everything
    ## except the $THETA, $OMEGA, $SIGMA and $ESTIMATE statements.

    ## Intended for use with NONMEM optimization methods
    
    scriptStub <- scan(scriptStubFile, what="", sep="\n", blank.lines.skip = F)

    method = paste("METHOD =", method)
    
    estimate <- c("$ESTIMATION",
                  method,
                  paste("AUTO=", as.integer(AUTO)),
                  paste("OLKJDF=", as.integer(OLKJDF)),
                  paste("OVARF=", as.integer(OVARF)),
                  paste("PRINT =", print, "NOABORT"),
                  ifelse(Lapl == TRUE, paste("LAPLACIAN INTERACTION"), ""),
                  paste("SIGL = 12"),
                  paste("SIG = 4"),
                  paste("MAXEVAL = 99999999"))
        
    script <- c(scriptStub, inits, priors, estimate, tables)
    if(file.exists(scriptFile)) file.remove(scriptFile)
    write(script, scriptFile)
    return(NULL)
}

runOpt <- function(modelName, modelDir, priors, tables, print = 100,
                     OLKJDF = 1, OVARF = -1, AUTO = 1,
                     method = 1, Lapl = FALSE,
                     grid = TRUE, pe = "orte 32", mode = "para"){
    inits <- geninit()
    
    scriptFileRoot <- modelName
    thisModelDir <- file.path(modelDir, modelName)
    makeNONMEMScriptOpt(scriptFile = file.path(thisModelDir,
                                            paste(scriptFileRoot, ".ctl", sep = "")),
                     scriptStubFile = file.path(modelDir,
                                                paste(modelName, "stub.ctl", sep = "")),
                     inits = inits, priors = priors, tables = tables,
                     print = print,
                     OLKJDF = OLKJDF, OVARF = OVARF, AUTO = AUTO,
                     method = method, Lapl = Lapl)
    if(mode == "para"){
    file.copy(file.path(modelDir, paste(modelName, ".pnm", sep = "")),
              file.path(thisModelDir, paste(scriptFileRoot, ".pnm", sep = "")),
              overwrite = TRUE)
    NONR(run = scriptFileRoot,
         command = "/opt/NONMEM/nm74gf/nmqual/autolog.pl",
         project = thisModelDir,
         grid = grid,
         wait = FALSE,
         diag = FALSE,
         fdata = FALSE,
         purge = TRUE,
         checkrunno = TRUE,
         pe = pe,
         mode = mode
         )
    }else{
      NONR(run = scriptFileRoot,
           command = "/opt/NONMEM/nm74gf/nmqual/autolog.pl",
           project = thisModelDir,
           grid = grid,
           wait = FALSE,
           diag = FALSE,
           fdata = FALSE,
           purge = TRUE,
           checkrunno = TRUE
      )
      
    }
}
