#----------------------------------------------------------------------------------------
# R code for running simulation in
# "Marginal analysis of multiple outcomes with informative cluster size"
# includes
# 1. Functions
# 2. Code for data generation
# 3. Fitting simulated data with
#      i. univariate GEE and CWGEE for each outcome
#      ii. multivariate GEE and CWGEE with common slope for MetS
#           - with unstructured, exchangeable, independent correlaiton structures
#      iii. multivariate GEE and CWGEE with unique slope for MetS
#             - with unstructured, exchangeable, independent correlaiton structures
# 4. Save results 
#----------------------------------------------------------------------------------------

#-----------------------------
# install packages
#-----------------------------

library(bridgedist)
library(gee)
library(MASS)
library(geepack)
library(Matrix)
library(fastDummies)
library(matrixcalc)


### parameters to set------------------------------------------------
#
# N -- number of subjects
# nteeth -- maximum number of teeth per subject
# ntime -- number of visits
# truetau -- correlation between teeth per subject
# truealpha -- correlation between viist per tooth
# truea1 -- intercept for outcome 1
# truea2 -- intercept for outcome 2
# truea3 -- intercept for outcome 3
# truebeta1 -- effect of MetS
# truebeta2 -- effect of visit
# truebeta3 -- interaction between MetS and visit
#
# RR -- random integer between 1 and 1000 to use to set the seed
#---------------------------------------------------------------------


##First read in the arguments listed at the command line
args=(commandArgs(TRUE))
print(args)

## args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.

if(length(args)==0){
  print("No arguments supplied.")
  ## supply default values
  N=100
  nteeth=28
  ntimemean=5
  truea1=-1
  truea2=-1.5
  truea3=-2
  truebeta1=0.5
  truetau=0
  truealpha12=0
  truealpha13=0.35
  truealpha23=0.7
  RR=241
}else{
  for (i in (1:length(args))) eval(parse(text=args[[i]]))
}

SS <- as.numeric(Sys.getenv("SGE_TASK_ID"))

if (is.na(SS)) SS <- 1

print(paste("N=",N))
print(paste("SS=",SS))
print(paste("RR=",RR))
print(paste("truealpha12=",truealpha12))
print(paste("truealpha13=",truealpha13))
print(paste("truealpha23=",truealpha23))
print(paste("truetau=",truetau))





#--------------------------------------------------------------------
# L Matrix, later used for multivariate Wald statistic
#--------------------------------------------------------------------
Lmat <- matrix(c(1,0,0,1), nrow=2, ncol=2)

# S=500

# -----------------------------------------------------------------
# function to create dummy variables for outcome and predictors
# -----------------------------------------------------------------



make.dummy <- function(data, resp.ind, nresponse, xvars.expand=NULL, xvars.common=NULL){
  
  data.dummy <- eval(data)
  resp.dummy <- dummy_cols(data.dummy$resp.ind)
  names(resp.dummy) <- c("resp.ind", paste0("resp.ind", 1:nresponse))
  resp.dummy1 <- resp.dummy[,-c(1:2)]
  
  if(is.null(xvars.expand)){
    
    commonxvars <- as.data.frame(data.dummy[,xvars.common])
    names(commonxvars) <- xvars.common
    cbind(resp.dummy, commonxvars)
    
  }else{
    
    nexpand <- length(xvars.expand)
    
    xvars.dummy.list <- list()
    for (i in 1:nexpand){
      xvars.dummy <- resp.dummy1 * data.dummy[,xvars.expand[i]]
      xvars.dummy <- cbind(data.dummy[,xvars.expand[i]], xvars.dummy)
      names(xvars.dummy) <- c(xvars.expand[i], paste0(xvars.expand[i], 2:nresponse))
      xvars.dummy.list[[i]] <- xvars.dummy
    }
    
    xvars.dummy.frame <- do.call("cbind", xvars.dummy.list)
    
    if (is.null(xvars.common)){
      cbind(resp.dummy, xvars.dummy.frame)
    }else{
      xvars.common.frame <- data.frame(data.dummy[,xvars.common])
      names(xvars.common.frame) <- xvars.common
      cbind(resp.dummy, xvars.dummy.frame, xvars.common.frame)
    }
  }
  
}



#------------------------------------------------------------
# function to estimate stage 1 alpha for exch structure
# for fixed number of categories
#------------------------------------------------------------


estalpha1_exch_fixed <- function(mdat, id, Z, nresponse){
  
  match.call()
  
  K <- nresponse
  G1 <- G2 <- 0
  
  for(i in id){
    csizei <- nlevels(as.factor(mdat[mdat$cluster == i,]$unit))
    wi <- 1/csizei
    cveci <- unique(mdat[mdat$cluster == i,]$unit)
    Z_i <- Z[mdat$cluster == i]
    matZ_i <- matrix(Z_i, nrow = K)
    
    ZZ1 <- 0
    for(j in 1:csizei){
      ZZ1 <- ZZ1 + t(matZ_i[,j]) %*% matZ_i[,j]
    }
    ZZ2 <- sum ( ( colSums(matZ_i) ) ^ 2 )
    
    G1 <- G1 + wi * ZZ1
    G2 <- G2 + wi * ZZ2
  }
  
  thenum <- - ( K - 1) * G1 + sqrt( ( K - 1 ) ^ 2 * G1 ^ 2 - ( G1 - G2 ) * ( ( K - 1 ) * ( G1 * ( K - 1 ) - G2 ) ) )
  theden <- ( K - 1 ) * ( ( K - 1 ) * G1 - G2 )
  
  alpha0 <- thenum / theden
  
  alpha <- ( alpha0 * ( alpha0 * ( K - 2 ) + 2 ) ) / ( 1 + alpha0 ^ 2 * ( K - 1 ) ) 
  
  R <- toeplitz(c(1, rep(alpha, K - 1)))
  
}



#------------------------------------------------------------
# function to estimate stage 1 alpha for exch structure
# for fixed number of categories (unweighted)
#------------------------------------------------------------


estalpha1_exch_fixed_unw <- function(mdat, id, Z, nresponse){
  
  match.call()
  
  K <- nresponse
  G1 <- G2 <- 0
  
  for(i in id){
    csizei <- nlevels(as.factor(mdat[mdat$cluster == i,]$unit))
    wi <- 1/csizei
    cveci <- unique(mdat[mdat$cluster == i,]$unit)
    Z_i <- Z[mdat$cluster == i]
    matZ_i <- matrix(Z_i, nrow = K)
    
    ZZ1 <- 0
    for(j in 1:csizei){
      ZZ1 <- ZZ1 + t(matZ_i[,j]) %*% matZ_i[,j]
    }
    ZZ2 <- sum ( ( colSums(matZ_i) ) ^ 2 )
    
    G1 <- G1 + ZZ1
    G2 <- G2 + ZZ2
  }
  
  thenum <- - ( K - 1) * G1 + sqrt( ( K - 1 ) ^ 2 * G1 ^ 2 - ( G1 - G2 ) * ( ( K - 1 ) * ( G1 * ( K - 1 ) - G2 ) ) )
  theden <- ( K - 1 ) * ( ( K - 1 ) * G1 - G2 )
  
  alpha0 <- thenum / theden
  
  alpha <- ( alpha0 * ( alpha0 * ( K - 2 ) + 2 ) ) / ( 1 + alpha0 ^ 2 * ( K - 1 ) ) 
  
  R <- toeplitz(c(1, rep(alpha, K - 1)))
  
}


#----------------------------------------------------------------------------
# function to estimate stage 1 alpha for unstructured correlation structure
# for fixed number of categories
#----------------------------------------------------------------------------


estalpha1_unstr_fixed <- function(mdat, id, Z, nresponse){
  
  Zmat <- 0
  
  for(i in id){
    
    csizei <- nlevels(as.factor(mdat[mdat$cluster == i,]$unit))
    wi <- 1/csizei
    cveci <- unique(mdat[mdat$cluster == i,]$unit)
    
    Zmati <- 0
    
    for (j in cveci){
      
      Z_ij <- Z[mdat$cluster == i & mdat$unit == j]
      Zmati <- Zmati + Z_ij %*% t(Z_ij)
      
    }
    
    Zmat <- Zmat + wi * Zmati
    
  }
  
  delta_k <- delta0 <- diag(rep(1, dim(Zmat)[1]))
  #iter <- 0
  diffR <- 1
  
  while (sum(abs(diffR)) > 0.000001){
    
    delta <- diag( delta_k ^ (1/2) %*% Zmat %*% delta_k ^ (1/2) ) ^ (1/2)
    delta_k1 <- diag( delta )
    diffR <- delta_k1 - delta_k
    delta_k <- delta_k1 
    
    #iter <- iter + 1
    #print(iter)
    
  }
  
  R1 <- solve( delta_k ^ (1/2) ) %*% ( delta_k ^ (1/2) %*% Zmat %*% delta_k ^ (1/2) ) ^ (1/2) %*% solve( delta_k ^ (1/2) )
  
  e <- rep(1, length = dim(R1)[1])
  v_m <- solve( hadamard.prod(R1, R1) ) %*% e
  
  diagv_m <- matrix(0, dim(R1)[1], dim(R1)[1])
  diag(diagv_m) <- v_m
  
  #R <- R1 %*% diagv_m %*% R1
  
  Zmathat <- Zmat/length(id)
  ZZ <- matrix(0, dim(Zmat)[1], dim(Zmat)[1])
  diag(ZZ) <- diag(Zmathat)
  ZZinv <- solve(ZZ)
  R <- (ZZinv) ^ (1/2) %*% Zmathat %*% (ZZinv) ^ (1/2)
  
}


#------------------------------------------------------------
# function for iterative CWGEE
#------------------------------------------------------------



mvoCWGEEgen <- function(mdat, beta0, nbeta, id, Y_var, nresponse, X_mat, corr.str, Z){
  
  match.call()
  beta <- beta0
  bdiff <- 1
  iter <- 0
  
  while(sum(abs(bdiff))>.00000001){
    DWZ <- matrix(0, nrow = nbeta)
    DWD <- matrix(0, ncol = nbeta, nrow = nbeta)
    DWZZWD <- matrix(0, ncol = nbeta, nrow = nbeta)
    
    if (corr.str == "ind")
      R <- diag(nresponse)
    if (corr.str == "unstr")
      R <- estalpha1_unstr_fixed(mdat, id, Z, nresponse)
    if (corr.str == "exch")
      R <- estalpha1_exch_fixed(mdat, id, Z, nresponse)
    
    for (i in id){
      
      csizei <- nlevels(as.factor(mdat[mdat$cluster == i,]$unit))
      wi <- 1/csizei
      unitveci <- unique(mdat[mdat$cluster == i,]$unit)
      
      DWZj <- matrix(0, nrow = nbeta)
      DWDj <- matrix(0, ncol = nbeta, nrow = nbeta)
      
      for(j in unitveci){
        
        y <- as.matrix( Y_var[mdat$cluster==i & mdat$unit == j] )
        x <- as.matrix( X_mat[mdat$cluster==i & mdat$unit == j,] )
        u <- exp( x %*% beta ) / ( 1 + exp( x %*% beta ) ) 
        dudb <- exp( x %*% beta ) / ( 1 + exp( x %*% beta ) ) ^ 2
        D <- x[,1] * dudb
        for(p in 2:nbeta) D <- cbind(D, x[,p] * dudb)
        vv <- u * ( 1 - u )
        V <- matrix(ncol = length(vv), nrow = length(vv), 0)
        diag(V) <- vv
        W <- V ^ (1/2) %*% R %*% V ^ (1/2)
        invW <- solve(W)
        DWZj <- DWZj + t(D) %*% invW %*% ( y - u ) 
        DWDj <- DWDj + t(D) %*% invW %*% D
        
      }
      
      DWZ <- DWZ + wi * DWZj
      DWD <- DWD + wi * DWDj
      DWZZWD <- DWZZWD + ( wi * DWZj ) %*% t( wi * DWZj )
      
    }
    
    invDWD <- solve(DWD)
    bdiff <- invDWD %*% DWZ
    beta <- beta + bdiff
    covbeta <- invDWD %*% DWZZWD %*% invDWD
    U <- exp(X_mat %*% beta) / (1 + exp(X_mat %*% beta))
    var <- U * (1 - U)
    Z <- (Y_var - U) / sqrt(var)
    iter <- iter + 1
  }
  
  fit <- list()
  fit$coefficients <- beta
  fit$robust.variance <- covbeta
  fit$robust.se <- sqrt(diag(covbeta))
  fit$wald.chisq <- (beta/sqrt(diag(covbeta)))^2
  fit$p.value <- 1 - pchisq(fit$wald.chisq, df = 1)
  fit$R <- R
  fit$niter <- iter
  fit
  
}

#------------------------------------------------------------
# function for iterative GEE
#------------------------------------------------------------



mvoGEEgen <- function(mdat, beta0, nbeta, id, Y_var, nresponse, X_mat, corr.str, Z){
  
  match.call()
  beta <- beta0
  bdiff <- 1
  iter <- 0
  
  while(sum(abs(bdiff))>.00000001){
    DWZ <- matrix(0, nrow = nbeta)
    DWD <- matrix(0, ncol = nbeta, nrow = nbeta)
    DWZZWD <- matrix(0, ncol = nbeta, nrow = nbeta)
    
    if (corr.str == "ind")
      R <- diag(nresponse)
    if (corr.str == "unstr")
      R <- estalpha1_unstr_fixed(mdat, id, Z, nresponse)
    if (corr.str == "exch")
      R <- estalpha1_exch_fixed_unw(mdat, id, Z, nresponse)
    
    for (i in id){
      
      csizei <- nlevels(as.factor(mdat[mdat$cluster == i,]$unit))
      wi <- 1/csizei
      unitveci <- unique(mdat[mdat$cluster == i,]$unit)
      
      DWZj <- matrix(0, nrow = nbeta)
      DWDj <- matrix(0, ncol = nbeta, nrow = nbeta)
      
      for(j in unitveci){
        
        y <- as.matrix( Y_var[mdat$cluster==i & mdat$unit == j] )
        x <- as.matrix( X_mat[mdat$cluster==i & mdat$unit == j,] )
        u <- exp( x %*% beta ) / ( 1 + exp( x %*% beta ) ) 
        dudb <- exp( x %*% beta ) / ( 1 + exp( x %*% beta ) ) ^ 2
        D <- x[,1] * dudb
        for(p in 2:nbeta) D <- cbind(D, x[,p] * dudb)
        vv <- u * ( 1 - u )
        V <- matrix(ncol = length(vv), nrow = length(vv), 0)
        diag(V) <- vv
        W <- V ^ (1/2) %*% R %*% V ^ (1/2)
        invW <- solve(W)
        DWZj <- DWZj + t(D) %*% invW %*% ( y - u ) 
        DWDj <- DWDj + t(D) %*% invW %*% D
        
      }
      
      DWZ <- DWZ + DWZj
      DWD <- DWD + DWDj
      DWZZWD <- DWZZWD + ( DWZj ) %*% t( DWZj )
      
    }
    
    invDWD <- solve(DWD)
    bdiff <- invDWD %*% DWZ
    beta <- beta + bdiff
    covbeta <- invDWD %*% DWZZWD %*% invDWD
    U <- exp(X_mat %*% beta) / (1 + exp(X_mat %*% beta))
    var <- U * (1 - U)
    Z <- (Y_var - U) / sqrt(var)
    iter <- iter + 1
  }
  
  fit <- list()
  fit$coefficients <- beta
  fit$robust.variance <- covbeta
  fit$robust.se <- sqrt(diag(covbeta))
  fit$wald.chisq <- (beta/sqrt(diag(covbeta)))^2
  fit$p.value <- 1 - pchisq(fit$wald.chisq, df = 1)
  fit$R <- R
  fit$niter <- iter
  fit
  
}


#------------------------------------------------------------
# main function for multivariate CWGEE
#------------------------------------------------------------



mvoCWGEE <- function(formula, data, cluster, resp.ind, unit, corr.str, common.slope = NULL){
  
  call <- match.call()
  mcall <- match.call(expand.dots = FALSE)
  mf <- match(c("formula", "data", "cluster", "resp.ind", "unit"), names(mcall), 0L)
  m <- mcall[c(1L, mf)]
  if (is.null(m$cluster))
    m$cluster <- as.name("id")
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  y <- model.response(m)
  
  corr.strs <- c("ind", "exch", "unstr")
  strcheck <- match(corr.str, corr.strs, -1)
  if (strcheck < 1){
    stop("unknown working correlation structure")
  }
  
  cluster <- model.extract(m, "cluster")
  resp.ind <- model.extract(m, "resp.ind")
  unit <- model.extract(m, "unit")
  
  nresponse <- max(resp.ind)
  
  mterms <- attr(m, "terms")
  xvars <- as.character(attr(mterms, "variables"))[-c(1:2)]
  names(m) <- c("y", xvars, "cluster", "resp.ind", "unit")
  
  if (identical(xvars, common.slope)){
    xvars.expand <- NULL
    xvars.common <- xvars
  }else if (is.null(common.slope)){
    xvars.expand <- xvars
    xvars.common <- NULL
  }else{
    xvars.expand <- setdiff(xvars, common.slope)
    xvars.common <- intersect(xvars, common.slope)
  }
  
  dummyout <- make.dummy(data=m, resp.ind, nresponse, xvars.expand, xvars.common)[,-1]
  newxvars <- names(dummyout)
  
  mdat <- as.data.frame(cbind(y, dummyout, resp.ind, cluster, unit))
  mformula <- as.formula(paste("y", "~", paste(c("-1", paste(newxvars, collapse="+")), collapse="+")))
  nbeta <- length(newxvars)
  X_mat <- model.matrix(mformula, data = mdat)
  mframe <- model.frame(mformula, data = mdat)
  Y_var <- model.response(mframe)
  id <- unique.cluster <- unique(cluster)
  
  # use GLM to get initial beta estimates
  beta0 <- as.vector(coef(glm(mformula, data = mdat, family = "binomial")))
  Z <- glm(mformula, data = mdat, family = "binomial")$residuals
  #list(Y, nresponse, m, common.slope, unit, cluster, xvars, xvars.expand, xvars.common)
  
  cwgeeout <- mvoCWGEEgen(mdat, beta0, nbeta, id, Y_var, nresponse, X_mat, corr.str, Z)
  #list(mdat, beta0, nbeta, id, Y_var, nresponse, X_mat, corr.str, Z)
  
  coef.names <- newxvars
  
  result <- list()
  result$call <- match.call()
  result$coefficients <- cwgeeout$coefficients
  result$coef.names <- coef.names
  result$robust.variance <- cwgeeout$robust.variance
  result$robust.se <- cwgeeout$robust.se
  result$wald.chisq <- cwgeeout$wald.chisq
  result$p.value <- cwgeeout$p.value
  result$corr.str <- corr.str
  result$corr.matrix <- cwgeeout$R
  result$niter <- cwgeeout$niter
  class(result) <- "cwgee"
  result
  
}


#------------------------------------------------------------
# main function for multivariate GEE
#------------------------------------------------------------


mvoGEE <- function(formula, data, cluster, resp.ind, unit, corr.str, common.slope = NULL){
  
  call <- match.call()
  mcall <- match.call(expand.dots = FALSE)
  mf <- match(c("formula", "data", "cluster", "resp.ind", "unit"), names(mcall), 0L)
  m <- mcall[c(1L, mf)]
  if (is.null(m$cluster))
    m$cluster <- as.name("id")
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  y <- model.response(m)
  
  corr.strs <- c("ind", "exch", "unstr")
  strcheck <- match(corr.str, corr.strs, -1)
  if (strcheck < 1){
    stop("unknown working correlation structure")
  }
  
  cluster <- model.extract(m, "cluster")
  resp.ind <- model.extract(m, "resp.ind")
  unit <- model.extract(m, "unit")
  
  nresponse <- max(resp.ind)
  
  mterms <- attr(m, "terms")
  xvars <- as.character(attr(mterms, "variables"))[-c(1:2)]
  names(m) <- c("y", xvars, "cluster", "resp.ind", "unit")
  
  if (identical(xvars, common.slope)){
    xvars.expand <- NULL
    xvars.common <- xvars
  }else if (is.null(common.slope)){
    xvars.expand <- xvars
    xvars.common <- NULL
  }else{
    xvars.expand <- setdiff(xvars, common.slope)
    xvars.common <- intersect(xvars, common.slope)
  }
  
  dummyout <- make.dummy(data=m, resp.ind, nresponse, xvars.expand, xvars.common)[,-1]
  newxvars <- names(dummyout)
  
  mdat <- as.data.frame(cbind(y, dummyout, resp.ind, cluster, unit))
  mformula <- as.formula(paste("y", "~", paste(c("-1", paste(newxvars, collapse="+")), collapse="+")))
  nbeta <- length(newxvars)
  X_mat <- model.matrix(mformula, data = mdat)
  mframe <- model.frame(mformula, data = mdat)
  Y_var <- model.response(mframe)
  id <- unique.cluster <- unique(cluster)
  
  # use GLM to get initial beta estimates
  beta0 <- as.vector(coef(glm(mformula, data = mdat, family = "binomial")))
  Z <- glm(mformula, data = mdat, family = "binomial")$residuals
  #list(Y, nresponse, m, common.slope, unit, cluster, xvars, xvars.expand, xvars.common)
  
  geeout <- mvoGEEgen(mdat, beta0, nbeta, id, Y_var, nresponse, X_mat, corr.str, Z)
  #list(mdat, beta0, nbeta, id, Y_var, nresponse, X_mat, corr.str, Z)
  
  coef.names <- newxvars
  
  result <- list()
  result$call <- match.call()
  result$coefficients <- geeout$coefficients
  result$coef.names <- coef.names
  result$robust.variance <- geeout$robust.variance
  result$robust.se <- geeout$robust.se
  result$wald.chisq <- geeout$wald.chisq
  result$p.value <- geeout$p.value
  result$corr.str <- corr.str
  result$corr.matrix <- geeout$R
  result$niter <- geeout$niter
  class(result) <- "cwgee"
  result
  
}


#------------------------------------------------------------
# for package: print output from mvoCWGEE
#------------------------------------------------------------


print.cwgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

#------------------------------------------------------------
# for package: print summar from mvoCWGEE
#------------------------------------------------------------


print.summary.cwgee <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients)
  cat("\nWorking correlation structure:\n")
  print(x$corr.str)
  cat("\nCorrelation matrix estimate:\n")
  print(x$corr.matrix)
  cat("\nTemporal working correlation structure:\n")
  print(x$time.str)
  cat("\nCorrelation estimate:\n")
  print(x$corr.coef)
  cat("\nNumber of iterations:\n")
  print(x$niter)
}

#------------------------------------------------------------
# for package: summary of mvoCWGEE
#------------------------------------------------------------


summary.cwgee <- function(object, ...) {
  robust.se <- sqrt(diag(object$robust.variance))
  robust.z <- coef(object)/robust.se
  pvalue <- 2 * (1 - pnorm(abs(robust.z)))
  TAB <- cbind(Estimate = coef(object), 
               Robust_S.E. = robust.se, 
               Robust_z = robust.z, 
               p.value = pvalue)
  TAB <- round(TAB, 5)
  colnames(TAB) <- c("Estimate", "san.se", "san.z", "Pr(>|san.z|)")
  rownames(TAB) <- object$coef.names
  res <- list(coefficients = TAB, 
              call = object$call,
              corr.str = object$corr.str,
              corr.matrix = object$corr.matrix,
              corr.coef = object$alpha,
              time.str = object$time.str,
              niter = object$niter)
  class(res) <- "summary.cwgee"
  res
}

#---------------------------------------------
# set seed for simulation
#---------------------------------------------

theseed <- SS * RR - SS + RR
set.seed(theseed)


truecatvec <- c(truea1, truea2, truea3)
numoutcome <- length(truecatvec)
numcat <- length(truecatvec) + 1
catvec <- c(1:numcat)

#-------------------------------------------------------------------------------
# simulate correlated binary outcomes with informative cluster size 
#-------------------------------------------------------------------------------

rho12_i <- sin( pi * truealpha12 / 2)
rho13_i <- sin( pi * truealpha13 / 2)
rho23_i <- sin( pi * truealpha23 / 2)
tau_i <- sin( pi * truetau / 2)
phi <- 0.5

therho12 <- rho12_i
therho13 <- rho13_i
therho23 <- rho23_i
thetau <- tau_i

## correlation between teeth
thesigma1 <- toeplitz(c(1, rep(tau_i, nteeth - 1)))
## correlation between outcomes
thesigma2 <- matrix(c(1, rho12_i, rho13_i,
                      rho12_i, 1, rho23_i,
                      rho13_i, rho23_i, 1), ncol = numoutcome)
## use kronecker product to join two correlations
fullsigma <- nearPD( thesigma1 %x% thesigma2 )$mat
themu <- rep(0, nrow(fullsigma))

## generate binary predictor 
x <- c(rep(0, N / 2), rep(1, N / 2))

ally <- vector("list", N)
nteethpp <- baserisk <- rep(NA, N)
teethpp <- vector("list", N)

nstay <- 0
maxnstayiter <- 1000


## simulate data for each patient i

for (i in 1:N){
  
  nstayiter <- maxnstayiter
  
  while(nstayiter == maxnstayiter){
    
    Z <- mvrnorm(n = 1, mu = themu, Sigma = fullsigma)
    u <- pnorm(Z)
    b <- ( 1 / phi ) * log( sin( phi * pi * u) / sin( phi * pi * ( 1 - u ) ) )
    
    ## generate baseline risk
    lambda <- exp(mean(-b)) / ( 1 + exp(mean(-b)) )
    baserisk[i] <- mean(-b)
    
    nstayiter <- 0
    
    ## generate cluster size 
    while(nstay < 2 & nstayiter < maxnstayiter){
      nstay <- rbinom(n = 1, size = nteeth, prob = lambda)
      nstayiter <- nstayiter + 1
    }
  }
  
  ## adjust vector of teeth to match cluster size
  nteethpp[i] <- nstay
  tstay <- sort(sample(c(1:nteeth), nstay))
  teethpp[[i]] <- tstay
  
  tstayvec <- c( tstay, tstay + nteeth, tstay + 2 * nteeth)
  tstayvec <- sort(tstayvec)
  
  bstay <- b[tstayvec]
  
  yvec <- vector("list", nstay)
  
  ## generate outcome
  for (t in 1:nstay){
    
    tt <- t + nstay
    ttt <- tt + nstay
    
    truexb <- truebeta1 * x[i] 
    
    p1 <- exp( bstay[t] + ( 1 / phi ) * ( truea1 + truexb ) ) / ( 1 + exp( bstay[t] + ( 1 / phi ) * ( truea1 + truexb ) ) )
    p2 <- exp( bstay[tt] + ( 1 / phi ) * ( truea2 + truexb ) ) / ( 1 + exp( bstay[tt] + ( 1 / phi ) * ( truea2 + truexb ) ) )
    p3 <- exp( bstay[ttt] + ( 1 / phi ) * ( truea3 + truexb ) ) / ( 1 + exp( bstay[ttt] + ( 1 / phi ) * ( truea3 + truexb ) ) )
    
    y1dummy <- rbinom( n = 1, size = 1, prob = p1 )
    y2dummy <- rbinom( n = 1, size = 1, prob = p2 )
    y3dummy <- rbinom( n = 1, size = 1, prob = p3 )
    
    yvec[[t]] <- c( y1dummy, y2dummy, y3dummy )
    
  }
  
  ally[[i]] <- unlist(yvec)
  
  nstay <- 0
}

## create data

y <- unlist(ally)
x <- rep(x, times = nteethpp * numoutcome)

subject <- rep(1:N, times = nteethpp * numoutcome)

outcomelist <- vector("list", N)
for(i in 1:N) outcomelist[[i]] <- rep(1:numoutcome, nteethpp[i])
outcome <- unlist(outcomelist)

toothlist <- vector("list", N)
for(i in 1:N) toothlist[[i]] <- rep(1:nteethpp[i], each = numoutcome)
tooth <- unlist(toothlist)

simdat <- as.data.frame(cbind(subject, tooth, outcome, x, y))

out1y <- simdat[which(simdat$outcome == 1),]$y
out2y <- simdat[which(simdat$outcome == 2),]$y
out3y <- simdat[which(simdat$outcome == 3),]$y
prop.table(table(out1y, out2y, out3y))

csize <- as.data.frame(aggregate(simdat$tooth, by = list(simdat$subject), max))
names(csize) <- c("subject", "csize")
simdat_wcsize <- merge(simdat, csize, by = "subject")


#--------------------------------------------------------------------------------------
# anlalyze one outcome at a time using GEE and CWGEE
#--------------------------------------------------------------------------------------

## outcome1
simdat1 <- subset(simdat_wcsize, simdat_wcsize$outcome == 1)

outGEE1 <- geeglm(y ~ x, data = simdat1, id = subject)
outGEE1_resp1est <- summary(outGEE1)$coef[1,1]
outGEE1_betaest <- summary(outGEE1)$coef[2,1]
outGEE1_resp1se <- summary(outGEE1)$coef[1,2]
outGEE1_betase <- summary(outGEE1)$coef[2,2]
outGEE1_betareject <- ifelse(summary(outGEE1)$coef[2,4] < 0.05, 1, 0)

outCWGEE1 <- geeglm(y ~ x, data = simdat1, id = subject, weights = csize)
outCWGEE1_resp1est <- summary(outCWGEE1)$coef[1,1]
outCWGEE1_betaest <- summary(outCWGEE1)$coef[2,1]
outCWGEE1_resp1se <- summary(outCWGEE1)$coef[1,2]
outCWGEE1_betase <- summary(outCWGEE1)$coef[2,2]
outCWGEE1_betareject <- ifelse(summary(outCWGEE1)$coef[2,4] < 0.05, 1, 0)


## outcome2
simdat2 <- subset(simdat_wcsize, simdat_wcsize$outcome == 2)

outGEE2 <- geeglm(y ~ x, data = simdat2, id = subject)
outGEE2_resp1est <- summary(outGEE2)$coef[1,1]
outGEE2_betaest <- summary(outGEE2)$coef[2,1]
outGEE2_resp1se <- summary(outGEE2)$coef[1,2]
outGEE2_betase <- summary(outGEE2)$coef[2,2]
outGEE2_betareject <- ifelse(summary(outGEE2)$coef[2,4] < 0.05, 1, 0)

outCWGEE2 <- geeglm(y ~ x, data = simdat2, id = subject, weights = csize)
outCWGEE2_resp1est <- summary(outCWGEE2)$coef[1,1]
outCWGEE2_betaest <- summary(outCWGEE2)$coef[2,1]
outCWGEE2_resp1se <- summary(outCWGEE2)$coef[1,2]
outCWGEE2_betase <- summary(outCWGEE2)$coef[2,2]
outCWGEE2_betareject <- ifelse(summary(outCWGEE2)$coef[2,4] < 0.05, 1, 0)


## outcome3
simdat3 <- subset(simdat_wcsize, simdat_wcsize$outcome == 3)

outGEE3 <- geeglm(y ~ x, data = simdat3, id = subject)
outGEE3_resp1est <- summary(outGEE3)$coef[1,1]
outGEE3_betaest <- summary(outGEE3)$coef[2,1]
outGEE3_resp1se <- summary(outGEE3)$coef[1,2]
outGEE3_betase <- summary(outGEE3)$coef[2,2]
outGEE3_betareject <- ifelse(summary(outGEE3)$coef[2,4] < 0.05, 1, 0)

outCWGEE3 <- geeglm(y ~ x, data = simdat3, id = subject, weights = csize)
outCWGEE3_resp1est <- summary(outCWGEE3)$coef[1,1]
outCWGEE3_betaest <- summary(outCWGEE3)$coef[2,1]
outCWGEE3_resp1se <- summary(outCWGEE3)$coef[1,2]
outCWGEE3_betase <- summary(outCWGEE3)$coef[2,2]
outCWGEE3_betareject <- ifelse(summary(outCWGEE3)$coef[2,4] < 0.05, 1, 0)


#-------------------------------------------------------------------
# analyze all outcomes toghether using GEE and CWGEE
#-------------------------------------------------------------------


## unstructured correlation with common slope for x

outGEEunstr_commonslope <- mvoGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "unstr", common.slope = c("x"))
outGEEunstr_commonslope_resp1est <- summary(outGEEunstr_commonslope)$coeff[1,1]
outGEEunstr_commonslope_resp2est <- summary(outGEEunstr_commonslope)$coeff[2,1]
outGEEunstr_commonslope_resp3est <- summary(outGEEunstr_commonslope)$coeff[3,1]
outGEEunstr_commonslope_betaest <- summary(outGEEunstr_commonslope)$coeff[4,1]
outGEEunstr_commonslope_resp1se <- summary(outGEEunstr_commonslope)$coeff[1,2]
outGEEunstr_commonslope_resp2se <- summary(outGEEunstr_commonslope)$coeff[2,2]
outGEEunstr_commonslope_resp3se <- summary(outGEEunstr_commonslope)$coeff[3,2]
outGEEunstr_commonslope_betase <- summary(outGEEunstr_commonslope)$coeff[4,2]
outGEEunstr_commonslope_betareject <- ifelse(summary(outGEEunstr_commonslope)$coeff[4,4] < 0.05, 1, 0)
outGEEunstr_commonslope_alpha12est <- summary(outGEEunstr_commonslope)$corr.matrix[1,2]
outGEEunstr_commonslope_alpha13est <- summary(outGEEunstr_commonslope)$corr.matrix[1,3]
outGEEunstr_commonslope_alpha23est <- summary(outGEEunstr_commonslope)$corr.matrix[2,3]


outCWGEEunstr_commonslope <- mvoCWGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "unstr", common.slope = c("x"))
outCWGEEunstr_commonslope_resp1est <- summary(outCWGEEunstr_commonslope)$coeff[1,1]
outCWGEEunstr_commonslope_resp2est <- summary(outCWGEEunstr_commonslope)$coeff[2,1]
outCWGEEunstr_commonslope_resp3est <- summary(outCWGEEunstr_commonslope)$coeff[3,1]
outCWGEEunstr_commonslope_betaest <- summary(outCWGEEunstr_commonslope)$coeff[4,1]
outCWGEEunstr_commonslope_resp1se <- summary(outCWGEEunstr_commonslope)$coeff[1,2]
outCWGEEunstr_commonslope_resp2se <- summary(outCWGEEunstr_commonslope)$coeff[2,2]
outCWGEEunstr_commonslope_resp3se <- summary(outCWGEEunstr_commonslope)$coeff[3,2]
outCWGEEunstr_commonslope_betase <- summary(outCWGEEunstr_commonslope)$coeff[4,2]
outCWGEEunstr_commonslope_betareject <- ifelse(summary(outCWGEEunstr_commonslope)$coeff[4,4] < 0.05, 1, 0)
outCWGEEunstr_commonslope_alpha12est <- summary(outCWGEEunstr_commonslope)$corr.matrix[1,2]
outCWGEEunstr_commonslope_alpha13est <- summary(outCWGEEunstr_commonslope)$corr.matrix[1,3]
outCWGEEunstr_commonslope_alpha23est <- summary(outCWGEEunstr_commonslope)$corr.matrix[2,3]


## exchangeable correlation with common slope for x

outGEEexch_commonslope <- mvoGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "exch", common.slope = c("x"))
outGEEexch_commonslope_resp1est <- summary(outGEEexch_commonslope)$coeff[1,1]
outGEEexch_commonslope_resp2est <- summary(outGEEexch_commonslope)$coeff[2,1]
outGEEexch_commonslope_resp3est <- summary(outGEEexch_commonslope)$coeff[3,1]
outGEEexch_commonslope_betaest <- summary(outGEEexch_commonslope)$coeff[4,1]
outGEEexch_commonslope_resp1se <- summary(outGEEexch_commonslope)$coeff[1,2]
outGEEexch_commonslope_resp2se <- summary(outGEEexch_commonslope)$coeff[2,2]
outGEEexch_commonslope_resp3se <- summary(outGEEexch_commonslope)$coeff[3,2]
outGEEexch_commonslope_betase <- summary(outGEEexch_commonslope)$coeff[4,2]
outGEEexch_commonslope_betareject <- ifelse(summary(outGEEexch_commonslope)$coeff[4,4] < 0.05, 1, 0)
outGEEexch_commonslope_alpha12est <- summary(outGEEexch_commonslope)$corr.matrix[1,2]
outGEEexch_commonslope_alpha13est <- summary(outGEEexch_commonslope)$corr.matrix[1,3]
outGEEexch_commonslope_alpha23est <- summary(outGEEexch_commonslope)$corr.matrix[2,3]


outCWGEEexch_commonslope <- mvoCWGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "exch", common.slope = c("x"))
outCWGEEexch_commonslope_resp1est <- summary(outCWGEEexch_commonslope)$coeff[1,1]
outCWGEEexch_commonslope_resp2est <- summary(outCWGEEexch_commonslope)$coeff[2,1]
outCWGEEexch_commonslope_resp3est <- summary(outCWGEEexch_commonslope)$coeff[3,1]
outCWGEEexch_commonslope_betaest <- summary(outCWGEEexch_commonslope)$coeff[4,1]
outCWGEEexch_commonslope_resp1se <- summary(outCWGEEexch_commonslope)$coeff[1,2]
outCWGEEexch_commonslope_resp2se <- summary(outCWGEEexch_commonslope)$coeff[2,2]
outCWGEEexch_commonslope_resp3se <- summary(outCWGEEexch_commonslope)$coeff[3,2]
outCWGEEexch_commonslope_betase <- summary(outCWGEEexch_commonslope)$coeff[4,2]
outCWGEEexch_commonslope_betareject <- ifelse(summary(outCWGEEexch_commonslope)$coeff[4,4] < 0.05, 1, 0)
outCWGEEexch_commonslope_alpha12est <- summary(outCWGEEexch_commonslope)$corr.matrix[1,2]
outCWGEEexch_commonslope_alpha13est <- summary(outCWGEEexch_commonslope)$corr.matrix[1,3]
outCWGEEexch_commonslope_alpha23est <- summary(outCWGEEexch_commonslope)$corr.matrix[2,3]


## independent correlation with common slope for x

outGEEind_commonslope <- mvoGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "ind", common.slope = c("x"))
outGEEind_commonslope_resp1est <- summary(outGEEind_commonslope)$coeff[1,1]
outGEEind_commonslope_resp2est <- summary(outGEEind_commonslope)$coeff[2,1]
outGEEind_commonslope_resp3est <- summary(outGEEind_commonslope)$coeff[3,1]
outGEEind_commonslope_betaest <- summary(outGEEind_commonslope)$coeff[4,1]
outGEEind_commonslope_resp1se <- summary(outGEEind_commonslope)$coeff[1,2]
outGEEind_commonslope_resp2se <- summary(outGEEind_commonslope)$coeff[2,2]
outGEEind_commonslope_resp3se <- summary(outGEEind_commonslope)$coeff[3,2]
outGEEind_commonslope_betase <- summary(outGEEind_commonslope)$coeff[4,2]
outGEEind_commonslope_betareject <- ifelse(summary(outGEEind_commonslope)$coeff[4,4] < 0.05, 1, 0)
outGEEind_commonslope_alpha12est <- summary(outGEEind_commonslope)$corr.matrix[1,2]
outGEEind_commonslope_alpha13est <- summary(outGEEind_commonslope)$corr.matrix[1,3]
outGEEind_commonslope_alpha23est <- summary(outGEEind_commonslope)$corr.matrix[2,3]


outCWGEEind_commonslope <- mvoCWGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "ind", common.slope = c("x"))
outCWGEEind_commonslope_resp1est <- summary(outCWGEEind_commonslope)$coeff[1,1]
outCWGEEind_commonslope_resp2est <- summary(outCWGEEind_commonslope)$coeff[2,1]
outCWGEEind_commonslope_resp3est <- summary(outCWGEEind_commonslope)$coeff[3,1]
outCWGEEind_commonslope_betaest <- summary(outCWGEEind_commonslope)$coeff[4,1]
outCWGEEind_commonslope_resp1se <- summary(outCWGEEind_commonslope)$coeff[1,2]
outCWGEEind_commonslope_resp2se <- summary(outCWGEEind_commonslope)$coeff[2,2]
outCWGEEind_commonslope_resp3se <- summary(outCWGEEind_commonslope)$coeff[3,2]
outCWGEEind_commonslope_betase <- summary(outCWGEEind_commonslope)$coeff[4,2]
outCWGEEind_commonslope_betareject <- ifelse(summary(outCWGEEind_commonslope)$coeff[4,4] < 0.05, 1, 0)
outCWGEEind_commonslope_alpha12est <- summary(outCWGEEind_commonslope)$corr.matrix[1,2]
outCWGEEind_commonslope_alpha13est <- summary(outCWGEEind_commonslope)$corr.matrix[1,3]
outCWGEEind_commonslope_alpha23est <- summary(outCWGEEind_commonslope)$corr.matrix[2,3]


## unstructured correlation with unique slopes for x

outGEEunstr <- mvoGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "unstr")
outGEEunstr_resp1est <- summary(outGEEunstr)$coeff[1,1]
outGEEunstr_resp2est <- summary(outGEEunstr)$coeff[2,1]
outGEEunstr_resp3est <- summary(outGEEunstr)$coeff[3,1]
outGEEunstr_betaest <- summary(outGEEunstr)$coeff[4,1]
outGEEunstr_beta2est <- summary(outGEEunstr)$coeff[5,1]
outGEEunstr_beta3est <- summary(outGEEunstr)$coeff[6,1]
outGEEunstr_resp1se <- summary(outGEEunstr)$coeff[1,2]
outGEEunstr_resp2se <- summary(outGEEunstr)$coeff[2,2]
outGEEunstr_resp3se <- summary(outGEEunstr)$coeff[3,2]
outGEEunstr_betase <- summary(outGEEunstr)$coeff[4,2]
outGEEunstr_beta2se <- summary(outGEEunstr)$coeff[5,2]
outGEEunstr_beta3se <- summary(outGEEunstr)$coeff[6,2]
outGEEunstr_betareject <- ifelse(summary(outGEEunstr)$coeff[4,4] < 0.05, 1, 0)
outGEEunstr_beta2reject <- ifelse(summary(outGEEunstr)$coeff[5,4] < 0.05, 1, 0)
outGEEunstr_beta3reject <- ifelse(summary(outGEEunstr)$coeff[6,4] < 0.05, 1, 0)
outGEEunstr_alpha12est <- summary(outGEEunstr)$corr.matrix[1,2]
outGEEunstr_alpha13est <- summary(outGEEunstr)$corr.matrix[1,3]
outGEEunstr_alpha23est <- summary(outGEEunstr)$corr.matrix[2,3]
### multivariate Wald statistic for H0:beta12=beta13=0
outGEEunstr_beta2beta3est <- summary(outGEEunstr)$coeff[5:6,1]
outGEEunstr_beta2beta3cov <- outGEEunstr$robust.variance[5:6,5:6]
outGEEunstr_wald <- as.numeric(t(Lmat %*% outGEEunstr_beta2beta3est) %*% solve(Lmat %*% outGEEunstr_beta2beta3cov %*% Lmat) %*% (Lmat %*% outGEEunstr_beta2beta3est))  
outGEEunstr_waldp <- 1-pchisq(outGEEunstr_wald, df=nrow(Lmat))
outGEEunstr_beta2beta3reject <- ifelse(outGEEunstr_waldp < 0.05, 1, 0)


outCWGEEunstr <- mvoCWGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "unstr")
outCWGEEunstr_resp1est <- summary(outCWGEEunstr)$coeff[1,1]
outCWGEEunstr_resp2est <- summary(outCWGEEunstr)$coeff[2,1]
outCWGEEunstr_resp3est <- summary(outCWGEEunstr)$coeff[3,1]
outCWGEEunstr_betaest <- summary(outCWGEEunstr)$coeff[4,1]
outCWGEEunstr_beta2est <- summary(outCWGEEunstr)$coeff[5,1]
outCWGEEunstr_beta3est <- summary(outCWGEEunstr)$coeff[6,1]
outCWGEEunstr_resp1se <- summary(outCWGEEunstr)$coeff[1,2]
outCWGEEunstr_resp2se <- summary(outCWGEEunstr)$coeff[2,2]
outCWGEEunstr_resp3se <- summary(outCWGEEunstr)$coeff[3,2]
outCWGEEunstr_betase <- summary(outCWGEEunstr)$coeff[4,2]
outCWGEEunstr_beta2se <- summary(outCWGEEunstr)$coeff[5,2]
outCWGEEunstr_beta3se <- summary(outCWGEEunstr)$coeff[6,2]
outCWGEEunstr_betareject <- ifelse(summary(outCWGEEunstr)$coeff[4,4] < 0.05, 1, 0)
outCWGEEunstr_beta2reject <- ifelse(summary(outCWGEEunstr)$coeff[5,4] < 0.05, 1, 0)
outCWGEEunstr_beta3reject <- ifelse(summary(outCWGEEunstr)$coeff[6,4] < 0.05, 1, 0)
outCWGEEunstr_alpha12est <- summary(outCWGEEunstr)$corr.matrix[1,2]
outCWGEEunstr_alpha13est <- summary(outCWGEEunstr)$corr.matrix[1,3]
outCWGEEunstr_alpha23est <- summary(outCWGEEunstr)$corr.matrix[2,3]
### multivariate Wald statistic for H0:beta12=beta13=0
outCWGEEunstr_beta2beta3est <- summary(outCWGEEunstr)$coeff[5:6,1]
outCWGEEunstr_beta2beta3cov <- outCWGEEunstr$robust.variance[5:6,5:6]
outCWGEEunstr_wald <- as.numeric(t(Lmat %*% outCWGEEunstr_beta2beta3est) %*% solve(Lmat %*% outCWGEEunstr_beta2beta3cov %*% Lmat) %*% (Lmat %*% outCWGEEunstr_beta2beta3est))  
outCWGEEunstr_waldp <- 1-pchisq(outCWGEEunstr_wald, df=nrow(Lmat))
outCWGEEunstr_beta2beta3reject <- ifelse(outCWGEEunstr_waldp < 0.05, 1, 0)


## exchangeable correlation with unique slopes for x

outGEEexch <- mvoGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "exch")
outGEEexch_resp1est <- summary(outGEEexch)$coeff[1,1]
outGEEexch_resp2est <- summary(outGEEexch)$coeff[2,1]
outGEEexch_resp3est <- summary(outGEEexch)$coeff[3,1]
outGEEexch_betaest <- summary(outGEEexch)$coeff[4,1]
outGEEexch_beta2est <- summary(outGEEexch)$coeff[5,1]
outGEEexch_beta3est <- summary(outGEEexch)$coeff[6,1]
outGEEexch_resp1se <- summary(outGEEexch)$coeff[1,2]
outGEEexch_resp2se <- summary(outGEEexch)$coeff[2,2]
outGEEexch_resp3se <- summary(outGEEexch)$coeff[3,2]
outGEEexch_betase <- summary(outGEEexch)$coeff[4,2]
outGEEexch_beta2se <- summary(outGEEexch)$coeff[5,2]
outGEEexch_beta3se <- summary(outGEEexch)$coeff[6,2]
outGEEexch_betareject <- ifelse(summary(outGEEexch)$coeff[4,4] < 0.05, 1, 0)
outGEEexch_beta2reject <- ifelse(summary(outGEEexch)$coeff[5,4] < 0.05, 1, 0)
outGEEexch_beta3reject <- ifelse(summary(outGEEexch)$coeff[6,4] < 0.05, 1, 0)
outGEEexch_alpha12est <- summary(outGEEexch)$corr.matrix[1,2]
outGEEexch_alpha13est <- summary(outGEEexch)$corr.matrix[1,3]
outGEEexch_alpha23est <- summary(outGEEexch)$corr.matrix[2,3]
### multivariate Wald statistic for H0:beta12=beta13=0
outGEEexch_beta2beta3est <- summary(outGEEexch)$coeff[5:6,1]
outGEEexch_beta2beta3cov <- outGEEexch$robust.variance[5:6,5:6]
outGEEexch_wald <- as.numeric(t(Lmat %*% outGEEexch_beta2beta3est) %*% solve(Lmat %*% outGEEexch_beta2beta3cov %*% Lmat) %*% (Lmat %*% outGEEexch_beta2beta3est))  
outGEEexch_waldp <- 1-pchisq(outGEEexch_wald, df=nrow(Lmat))
outGEEexch_beta2beta3reject <- ifelse(outGEEexch_waldp < 0.05, 1, 0)

outCWGEEexch <- mvoCWGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "exch")
outCWGEEexch_resp1est <- summary(outCWGEEexch)$coeff[1,1]
outCWGEEexch_resp2est <- summary(outCWGEEexch)$coeff[2,1]
outCWGEEexch_resp3est <- summary(outCWGEEexch)$coeff[3,1]
outCWGEEexch_betaest <- summary(outCWGEEexch)$coeff[4,1]
outCWGEEexch_beta2est <- summary(outCWGEEexch)$coeff[5,1]
outCWGEEexch_beta3est <- summary(outCWGEEexch)$coeff[6,1]
outCWGEEexch_resp1se <- summary(outCWGEEexch)$coeff[1,2]
outCWGEEexch_resp2se <- summary(outCWGEEexch)$coeff[2,2]
outCWGEEexch_resp3se <- summary(outCWGEEexch)$coeff[3,2]
outCWGEEexch_betase <- summary(outCWGEEexch)$coeff[4,2]
outCWGEEexch_beta2se <- summary(outCWGEEexch)$coeff[5,2]
outCWGEEexch_beta3se <- summary(outCWGEEexch)$coeff[6,2]
outCWGEEexch_betareject <- ifelse(summary(outCWGEEexch)$coeff[4,4] < 0.05, 1, 0)
outCWGEEexch_beta2reject <- ifelse(summary(outCWGEEexch)$coeff[5,4] < 0.05, 1, 0)
outCWGEEexch_beta3reject <- ifelse(summary(outCWGEEexch)$coeff[6,4] < 0.05, 1, 0)
outCWGEEexch_alpha12est <- summary(outCWGEEexch)$corr.matrix[1,2]
outCWGEEexch_alpha13est <- summary(outCWGEEexch)$corr.matrix[1,3]
outCWGEEexch_alpha23est <- summary(outCWGEEexch)$corr.matrix[2,3]
### multivariate Wald statistic for H0:beta12=beta13=0
outCWGEEexch_beta2beta3est <- summary(outCWGEEexch)$coeff[5:6,1]
outCWGEEexch_beta2beta3cov <- outCWGEEexch$robust.variance[5:6,5:6]
outCWGEEexch_wald <- as.numeric(t(Lmat %*% outCWGEEexch_beta2beta3est) %*% solve(Lmat %*% outCWGEEexch_beta2beta3cov %*% Lmat) %*% (Lmat %*% outCWGEEexch_beta2beta3est))  
outCWGEEexch_waldp <- 1-pchisq(outCWGEEexch_wald, df=nrow(Lmat))
outCWGEEexch_beta2beta3reject <- ifelse(outCWGEEexch_waldp < 0.05, 1, 0)


## independent correlation with unique slopes for x

outGEEind <- mvoGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "ind")
outGEEind_resp1est <- summary(outGEEind)$coeff[1,1]
outGEEind_resp2est <- summary(outGEEind)$coeff[2,1]
outGEEind_resp3est <- summary(outGEEind)$coeff[3,1]
outGEEind_betaest <- summary(outGEEind)$coeff[4,1]
outGEEind_beta2est <- summary(outGEEind)$coeff[5,1]
outGEEind_beta3est <- summary(outGEEind)$coeff[6,1]
outGEEind_resp1se <- summary(outGEEind)$coeff[1,2]
outGEEind_resp2se <- summary(outGEEind)$coeff[2,2]
outGEEind_resp3se <- summary(outGEEind)$coeff[3,2]
outGEEind_betase <- summary(outGEEind)$coeff[4,2]
outGEEind_beta2se <- summary(outGEEind)$coeff[5,2]
outGEEind_beta3se <- summary(outGEEind)$coeff[6,2]
outGEEind_betareject <- ifelse(summary(outGEEind)$coeff[4,4] < 0.05, 1, 0)
outGEEind_beta2reject <- ifelse(summary(outGEEind)$coeff[5,4] < 0.05, 1, 0)
outGEEind_beta3reject <- ifelse(summary(outGEEind)$coeff[6,4] < 0.05, 1, 0)
outGEEind_alpha12est <- summary(outGEEind)$corr.matrix[1,2]
outGEEind_alpha13est <- summary(outGEEind)$corr.matrix[1,3]
outGEEind_alpha23est <- summary(outGEEind)$corr.matrix[2,3]
### multivariate Wald statistic for H0:beta12=beta13=0
outGEEind_beta2beta3est <- summary(outGEEind)$coeff[5:6,1]
outGEEind_beta2beta3cov <- outGEEind$robust.variance[5:6,5:6]
outGEEind_wald <- as.numeric(t(Lmat %*% outGEEind_beta2beta3est) %*% solve(Lmat %*% outGEEind_beta2beta3cov %*% Lmat) %*% (Lmat %*% outGEEind_beta2beta3est))  
outGEEind_waldp <- 1-pchisq(outGEEind_wald, df=nrow(Lmat))
outGEEind_beta2beta3reject <- ifelse(outGEEind_waldp < 0.05, 1, 0)


outCWGEEind <- mvoCWGEE(formula = y ~ x, data = simdat, cluster = subject, resp.ind = outcome, unit = tooth, corr.str = "ind")
outCWGEEind_resp1est <- summary(outCWGEEind)$coeff[1,1]
outCWGEEind_resp2est <- summary(outCWGEEind)$coeff[2,1]
outCWGEEind_resp3est <- summary(outCWGEEind)$coeff[3,1]
outCWGEEind_betaest <- summary(outCWGEEind)$coeff[4,1]
outCWGEEind_beta2est <- summary(outCWGEEind)$coeff[5,1]
outCWGEEind_beta3est <- summary(outCWGEEind)$coeff[6,1]
outCWGEEind_resp1se <- summary(outCWGEEind)$coeff[1,2]
outCWGEEind_resp2se <- summary(outCWGEEind)$coeff[2,2]
outCWGEEind_resp3se <- summary(outCWGEEind)$coeff[3,2]
outCWGEEind_betase <- summary(outCWGEEind)$coeff[4,2]
outCWGEEind_beta2se <- summary(outCWGEEind)$coeff[5,2]
outCWGEEind_beta3se <- summary(outCWGEEind)$coeff[6,2]
outCWGEEind_betareject <- ifelse(summary(outCWGEEind)$coeff[4,4] < 0.05, 1, 0)
outCWGEEind_beta2reject <- ifelse(summary(outCWGEEind)$coeff[5,4] < 0.05, 1, 0)
outCWGEEind_beta3reject <- ifelse(summary(outCWGEEind)$coeff[6,4] < 0.05, 1, 0)
outCWGEEind_alpha12est <- summary(outCWGEEind)$corr.matrix[1,2]
outCWGEEind_alpha13est <- summary(outCWGEEind)$corr.matrix[1,3]
outCWGEEind_alpha23est <- summary(outCWGEEind)$corr.matrix[2,3]
### multivariate Wald statistic for H0:beta12=beta13=0
outCWGEEind_beta2beta3est <- summary(outCWGEEind)$coeff[5:6,1]
outCWGEEind_beta2beta3cov <- outCWGEEind$robust.variance[5:6,5:6]
outCWGEEind_wald <- as.numeric(t(Lmat %*% outCWGEEind_beta2beta3est) %*% solve(Lmat %*% outCWGEEind_beta2beta3cov %*% Lmat) %*% (Lmat %*% outCWGEEind_beta2beta3est))  
outCWGEEind_waldp <- 1-pchisq(outCWGEEind_wald, df=nrow(Lmat))
outCWGEEind_beta2beta3reject <- ifelse(outCWGEEind_waldp < 0.05, 1, 0)




output <- cbind(RR, SS, theseed, N, nteeth,
                truea1, truea2, truea3, truebeta1, truetau, truealpha12, truealpha13, truealpha23, therho12, therho13, therho23, thetau,
                outGEE1_resp1est, outGEE1_resp1se, outGEE1_betaest, outGEE1_betase, outGEE1_betareject,
                outCWGEE1_resp1est, outCWGEE1_resp1se, outCWGEE1_betaest, outCWGEE1_betase, outCWGEE1_betareject,
                outGEE2_resp1est, outGEE2_resp1se, outGEE2_betaest, outGEE2_betase, outGEE2_betareject,
                outCWGEE2_resp1est, outCWGEE2_resp1se, outCWGEE2_betaest, outCWGEE2_betase, outCWGEE2_betareject,
                outGEE3_resp1est, outGEE3_resp1se, outGEE3_betaest, outGEE3_betase, outGEE3_betareject,
                outCWGEE3_resp1est, outCWGEE3_resp1se, outCWGEE3_betaest, outCWGEE3_betase, outCWGEE3_betareject,
                outGEEunstr_commonslope_resp1est, outGEEunstr_commonslope_resp1se, outGEEunstr_commonslope_resp2est, outGEEunstr_commonslope_resp2se, outGEEunstr_commonslope_resp3est, outGEEunstr_commonslope_resp3se, outGEEunstr_commonslope_betaest, outGEEunstr_commonslope_betase, outGEEunstr_commonslope_betareject, outGEEunstr_commonslope_alpha12est, outGEEunstr_commonslope_alpha13est, outGEEunstr_commonslope_alpha23est,
                outCWGEEunstr_commonslope_resp1est, outCWGEEunstr_commonslope_resp1se, outCWGEEunstr_commonslope_resp2est, outCWGEEunstr_commonslope_resp2se, outCWGEEunstr_commonslope_resp3est, outCWGEEunstr_commonslope_resp3se, outCWGEEunstr_commonslope_betaest, outCWGEEunstr_commonslope_betase, outCWGEEunstr_commonslope_betareject, outCWGEEunstr_commonslope_alpha12est, outCWGEEunstr_commonslope_alpha13est, outCWGEEunstr_commonslope_alpha23est,
                outGEEexch_commonslope_resp1est, outGEEexch_commonslope_resp1se, outGEEexch_commonslope_resp2est, outGEEexch_commonslope_resp2se, outGEEexch_commonslope_resp3est, outGEEexch_commonslope_resp3se, outGEEexch_commonslope_betaest, outGEEexch_commonslope_betase, outGEEexch_commonslope_betareject, outGEEexch_commonslope_alpha12est, outGEEexch_commonslope_alpha13est, outGEEexch_commonslope_alpha23est,
                outCWGEEexch_commonslope_resp1est, outCWGEEexch_commonslope_resp1se, outCWGEEexch_commonslope_resp2est, outCWGEEexch_commonslope_resp2se, outCWGEEexch_commonslope_resp3est, outCWGEEexch_commonslope_resp3se, outCWGEEexch_commonslope_betaest, outCWGEEexch_commonslope_betase, outCWGEEexch_commonslope_betareject, outCWGEEexch_commonslope_alpha12est, outCWGEEexch_commonslope_alpha13est, outCWGEEexch_commonslope_alpha23est,
                outGEEind_commonslope_resp1est, outGEEind_commonslope_resp1se, outGEEind_commonslope_resp2est, outGEEind_commonslope_resp2se, outGEEind_commonslope_resp3est, outGEEind_commonslope_resp3se, outGEEind_commonslope_betaest, outGEEind_commonslope_betase, outGEEind_commonslope_betareject, outGEEind_commonslope_alpha12est, outGEEind_commonslope_alpha13est, outGEEind_commonslope_alpha23est,
                outCWGEEind_commonslope_resp1est, outCWGEEind_commonslope_resp1se, outCWGEEind_commonslope_resp2est, outCWGEEind_commonslope_resp2se, outCWGEEind_commonslope_resp3est, outCWGEEind_commonslope_resp3se, outCWGEEind_commonslope_betaest, outCWGEEind_commonslope_betase, outCWGEEind_commonslope_betareject, outCWGEEind_commonslope_alpha12est, outCWGEEind_commonslope_alpha13est, outCWGEEind_commonslope_alpha23est,
                outGEEunstr_resp1est, outGEEunstr_resp1se, outGEEunstr_resp2est, outGEEunstr_resp2se, outGEEunstr_resp3est, outGEEunstr_resp3se, outGEEunstr_betaest, outGEEunstr_betase, outGEEunstr_betareject, outGEEunstr_beta2est, outGEEunstr_beta2se, outGEEunstr_beta2reject, outGEEunstr_beta3est, outGEEunstr_beta3se, outGEEunstr_beta3reject, outGEEunstr_alpha12est, outGEEunstr_alpha13est, outGEEunstr_alpha23est,
                outGEEunstr_wald, outGEEunstr_waldp, outGEEunstr_beta2beta3reject,
                outCWGEEunstr_resp1est, outCWGEEunstr_resp1se, outCWGEEunstr_resp2est, outCWGEEunstr_resp2se, outCWGEEunstr_resp3est, outCWGEEunstr_resp3se, outCWGEEunstr_betaest, outCWGEEunstr_betase, outCWGEEunstr_betareject, outCWGEEunstr_beta2est, outCWGEEunstr_beta2se, outCWGEEunstr_beta2reject, outCWGEEunstr_beta3est, outCWGEEunstr_beta3se, outCWGEEunstr_beta3reject, outCWGEEunstr_alpha12est, outCWGEEunstr_alpha13est, outCWGEEunstr_alpha23est,
                outCWGEEunstr_wald, outCWGEEunstr_waldp, outCWGEEunstr_beta2beta3reject,
                outGEEexch_resp1est, outGEEexch_resp1se, outGEEexch_resp2est, outGEEexch_resp2se, outGEEexch_resp3est, outGEEexch_resp3se, outGEEexch_betaest, outGEEexch_betase, outGEEexch_betareject, outGEEexch_beta2est, outGEEexch_beta2se, outGEEexch_beta2reject, outGEEexch_beta3est, outGEEexch_beta3se, outGEEexch_beta3reject, outGEEexch_alpha12est, outGEEexch_alpha13est, outGEEexch_alpha23est,
                outGEEexch_wald, outGEEexch_waldp, outGEEexch_beta2beta3reject,
                outCWGEEexch_resp1est, outCWGEEexch_resp1se, outCWGEEexch_resp2est, outCWGEEexch_resp2se, outCWGEEexch_resp3est, outCWGEEexch_resp3se, outCWGEEexch_betaest, outCWGEEexch_betase, outCWGEEexch_betareject, outCWGEEexch_beta2est, outCWGEEexch_beta2se, outCWGEEexch_beta2reject, outCWGEEexch_beta3est, outCWGEEexch_beta3se, outCWGEEexch_beta3reject, outCWGEEexch_alpha12est, outCWGEEexch_alpha13est, outCWGEEexch_alpha23est,
                outCWGEEexch_wald, outCWGEEexch_waldp, outCWGEEexch_beta2beta3reject,
                outGEEind_resp1est, outGEEind_resp1se, outGEEind_resp2est, outGEEind_resp2se, outGEEind_resp3est, outGEEind_resp3se, outGEEind_betaest, outGEEind_betase, outGEEind_betareject, outGEEind_beta2est, outGEEind_beta2se, outGEEind_beta2reject, outGEEind_beta3est, outGEEind_beta3se, outGEEind_beta3reject, outGEEind_alpha12est, outGEEind_alpha13est, outGEEind_alpha23est,
                outGEEind_wald, outGEEind_waldp, outGEEind_beta2beta3reject,
                outCWGEEind_resp1est, outCWGEEind_resp1se, outCWGEEind_resp2est, outCWGEEind_resp2se, outCWGEEind_resp3est, outCWGEEind_resp3se, outCWGEEind_betaest, outCWGEEind_betase, outCWGEEind_betareject, outCWGEEind_beta2est, outCWGEEind_beta2se, outCWGEEind_beta2reject, outCWGEEind_beta3est, outCWGEEind_beta3se, outCWGEEind_beta3reject, outCWGEEind_alpha12est, outCWGEEind_alpha13est, outCWGEEind_alpha23est,
                outCWGEEind_wald, outCWGEEind_waldp, outCWGEEind_beta2beta3reject
)




### for first simulation, output the column names as well

dir.name <- paste0("diroutN",N,"_nteeth",nteeth,"_a1",truea1,"_a2",truea2,"_a3",truea3,
                   "_beta1", truebeta1,
                   "_tau", truetau, "_alpha12", truealpha12, "_alpha13", truealpha13, "_alpha23", truealpha23)
if (!dir.exists (dir.name )) dir.create(dir.name)

out.file <-paste(dir.name, paste0("simoutN",N,"_nteeth",nteeth,"_a1",truea1,"_a2",truea2,"_a3",truea3,
                                  "_beta1", truebeta1,
                                  "_tau", truetau, "_alpha12", truealpha12, "_alpha13", truealpha13, "_alpha23", truealpha23, "_", SS, ".txt"), sep="/")
if(SS==1){
  write.table(output, out.file, col.names = T, row.names = F)
}else{
  write.table(output, out.file, col.names = F, row.names = F)
}

