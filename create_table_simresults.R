#----------------------------------------------------------------------------------------
# R code for creating tables for simulation results in
# "Marginal analysis of multiple outcomes with informative cluster size"
#----------------------------------------------------------------------------------------


#--------------------------------
# load packages
#--------------------------------
library(xtable)


## set working directory where the simulation results exist
#setwd("")

## assign parameters to choose simulation scenario
Nsubj <- 750
maxNteeth <- 28
Noutcome <- 3
truetau <- "0.5"
truealpha12 <- "0.4"
truealpha13 <- "0.35"
truealpha23 <- "0.3"

outsim0 <- read.table(paste0("alloutN", Nsubj, "_nteeth", maxNteeth, "_a1-1_a2-1.5_a3-2_beta10.5_tau", truetau, "_alpha12", truealpha12, "_alpha13", truealpha13, "_alpha23", truealpha23, ".txt"), header = TRUE)


#dim(outsim0)

#outsim <- outsim0[1:1000,]


## Function to compute bias, MSE, mean and coverage for all parameter estimates
simresults <- function(truebetavec, estbetavec, estsebetavec, covervec){
  meanbetavec <- apply(estbetavec, 2, mean, na.rm = TRUE)
  meansebetavec <- apply(estsebetavec, 2, mean, na.rm = TRUE)
  sdbetavec <- apply(estbetavec, 2, sd, na.rm = TRUE)
  #bias <- apply(estbetavec - truebetavec, 2, mean, na.rm = TRUE)
  #relbias <- apply((estbetavec - truebetavec)/truebetavec, 2, mean, na.rm = TRUE)
  relbias <- apply((estbetavec - truebetavec), 2, mean, na.rm = TRUE) # if truth = 0 then use bias instead of relative bias
  #mse <- apply((estbetavec - truebetavec) ^ 2, 2, mean, na.rm = TRUE)
  covpvec <- covervec * 100
  #simresultslist <- list(bias, relbias, mse, meanbetavec, meansebetavec, sdbetavec, covpvec)
  simresultslist <- list(meanbetavec, meansebetavec, sdbetavec, relbias, covpvec)
  return(simresultslist)
}

## Functon to compute mean type I error for beta
type1errorfunc <- function(rejectvec){
  type1errorvec <- apply(rejectvec, 2, mean, na.rm = TRUE)
}

## Function to compute mean correlation
meancorrfunc <- function(corrmatrixvalues){
  meancorrvec <- apply(corrmatrixvalues, 2, mean, na.rm = TRUE)
}

## Function to compute coverage probability
covpfunc <- function(estbetavec, estsebetavec, truebetavec){
  covprob <- c()
  for(i in 1:length(estbetavec)){
    lcl <- estbetavec[i] - qnorm(0.975) * estsebetavec[i]
    ucl <- estbetavec[i] + qnorm(0.975) * estsebetavec[i]
    covprob[i] <- mean(ifelse(lcl[,1] < truebetavec[i] & truebetavec[i] < ucl[,1], 1, 0), na.rm = TRUE)
  }
  covprob
}

colname_results <- c("GEE Unstr", "GEE Exch", "GEE Ind", "CWGEE Unstr", "CWGEE Exch", "CWGEE Ind")

truebeta2 <- rep(0, 1000)
truebeta3 <- rep(0, 1000)
truebetavec <- cbind(outsim[c("truea1", "truea2", "truea3", "truebeta1")], truebeta2, truebeta3)
truebetavec_commonslope <- outsim[c("truea1", "truea2", "truea3", "truebeta1")]
truebetas <- c(truebetavec[1,1], truebetavec[1,2], truebetavec[1,3], truebetavec[1,4], as.numeric(truealpha12), as.numeric(truealpha13), as.numeric(truealpha23))



GEEresp1_betaest <- outsim[c("outGEE1_resp1est", "outGEE1_betaest")]
GEEresp1_seest <- outsim[c("outGEE1_resp1se", "outGEE1_betase")]

GEEresp2_betaest <- outsim[c("outGEE2_resp1est", "outGEE2_betaest")]
GEEresp2_seest <- outsim[c("outGEE2_resp1se", "outGEE2_betase")]

GEEresp3_betaest <- outsim[c("outGEE3_resp1est", "outGEE3_betaest")]
GEEresp3_seest <- outsim[c("outGEE3_resp1se", "outGEE3_betase")]

CWGEEresp1_betaest <- outsim[c("outCWGEE1_resp1est", "outCWGEE1_betaest")]
CWGEEresp1_seest <- outsim[c("outCWGEE1_resp1se", "outCWGEE1_betase")]

CWGEEresp2_betaest <- outsim[c("outCWGEE2_resp1est", "outCWGEE2_betaest")]
CWGEEresp2_seest <- outsim[c("outCWGEE2_resp1se", "outCWGEE2_betase")]

CWGEEresp3_betaest <- outsim[c("outCWGEE3_resp1est", "outCWGEE3_betaest")]
CWGEEresp3_seest <- outsim[c("outCWGEE3_resp1se", "outCWGEE3_betase")]

GEEunstr_betaest <- outsim[c("outGEEunstr_resp1est", "outGEEunstr_resp2est", "outGEEunstr_resp3est", "outGEEunstr_betaest", "outGEEunstr_beta2est", "outGEEunstr_beta3est")]
GEEunstr_betase <- outsim[c("outGEEunstr_resp1se", "outGEEunstr_resp2se", "outGEEunstr_resp3se", "outGEEunstr_betase", "outGEEunstr_beta2se", "outGEEunstr_beta3se")]
GEEunstr_alphaest <- outsim[c("outGEEunstr_alpha12est", "outGEEunstr_alpha13est", "outGEEunstr_alpha23est")]
GEEunstr_covervec <- covpfunc(estbetavec=GEEunstr_betaest, estsebetavec=GEEunstr_betase, truebetavec=truebetavec)
GEEunstr_simresults <- simresults(truebetavec=truebetavec, estbetavec=GEEunstr_betaest, estsebetavec=GEEunstr_betase, covervec=GEEunstr_covervec)
GEEunstr_type1error <- type1errorfunc(outsim[c("outGEEunstr_beta2beta3reject", "outGEEunstr_beta2reject", "outGEEunstr_beta3reject")])
GEEunstr_type1error <- c(rep(NA, 4), GEEunstr_type1error)
GEEunstr_meancorr <- meancorrfunc(GEEunstr_alphaest)

GEEexch_betaest <- outsim[c("outGEEexch_resp1est", "outGEEexch_resp2est", "outGEEexch_resp3est", "outGEEexch_betaest", "outGEEexch_beta2est", "outGEEexch_beta3est")]
GEEexch_betase <- outsim[c("outGEEexch_resp1se", "outGEEexch_resp2se", "outGEEexch_resp3se", "outGEEexch_betase", "outGEEexch_beta2se", "outGEEexch_beta3se")]
GEEexch_alphaest <- outsim[c("outGEEexch_alpha12est", "outGEEexch_alpha13est", "outGEEexch_alpha23est")]
GEEexch_covervec <- covpfunc(estbetavec=GEEexch_betaest, estsebetavec=GEEexch_betase, truebetavec=truebetavec)
GEEexch_simresults <- simresults(truebetavec=truebetavec, estbetavec=GEEexch_betaest, estsebetavec=GEEexch_betase, covervec=GEEexch_covervec)
GEEexch_type1error <- type1errorfunc(outsim[c("outGEEexch_beta2beta3reject", "outGEEexch_beta2reject", "outGEEexch_beta3reject")])
GEEexch_type1error <- c(rep(NA, 4), GEEexch_type1error)
GEEexch_meancorr <- meancorrfunc(GEEexch_alphaest)

GEEind_betaest <- outsim[c("outGEEind_resp1est", "outGEEind_resp2est", "outGEEind_resp3est", "outGEEind_betaest", "outGEEind_beta2est", "outGEEind_beta3est")]
GEEind_betase <- outsim[c("outGEEind_resp1se", "outGEEind_resp2se", "outGEEind_resp3se", "outGEEind_betase", "outGEEind_beta2se", "outGEEind_beta3se")]
GEEind_alphaest <- outsim[c("outGEEind_alpha12est", "outGEEind_alpha13est", "outGEEind_alpha23est")]
GEEind_covervec <- covpfunc(estbetavec=GEEind_betaest, estsebetavec=GEEind_betase, truebetavec=truebetavec)
GEEind_simresults <- simresults(truebetavec=truebetavec, estbetavec=GEEind_betaest, estsebetavec=GEEind_betase, covervec=GEEind_covervec)
GEEind_type1error <- type1errorfunc(outsim[c("outGEEind_beta2beta3reject", "outGEEind_beta2reject", "outGEEind_beta3reject")])
GEEind_type1error <- c(rep(NA, 4), GEEind_type1error)
GEEind_meancorr <- meancorrfunc(GEEind_alphaest)

CWGEEunstr_betaest <- outsim[c("outCWGEEunstr_resp1est", "outCWGEEunstr_resp2est", "outCWGEEunstr_resp3est", "outCWGEEunstr_betaest", "outCWGEEunstr_beta2est", "outCWGEEunstr_beta3est")]
CWGEEunstr_betase <- outsim[c("outCWGEEunstr_resp1se", "outCWGEEunstr_resp2se", "outCWGEEunstr_resp3se", "outCWGEEunstr_betase", "outCWGEEunstr_beta2se", "outCWGEEunstr_beta3se")]
CWGEEunstr_alphaest <- outsim[c("outCWGEEunstr_alpha12est", "outCWGEEunstr_alpha13est", "outCWGEEunstr_alpha23est")]
CWGEEunstr_covervec <- covpfunc(estbetavec=CWGEEunstr_betaest, estsebetavec=CWGEEunstr_betase, truebetavec=truebetavec)
CWGEEunstr_simresults <- simresults(truebetavec=truebetavec, estbetavec=CWGEEunstr_betaest, estsebetavec=CWGEEunstr_betase, covervec=CWGEEunstr_covervec)
CWGEEunstr_type1error <- type1errorfunc(outsim[c("outCWGEEunstr_beta2beta3reject", "outCWGEEunstr_beta2reject", "outCWGEEunstr_beta3reject")])
CWGEEunstr_type1error <- c(rep(NA, 4), CWGEEunstr_type1error)
CWGEEunstr_meancorr <- meancorrfunc(CWGEEunstr_alphaest)

CWGEEexch_betaest <- outsim[c("outCWGEEexch_resp1est", "outCWGEEexch_resp2est", "outCWGEEexch_resp3est", "outCWGEEexch_betaest", "outCWGEEexch_beta2est", "outCWGEEexch_beta3est")]
CWGEEexch_betase <- outsim[c("outCWGEEexch_resp1se", "outCWGEEexch_resp2se", "outCWGEEexch_resp3se", "outCWGEEexch_betase", "outCWGEEexch_beta2se", "outCWGEEexch_beta3se")]
CWGEEexch_alphaest <- outsim[c("outCWGEEexch_alpha12est", "outCWGEEexch_alpha13est", "outCWGEEexch_alpha23est")]
CWGEEexch_covervec <- covpfunc(estbetavec=CWGEEexch_betaest, estsebetavec=CWGEEexch_betase, truebetavec=truebetavec)
CWGEEexch_simresults <- simresults(truebetavec=truebetavec, estbetavec=CWGEEexch_betaest, estsebetavec=CWGEEexch_betase, covervec=CWGEEexch_covervec)
CWGEEexch_type1error <- type1errorfunc(outsim[c("outCWGEEexch_beta2beta3reject", "outCWGEEexch_beta2reject", "outCWGEEexch_beta3reject")])
CWGEEexch_type1error <- c(rep(NA, 4), CWGEEexch_type1error)
CWGEEexch_meancorr <- meancorrfunc(CWGEEexch_alphaest)

CWGEEind_betaest <- outsim[c("outCWGEEind_resp1est", "outCWGEEind_resp2est", "outCWGEEind_resp3est", "outCWGEEind_betaest", "outCWGEEind_beta2est", "outCWGEEind_beta3est")]
CWGEEind_betase <- outsim[c("outCWGEEind_resp1se", "outCWGEEind_resp2se", "outCWGEEind_resp3se", "outCWGEEind_betase", "outCWGEEind_beta2se", "outCWGEEind_beta3se")]
CWGEEind_alphaest <- outsim[c("outCWGEEind_alpha12est", "outCWGEEind_alpha13est", "outCWGEEind_alpha23est")]
CWGEEind_covervec <- covpfunc(estbetavec=CWGEEind_betaest, estsebetavec=CWGEEind_betase, truebetavec=truebetavec)
CWGEEind_simresults <- simresults(truebetavec=truebetavec, estbetavec=CWGEEind_betaest, estsebetavec=CWGEEind_betase, covervec=CWGEEind_covervec)
CWGEEind_type1error <- type1errorfunc(outsim[c("outCWGEEind_beta2beta3reject", "outCWGEEind_beta2reject", "outCWGEEind_beta3reject")])
CWGEEind_type1error <- c(rep(NA, 4), CWGEEind_type1error)
CWGEEind_meancorr <- meancorrfunc(CWGEEind_alphaest)


GEEunstr_results <- as.data.frame(rbind(do.call("rbind", GEEunstr_simresults), GEEunstr_type1error), row.names = FALSE, col.names = FALSE)
GEEunstr_resultslong <- unlist(GEEunstr_results, use.names = FALSE)

GEEexch_results <- as.data.frame(rbind(do.call("rbind", GEEexch_simresults), GEEexch_type1error), row.names = FALSE, col.names = FALSE)
GEEexch_resultslong <- unlist(GEEexch_results, use.names = FALSE)

GEEind_results <- as.data.frame(rbind(do.call("rbind", GEEind_simresults), GEEind_type1error), row.names = FALSE, col.names = FALSE)
GEEind_resultslong <- unlist(GEEind_results, use.names = FALSE)


CWGEEunstr_results <- as.data.frame(rbind(do.call("rbind", CWGEEunstr_simresults), CWGEEunstr_type1error), row.names = FALSE, col.names = FALSE)
CWGEEunstr_resultslong <- unlist(CWGEEunstr_results, use.names = FALSE)

CWGEEexch_results <- as.data.frame(rbind(do.call("rbind", CWGEEexch_simresults), CWGEEexch_type1error), row.names = FALSE, col.names = FALSE)
CWGEEexch_resultslong <- unlist(CWGEEexch_results, use.names = FALSE)

CWGEEind_results <- as.data.frame(rbind(do.call("rbind", CWGEEind_simresults), CWGEEind_type1error), row.names = FALSE, col.names = FALSE)
CWGEEind_resultslong <- unlist(CWGEEind_results, use.names = FALSE)

allresults0 <- cbind(GEEunstr_resultslong, GEEexch_resultslong, GEEind_resultslong, CWGEEunstr_resultslong, CWGEEexch_resultslong, CWGEEind_resultslong)



#simresults0 <- c("Mean Est", "Mean SE", "SD Est", "Rel Bias", "Covp")
simresults1 <- c("Mean Est", "Mean SE", "SD Est", "Rel Bias", "Covp", "Type I Err")
simresultsheader <- c(rep(simresults1, 6), c("Mean Est", "Mean Est", "Mean Est"))

truebetas <- c(truebetavec[1,1], truebetavec[1,2], truebetavec[1,3], truebetavec[1,4], truebetavec[1,5], truebetavec[1,6], as.numeric(truealpha12), as.numeric(truealpha13), as.numeric(truealpha23))

truevalues <- Truth <- c(truebetas[1], rep("", 5),
                         truebetas[2], rep("", 5),
                         truebetas[3], rep("", 5),
                         truebetas[4], rep("", 5),
                         truebetas[5], rep("", 5),
                         truebetas[6], rep("", 5),
                         truebetas[7], truebetas[8], truebetas[9]
)

Parameter <- c("beta0", "", "", "", "", "",
               "beta1", "", "", "", "", "",
               "beta2", "", "", "", "", "",
               "beta3", "", "", "", "", "",
               "beta4", "", "", "", "", "",
               "beta5", "", "", "", "", "",
               "alpha12", "alpha13", "alpha23"
)

allresults0 <- round(cbind(GEEunstr_resultslong, GEEexch_resultslong, GEEind_resultslong, CWGEEunstr_resultslong, CWGEEexch_resultslong, CWGEEind_resultslong), 3)


alphaests <- round(cbind(GEEunstr_meancorr, GEEexch_meancorr, GEEind_meancorr, CWGEEunstr_meancorr, CWGEEexch_meancorr, CWGEEind_meancorr), 3)
allresults <- rbind(allresults0, alphaests)

allresults_final <- as.data.frame(cbind(Parameter, Truth, simresultsheader, allresults), row.names = FALSE, col.names = FALSE)

names(allresults_final) <- c("Parameter", "Truth", "SimResults", colname_results)


descr <- paste0("Nsubj = ", Nsubj, ", ",
                "maxNteeth = ", maxNteeth, ", ",
                "truealpha = (", truealpha12, ", ", truealpha13, ", ", truealpha23, "), ",
                "truetau = ", truetau, ", ",
                "Nsimulation = ", 1000)

thextable <- xtable(allresults_final, caption = descr)
print(thextable, include.rownames = FALSE)




#--------------------------------------------------------
# Simulation results for common slope scenarios
#--------------------------------------------------------

GEEunstr_commonslope_betaest <- outsim[c("outGEEunstr_commonslope_resp1est", "outGEEunstr_commonslope_resp2est", "outGEEunstr_commonslope_resp3est", "outGEEunstr_commonslope_betaest")]
GEEunstr_commonslope_betase <- outsim[c("outGEEunstr_commonslope_resp1se", "outGEEunstr_commonslope_resp2se", "outGEEunstr_commonslope_resp3se", "outGEEunstr_commonslope_betase")]
GEEunstr_commonslope_alphaest <- outsim[c("outGEEunstr_commonslope_alpha12est", "outGEEunstr_commonslope_alpha13est", "outGEEunstr_commonslope_alpha23est")]
GEEunstr_commonslope_covervec <- covpfunc(estbetavec=GEEunstr_commonslope_betaest, estsebetavec=GEEunstr_commonslope_betase, truebetavec=truebetavec_commonslope)
GEEunstr_commonslope_simresults <- simresults(truebetavec=truebetavec_commonslope, estbetavec=GEEunstr_commonslope_betaest, estsebetavec=GEEunstr_commonslope_betase, covervec=GEEunstr_commonslope_covervec)
GEEunstr_commonslope_meancorr <- meancorrfunc(GEEunstr_commonslope_alphaest)

GEEexch_commonslope_betaest <- outsim[c("outGEEexch_commonslope_resp1est", "outGEEexch_commonslope_resp2est", "outGEEexch_commonslope_resp3est", "outGEEexch_commonslope_betaest")]
GEEexch_commonslope_betase <- outsim[c("outGEEexch_commonslope_resp1se", "outGEEexch_commonslope_resp2se", "outGEEexch_commonslope_resp3se", "outGEEexch_commonslope_betase")]
GEEexch_commonslope_alphaest <- outsim[c("outGEEexch_commonslope_alpha12est", "outGEEexch_commonslope_alpha13est", "outGEEexch_commonslope_alpha23est")]
GEEexch_commonslope_covervec <- covpfunc(estbetavec=GEEexch_commonslope_betaest, estsebetavec=GEEexch_commonslope_betase, truebetavec=truebetavec_commonslope)
GEEexch_commonslope_simresults <- simresults(truebetavec=truebetavec_commonslope, estbetavec=GEEexch_commonslope_betaest, estsebetavec=GEEexch_commonslope_betase, covervec=GEEexch_commonslope_covervec)
GEEexch_commonslope_meancorr <- meancorrfunc(GEEexch_commonslope_alphaest)

GEEind_commonslope_betaest <- outsim[c("outGEEind_commonslope_resp1est", "outGEEind_commonslope_resp2est", "outGEEind_commonslope_resp3est", "outGEEind_commonslope_betaest")]
GEEind_commonslope_betase <- outsim[c("outGEEind_commonslope_resp1se", "outGEEind_commonslope_resp2se", "outGEEind_commonslope_resp3se", "outGEEind_commonslope_betase")]
GEEind_commonslope_alphaest <- outsim[c("outGEEind_commonslope_alpha12est", "outGEEind_commonslope_alpha13est", "outGEEind_commonslope_alpha23est")]
GEEind_commonslope_covervec <- covpfunc(estbetavec=GEEind_commonslope_betaest, estsebetavec=GEEind_commonslope_betase, truebetavec=truebetavec_commonslope)
GEEind_commonslope_simresults <- simresults(truebetavec=truebetavec_commonslope, estbetavec=GEEind_commonslope_betaest, estsebetavec=GEEind_commonslope_betase, covervec=GEEind_commonslope_covervec)
GEEind_commonslope_meancorr <- meancorrfunc(GEEind_commonslope_alphaest)

CWGEEunstr_commonslope_betaest <- outsim[c("outCWGEEunstr_commonslope_resp1est", "outCWGEEunstr_commonslope_resp2est", "outCWGEEunstr_commonslope_resp3est", "outCWGEEunstr_commonslope_betaest")]
CWGEEunstr_commonslope_betase <- outsim[c("outCWGEEunstr_commonslope_resp1se", "outCWGEEunstr_commonslope_resp2se", "outCWGEEunstr_commonslope_resp3se", "outCWGEEunstr_commonslope_betase")]
CWGEEunstr_commonslope_alphaest <- outsim[c("outCWGEEunstr_commonslope_alpha12est", "outCWGEEunstr_commonslope_alpha13est", "outCWGEEunstr_commonslope_alpha23est")]
CWGEEunstr_commonslope_covervec <- covpfunc(estbetavec=CWGEEunstr_commonslope_betaest, estsebetavec=CWGEEunstr_commonslope_betase, truebetavec=truebetavec_commonslope)
CWGEEunstr_commonslope_simresults <- simresults(truebetavec=truebetavec_commonslope, estbetavec=CWGEEunstr_commonslope_betaest, estsebetavec=CWGEEunstr_commonslope_betase, covervec=CWGEEunstr_commonslope_covervec)
CWGEEunstr_commonslope_meancorr <- meancorrfunc(CWGEEunstr_commonslope_alphaest)

CWGEEexch_commonslope_betaest <- outsim[c("outCWGEEexch_commonslope_resp1est", "outCWGEEexch_commonslope_resp2est", "outCWGEEexch_commonslope_resp3est", "outCWGEEexch_commonslope_betaest")]
CWGEEexch_commonslope_betase <- outsim[c("outCWGEEexch_commonslope_resp1se", "outCWGEEexch_commonslope_resp2se", "outCWGEEexch_commonslope_resp3se", "outCWGEEexch_commonslope_betase")]
CWGEEexch_commonslope_alphaest <- outsim[c("outCWGEEexch_commonslope_alpha12est", "outCWGEEexch_commonslope_alpha13est", "outCWGEEexch_commonslope_alpha23est")]
CWGEEexch_commonslope_covervec <- covpfunc(estbetavec=CWGEEexch_commonslope_betaest, estsebetavec=CWGEEexch_commonslope_betase, truebetavec=truebetavec_commonslope)
CWGEEexch_commonslope_simresults <- simresults(truebetavec=truebetavec_commonslope, estbetavec=CWGEEexch_commonslope_betaest, estsebetavec=CWGEEexch_commonslope_betase, covervec=CWGEEexch_commonslope_covervec)
CWGEEexch_commonslope_meancorr <- meancorrfunc(CWGEEexch_commonslope_alphaest)

CWGEEind_commonslope_betaest <- outsim[c("outCWGEEind_commonslope_resp1est", "outCWGEEind_commonslope_resp2est", "outCWGEEind_commonslope_resp3est", "outCWGEEind_commonslope_betaest")]
CWGEEind_commonslope_betase <- outsim[c("outCWGEEind_commonslope_resp1se", "outCWGEEind_commonslope_resp2se", "outCWGEEind_commonslope_resp3se", "outCWGEEind_commonslope_betase")]
CWGEEind_commonslope_alphaest <- outsim[c("outCWGEEind_commonslope_alpha12est", "outCWGEEind_commonslope_alpha13est", "outCWGEEind_commonslope_alpha23est")]
CWGEEind_commonslope_covervec <- covpfunc(estbetavec=CWGEEind_commonslope_betaest, estsebetavec=CWGEEind_commonslope_betase, truebetavec=truebetavec_commonslope)
CWGEEind_commonslope_simresults <- simresults(truebetavec=truebetavec_commonslope, estbetavec=CWGEEind_commonslope_betaest, estsebetavec=CWGEEind_commonslope_betase, covervec=CWGEEind_commonslope_covervec)
CWGEEind_commonslope_meancorr <- meancorrfunc(CWGEEind_commonslope_alphaest)


GEEunstr_commonslope_results <- as.data.frame(rbind(do.call("rbind", GEEunstr_commonslope_simresults)), row.names = FALSE, col.names = FALSE)
GEEunstr_commonslope_resultslong <- unlist(GEEunstr_commonslope_results, use.names = FALSE)

GEEexch_commonslope_results <- as.data.frame(rbind(do.call("rbind", GEEexch_commonslope_simresults)), row.names = FALSE, col.names = FALSE)
GEEexch_commonslope_resultslong <- unlist(GEEexch_commonslope_results, use.names = FALSE)

GEEind_commonslope_results <- as.data.frame(rbind(do.call("rbind", GEEind_commonslope_simresults)), row.names = FALSE, col.names = FALSE)
GEEind_commonslope_resultslong <- unlist(GEEind_commonslope_results, use.names = FALSE)


CWGEEunstr_commonslope_results <- as.data.frame(rbind(do.call("rbind", CWGEEunstr_commonslope_simresults)), row.names = FALSE, col.names = FALSE)
CWGEEunstr_commonslope_resultslong <- unlist(CWGEEunstr_commonslope_results, use.names = FALSE)

CWGEEexch_commonslope_results <- as.data.frame(rbind(do.call("rbind", CWGEEexch_commonslope_simresults)), row.names = FALSE, col.names = FALSE)
CWGEEexch_commonslope_resultslong <- unlist(CWGEEexch_commonslope_results, use.names = FALSE)

CWGEEind_commonslope_results <- as.data.frame(rbind(do.call("rbind", CWGEEind_commonslope_simresults)), row.names = FALSE, col.names = FALSE)
CWGEEind_commonslope_resultslong <- unlist(CWGEEind_commonslope_results, use.names = FALSE)

simresults_commonslope0 <- c("Mean Est", "Mean SE", "SD Est", "Rel Bias", "Covp")
simresults_commonslope <- c(rep(simresults_commonslope0, 4), c("Mean Est", "Mean Est", "Mean Est"))

truevalues <- Truth <- c(truebetas[1], rep("", 4),
                         truebetas[2], rep("", 4),
                         truebetas[3], rep("", 4),
                         truebetas[4], rep("", 4),
                         truebetas[5], truebetas[6], truebetas[7]
)

Parameter <- c("beta0", "", "", "", "", 
               "beta1", "", "", "", "", 
               "beta2", "", "", "", "", 
               "beta3", "", "", "", "",
               "alpha12", "alpha13", "alpha23"
)

allresults0_commonslope <- round(cbind(GEEunstr_commonslope_resultslong, GEEexch_commonslope_resultslong, GEEind_commonslope_resultslong, CWGEEunstr_commonslope_resultslong, CWGEEexch_commonslope_resultslong, CWGEEind_commonslope_resultslong), 3)


alphaests_commonslope <- round(cbind(GEEunstr_commonslope_meancorr, GEEexch_commonslope_meancorr, GEEind_commonslope_meancorr, CWGEEunstr_commonslope_meancorr, CWGEEexch_commonslope_meancorr, CWGEEind_commonslope_meancorr), 3)
allresults_commonslope <- rbind(allresults0_commonslope, alphaests_commonslope)

allresults_commonslope_final <- as.data.frame(cbind(Parameter, Truth, simresults_commonslope, allresults_commonslope), row.names = FALSE, col.names = FALSE)

names(allresults_commonslope_final) <- c("Parameter", "Truth", "SimResults", colname_results)


descr <- paste0("Nsubj = ", Nsubj, ", ",
                "maxNteeth = ", maxNteeth, ", ",
                "truealpha = (", truealpha12, ", ", truealpha13, ", ", truealpha23, "), ",
                "truetau = ", truetau, ", ",
                "Nsimulation = ", 1000)

thextable_mult <- xtable(allresults_commonslope_final, caption = descr)


type1errors <- as.data.frame(cbind(GEEunstr_type1error, GEEexch_type1error, GEEind_type1error, CWGEEunstr_type1error, CWGEEexch_type1error, CWGEEind_type1error)[5:7, ])
rownames(type1errors) <- c("beta2", "beta3", "beta2beta3")
colnames(type1errors) <- colname_results
thextable_type1 <- xtable(type1errors, caption = descr)

print(thextable_type1)
print(thextable_mult, include.rownames = FALSE)

