#----------------------------------------------------------------------------------------
# R code for creating Type I error plots for simulation results in
# "Marginal analysis of multiple outcomes with informative cluster size"
#----------------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)

#setwd("")

### Function to read in simulation output and extract Type I error for GEE and CWGEE unstructured 

toplotfunc <- function(Nsubj, truetau, truealpha12, truealpha13, truealpha23){
  
  outsim0 <- read.table(paste0("alloutN", Nsubj, "_nteeth28", "_a1-1_a2-1.5_a3-2_beta10.5_tau", truetau, "_alpha12", truealpha12, "_alpha13", truealpha13, "_alpha23", truealpha23, ".txt"), header = TRUE)
  
  dim(outsim0)
  
  outsim <- outsim0[1:1000,]
  
  thenames <- c("tau", "alpha", "Parameter", "type1err", "Model")
  
  outsimGEE2 <- outsim[c("truetau", "truealpha12", "outGEEunstr_beta2beta3reject")]
  outsimGEE2$type1err <- mean(outsimGEE2$outGEEunstr_beta2beta3reject)
  outsimGEE2$model <- "GEE"
  outsimGEE2$parameter <- "2"
  outGEE2 <- outsimGEE2[1,c("truetau", "truealpha12", "parameter", "type1err", "model")]
  
  # outsimGEE3 <- outsim[c("truetau", "truealpha12", "outGEEunstr_beta3reject")]
  # outsimGEE3$type1err <- mean(outsimGEE3$outGEEunstr_beta3reject)
  # outsimGEE3$model <- "GEE"
  # outsimGEE3$parameter <- "3"
  # outGEE3 <- outsimGEE3[1,c("truetau", "truealpha12", "parameter", "type1err", "model")]
  
  outsimCWGEE2 <- outsim[c("truetau", "truealpha12", "outCWGEEunstr_beta2beta3reject")]
  outsimCWGEE2$type1err <- mean(outsimCWGEE2$outCWGEEunstr_beta2beta3reject)
  outsimCWGEE2$model <- "CWGEE"
  outsimCWGEE2$parameter <- "2"
  outCWGEE2 <- outsimCWGEE2[1,c("truetau", "truealpha12", "parameter", "type1err", "model")]
  
  # outsimCWGEE3 <- outsim[c("truetau", "truealpha12", "outCWGEEunstr_beta3reject")]
  # outsimCWGEE3$type1err <- mean(outsimCWGEE3$outCWGEEunstr_beta3reject)
  # outsimCWGEE3$model <- "CWGEE"
  # outsimCWGEE3$parameter <- "3"
  # outCWGEE3 <- outsimCWGEE3[1,c("truetau", "truealpha12", "parameter", "type1err", "model")]
  
  out <- rbind(outGEE2, outCWGEE2)
  names(out) <- thenames
  out$tau <- as.character(out$tau)
  out$alpha <- as.character(out$alpha)
  out
  
}

Nsubj <- 750
truealpha12 <- "0"
truealpha13 <- "0"
truealpha23 <- "0"
toplotfunc_no_tau0 <- toplotfunc(Nsubj, truetau = "0", truealpha12, truealpha13, truealpha23)
toplotfunc_no_tau0.25 <- toplotfunc(Nsubj, truetau = "0.25", truealpha12, truealpha13, truealpha23)
toplotfunc_no_tau0.5 <- toplotfunc(Nsubj, truetau = "0.5", truealpha12, truealpha13, truealpha23)
toplotfunc_no_tau0.75 <- toplotfunc(Nsubj, truetau = "0.75", truealpha12, truealpha13, truealpha23)

toplotfunc_no <- rbind(toplotfunc_no_tau0, toplotfunc_no_tau0.25, toplotfunc_no_tau0.5, toplotfunc_no_tau0.75)

truealpha12 <- "0.4"
truealpha13 <- "0.35"
truealpha23 <- "0.3"
toplotfunc_low_tau0 <- toplotfunc(Nsubj, truetau = "0", truealpha12, truealpha13, truealpha23)
toplotfunc_low_tau0.25 <- toplotfunc(Nsubj, truetau = "0.25", truealpha12, truealpha13, truealpha23)
toplotfunc_low_tau0.5 <- toplotfunc(Nsubj, truetau = "0.5", truealpha12, truealpha13, truealpha23)
toplotfunc_low_tau0.75 <- toplotfunc(Nsubj, truetau = "0.75", truealpha12, truealpha13, truealpha23)

toplotfunc_low <- rbind(toplotfunc_low_tau0, toplotfunc_low_tau0.25, toplotfunc_low_tau0.5, toplotfunc_low_tau0.75)

truealpha12 <- "0.6"
truealpha13 <- "0.55"
truealpha23 <- "0.5"
toplotfunc_med_tau0 <- toplotfunc(Nsubj, truetau = "0", truealpha12, truealpha13, truealpha23)
toplotfunc_med_tau0.25 <- toplotfunc(Nsubj, truetau = "0.25", truealpha12, truealpha13, truealpha23)
toplotfunc_med_tau0.5 <- toplotfunc(Nsubj, truetau = "0.5", truealpha12, truealpha13, truealpha23)
toplotfunc_med_tau0.75 <- toplotfunc(Nsubj, truetau = "0.75", truealpha12, truealpha13, truealpha23)

toplotfunc_med <- rbind(toplotfunc_med_tau0, toplotfunc_med_tau0.25, toplotfunc_med_tau0.5, toplotfunc_med_tau0.75)

truealpha12 <- "0.8"
truealpha13 <- "0.75"
truealpha23 <- "0.7"
toplotfunc_hi_tau0 <- toplotfunc(Nsubj, truetau = "0", truealpha12, truealpha13, truealpha23)
toplotfunc_hi_tau0.25 <- toplotfunc(Nsubj, truetau = "0.25", truealpha12, truealpha13, truealpha23)
toplotfunc_hi_tau0.5 <- toplotfunc(Nsubj, truetau = "0.5", truealpha12, truealpha13, truealpha23)
toplotfunc_hi_tau0.75 <- toplotfunc(Nsubj, truetau = "0.75", truealpha12, truealpha13, truealpha23)

toplotfunc_hi <- rbind(toplotfunc_hi_tau0, toplotfunc_hi_tau0.25, toplotfunc_hi_tau0.5, toplotfunc_hi_tau0.75)


#ggplot(toplotfunc_no, aes(x=tau, y=relbias, fill = Model)) + geom_split_violin() + theme_minimal() 
#ggplot(toplotfunc_low, aes(x=tau, y=relbias, fill = Model)) + geom_split_violin() + theme_minimal() 
#ggplot(toplotfunc_med, aes(x=tau, y=relbias, fill = Model)) + geom_split_violin() + theme_minimal() 
#ggplot(toplotfunc_hi, aes(x=tau, y=relbias, fill = Model)) + geom_split_violin() + theme_minimal() 


toplot <- rbind(toplotfunc_no, toplotfunc_low, toplotfunc_med, toplotfunc_hi)
toplot$taulabel <- factor(toplot$tau, labels = c("tau==0", "tau==0.25", "tau==0.5", "tau==0.75"))
toplot$paramlabel <- factor(toplot$Parameter, labels = c("H[0]:beta[12]==beta[13]==0"))


plotN750 <- 
  ggplot(toplot, aes(x=alpha, y=type1err, group = Model)) + 
  geom_point(size=3, aes(shape = Model, color = Model)) +
  theme_bw() +
  ylab("Type I error") + 
  xlab("Correlation among outcomes") +
  scale_x_discrete(breaks=c("0", "0.4", "0.6", "0.8"), labels=c("None", "Low", "Medium", "High")) + 
  geom_hline(yintercept = 0.05, linetype = "longdash") +
  #geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0.04, ymax = 0.06), alpha = 0.03) + 
  facet_wrap(~taulabel, labeller = label_parsed) 
plotN750


### in black and white
plotN750bw <- 
  ggplot(toplot, aes(x=alpha, y=type1err, group = Model)) + 
  geom_point(size=3, aes(shape=Model)) +
  scale_shape_manual(values=c(16, 21)) +
  #scale_color_grey(start = 0.8, end = 1) + 
  theme_bw(base_size = 15) +
  ylab("Type I error") + 
  xlab("Correlation among outcomes") +
  scale_x_discrete(breaks=c("0", "0.4", "0.6", "0.8"), labels=c("None", "Low", "Medium", "High")) + 
  geom_hline(yintercept = 0.05, linetype = "longdash") +
  facet_wrap(~taulabel, labeller = label_parsed) 
plotN750bw

