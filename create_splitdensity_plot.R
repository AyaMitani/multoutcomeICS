#----------------------------------------------------------------------------------------
# R code for creating split density plots for simulation results in
# "Marginal analysis of multiple outcomes with informative cluster size"
# Referenced from: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
#----------------------------------------------------------------------------------------


library(dplyr)
library(ggplot2)



#setwd("")


GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# ggplot(my_data, aes(x, y, fill = m)) + geom_split_violin() + theme_minimal() 



### Function to read in simulation output and extract relative risks for GEE and CWGEE unstructured 

toplotfunc <- function(Nsubj, truetau, truealpha12, truealpha13, truealpha23){
  
  outsim0 <- read.table(paste0("alloutN", Nsubj, "_nteeth28", "_a1-1_a2-1.5_a3-2_beta10.5_tau", truetau, "_alpha12", truealpha12, "_alpha13", truealpha13, "_alpha23", truealpha23, ".txt"), header = TRUE)
  
  dim(outsim0)
  
  outsim <- outsim0[1:1000,]
  
  thenames <- c("tau", "alpha", "relbias", "Model")
  outsimGEE <- outsim[c("truetau", "truealpha12", "outGEEunstr_commonslope_betaest", "truebeta1")]
  outsimGEE$relbias <- (outsimGEE$outGEEunstr_commonslope_betaest - outsimGEE$truebeta1) / outsimGEE$truebeta1
  outsimGEE$model <- "GEE"
  outGEE <- outsimGEE[c("truetau", "truealpha12", "relbias", "model")]
  
  outsimCWGEE <- outsim[c("truetau", "truealpha12", "outCWGEEunstr_commonslope_betaest", "truebeta1")]
  outsimCWGEE$relbias <- (outsimCWGEE$outCWGEEunstr_commonslope_betaest - outsimCWGEE$truebeta1) / outsimCWGEE$truebeta1
  outsimCWGEE$model <- "CWGEE"
  outCWGEE <- outsimCWGEE[c("truetau", "truealpha12", "relbias", "model")]
  
  out <- rbind(outGEE, outCWGEE)
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


toplot <- rbind(toplotfunc_no, toplotfunc_low, toplotfunc_med, toplotfunc_hi)
toplot$taulabel <- factor(toplot$tau, labels = c("tau==0", "tau==0.25", "tau==0.5", "tau==0.75"))

plotN750 <- 
ggplot(toplot, aes(x=alpha, y=relbias, fill = Model, color = Model)) + 
  geom_split_violin() + 
  theme_bw() +
  ylab("Relative Bias"~(beta)) + 
  xlab("Correlation among outcomes") +
  scale_x_discrete(breaks=c("0", "0.4", "0.6", "0.8"), labels=c("None", "Low", "Medium", "High")) + 
  geom_hline(yintercept = 0, linetype = "longdash", color = "black") +
  facet_wrap(~taulabel, labeller = label_parsed)
plotN750


### in black and white
plotN750bw <- 
  ggplot(toplot, aes(x=alpha, y=relbias, fill = Model)) + 
  scale_fill_grey(start = 0.4, end = 1) + 
  geom_split_violin() + 
  #theme_minimal(base_size = 14) +
  theme_bw() + 
  ylab("Relative Bias "~(beta)) + 
  xlab("Correlation among outcomes") +
  scale_x_discrete(breaks=c("0", "0.4", "0.6", "0.8"), labels=c("None", "Low", "Medium", "High")) + 
  geom_hline(yintercept = 0, linetype = "longdash") +
  facet_wrap(~taulabel, labeller = label_parsed) + 
  theme(legend.position="right")
plotN750bw
