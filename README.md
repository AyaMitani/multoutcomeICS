# multoutcomeICS
R codes for simulation study in "Marginal analysis of multiple outcomes with informative cluster size"

## simulation.R
R code for running simulation that includes
1. Functions
2. Code for data generation
3. Fitting simulated data with<br/>
      i. univariate GEE and CWGEE for each outcome<br/>
      ii. multivariate GEE and CWGEE with common slope for MetS<br/>
           - with unstructured, exchangeable, independent correlaiton structures<br/>
      iii. multivariate GEE and CWGEE with unique slope for MetS<br/>
             - with unstructured, exchangeable, independent correlaiton structures
4. Save results 

## create_table_simresults.R
R code for creating tables for simulation results

## create_splitdensity_plot.R
R code for creating split density plots for simulation results in
Referenced from: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2

## create_plot_for_type1error_mvWald.R
R code for creating Type I error plots for simulation results
