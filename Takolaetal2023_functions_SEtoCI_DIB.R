# Author: Elina Takola
# Date: 04 Oct 2023
# Publication: Takola et al. (in prep.) An open-access global database of meta-analyses 
# investigating yield and biodiversity responses to different management practices.

# All datasets are available in https://osf.io/7s9r4/ 

# This script contains the functions that were used to convert standard errors into 95% confidence intervals

# Functions for upper and lower CI
conf_interval_upper <- function(mean, n, se) {
  t_crit <- qt(0.95, n-1)
  upper_bound <- mean + t_crit * se
  return(upper_bound)
}


conf_interval_lower <- function(mean, n, se) {
  t_crit <- qt(0.95, n-1)
  lower_bound <- mean - t_crit * se
  return(lower_bound)
}


#############################
#> sessionInfo()
#R version 4.3.0 (2023-04-21 ucrt)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19045)
#Matrix products: default
#locale:
#  [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
#[4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
#time zone: Europe/Berlin
#tzcode source: internal
#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#loaded via a namespace (and not attached):
#  [1] compiler_4.3.0 tools_4.3.0    tinytex_0.45   xfun_0.39 
#############################
