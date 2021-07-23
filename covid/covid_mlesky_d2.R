# Maximum likeihood estimation of Ne(t) & R(t) from trees dated with treedater
# Required: 
# - RDS files containing dated trees from covid_dating_trees_d1.R

# devtools::install_github('emvolz-phylodynamics/mlesky')


library( ape )
library( lubridate )
library( glue )
library( mlesky )
library( treedater ) 
library( sarscov2 ) 
library( ggplot2 )
library( grid )
library( gridExtra )
library( ggtree )
library( alakazam )
library( stringi )


#' Runs mlesky on dated trees from treedater
#' @param tds_list A list of a list of dated trees. Output from d1_date_trees.R
#' @param ofn Name of mlesky output RDS file to save
#' @param taxis Time period over which to plot mlesky
#' @return Output from mlesky: Effective population size over time (Ne(t))
run_mlesky_covar <- function(tds_list_fn, ofn, taxis, formula, data, formula_order, mc.cores) {
  
  tds_list = readRDS(tds_list_fn)
  
  res_mlesky_list = lapply(tds_list, function(tds) {
    
    weeks = round(as.numeric((date_decimal(max(tds[[1]]$sts))-date_decimal(min(tds[[1]]$sts)))/7))

    class( tds ) <- 'multiPhylo' 
    
    tds = lapply( tds , 
                  function(x) {x$tip.label = unlist(lapply(strsplit(x$tip.label, '[|]'), function(y) paste0(y[1])))
                  return(x)} )
    
    TimeBeforePresent =  max( tds[[1]]$sts ) - min( tds[[1]]$sts )
    
    sgs = parallel::mclapply( tds, function(td) {
      mlskygrid(td, tau = NULL, tau_lower=.001, tau_upper = 10 , sampleTimes = td$sts , res = res, ncpu = 3, NeStartTimeBeforePresent = NeStartTimeBeforePresent,
                formula_order = formula_order,
                model = 1,
                formula = formula,
                data = data
                
      )
    }, mc.cores = mc.cores)
    
    
    out = lapply(sgs, function(sg) {
      with( sg, approx( time, ne, rule = 1, xout = taxis )$y )
    })
    
    return(list(out = out, beta = lapply(sgs, function(sg) sg$beta)))
  })
  
  
  res_mlesky <- 
    do.call( cbind, lapply(res_mlesky_list, function(x) do.call( cbind, x$out ) ))
  
  saveRDS( list( time = taxis, ne = res_mlesky, beta = lapply(res_mlesky_list, function(x) x$beta))  , file=paste0(ofn, "_mlesky", '.rds' ))
  
  res_mlesky
}

# read in OxCGRT data
covars_data <- readRDS("covars_transformed_oxcgrt_loess_4pc.rds")


# Run mlesky; looping for different back-shifts
for(lag in seq(-16, 28, 2)) {
  covars_data$time = decimal_date( as.Date(lubridate::ymd(covars_data$date) + lubridate::days(lag)))
  run_mlesky = run_mlesky_covar(tds_list_fn = tds_list_fn, formula = ~ national_ContainmentHealthIndex, data = covars_data, formula_order = 1, mc.cores = 7,
                                ofn = paste0('Sample_England_', Sys.Date(), '_', "_lineage_", lineage, "_covar_", "national_ContainmentHealthIndex", "_formula_order1_lag_", lag, "number_res_", "2perweek_loess4pc"), taxis = taxis)
}


# Run mlesky no covar
run_mlesky_no_covar = run_mlesky_covar(tds_list_fn = tds_list_fn, formula = ~ national_ContainmentHealthIndex, data = NULL, formula_order = 1, mc.cores = 7,
                              ofn = paste0('analysis_no_covar'), taxis = taxis)
