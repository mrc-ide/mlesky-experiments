library(ape)
library(lubridate)
library(glue)
library(mlesky)
library(treedater) 
require(ggplot2)
require(lubridate)
library(cowplot)
library(grid)
library(gridExtra)


LAGS = seq(-28, 28, by=7)
NCPU = 6
MODEL = 2 #skygrid 

# load data 
tres = read.tree( 'Sample_England_2021-02-03_n_5000_nreps_10_n_tree_dating_10_dated.nwk' ) 
cdata = readRDS( 'covars_transformed_oxcgrt_loess_4pc.rds' )
cdata$time <- decimal_date( cdata$date )
cdata$sctime <- scale( cdata$time )
cdata$intercept<- 1
# scale and smooth covariate 
cdata$national_ContainmentHealthIndex <- scale( cdata$national_ContainmentHealthIndex )
cdata$smooth_national_ContainmentHealthIndex <- predict( loess( national_ContainmentHealthIndex ~ time, data = cdata )  )

# optimise res and tau based on first tree
tre = tres[[1]] 
sts <- as.numeric( sapply( strsplit( tre$tip.label, split = '\\|' ), '[', 3 ) )
tre$tip.label <- sapply( strsplit(tre$tip.label, '[|]') , '[' , 1 )
names( sts ) <- tre$tip.label 
RES=optim_res_aic(tre
	, res = round( (max( sts ) -  decimal_date( as.Date( '2020-01-30' )) )*365 / c(7,14,21,28) ) 
	, ncpu = 12
	, tau = NULL, tau_lower = .001, tau_upper = 10
	, sampleTimes = sts
	, NeStartTimeBeforePresent = max( sts ) -  decimal_date( as.Date( '2020-01-30' ))
	, model = MODEL
)
sg = mlskygrid( tre
	, tau = NULL, tau_lower = .001, tau_upper = 10
	, sampleTimes = sts
	, res = RES
	, ncpu = NCPU
	, NeStartTimeBeforePresent = max( sts ) -  decimal_date( as.Date( '2020-01-30' ))
	, model = MODEL
	, ncross = 5
)
TAU = sg$tau

sgc = mlskygrid( tre
	, tau = TAU
	, res = RES
	, sampleTimes = sts
	, NeStartTimeBeforePresent = max( sts ) - decimal_date( as.Date( '2020-01-30' ))# min( cdata$time )
	, model = MODEL
	, data = cdata
	, formula = diffLogNe ~ national_ContainmentHealthIndex + 1
)



#' function to run on all trees 
#' @param subn subsample size; if NULL use entire tree
#' @param lag Days to shift covariate
#' @param ktre id for output files 
run_mlesky <- function(tre, subn = NULL, lag = 0, ktre = 0, fo = diffLogNe ~ national_ContainmentHealthIndex + 1, path = NULL )
{
	if ( !is.null( subn )){
		keep <- sample( tre$tip.label, replace=FALSE, size = subn )
		tre = keep.tip( tre, keep )
	}
	sts <- as.numeric( sapply( strsplit( tre$tip.label, split = '\\|' ), '[', 3 ) )
	tre$tip.label <- sapply( strsplit(tre$tip.label, '[|]') , '[' , 1 )
	names( sts ) <- tre$tip.label 
	if ( is.na( lag ) ){
		sg = mlskygrid( tre
			, tau = TAU
			, res = RES
			, sampleTimes = sts
			, NeStartTimeBeforePresent = max( sts ) -  decimal_date( as.Date( '2020-01-30' ))# min( cdata$time )
			, model = MODEL
		)
		if ( !is.null( subn )){
			saveRDS(sg, file = paste0('g3.', subn, '/', ktre, '.', 'nocovar', '.rds' ) )
		} else{
			saveRDS(sg, file = paste0('g3' , '/', ktre, '.', 'nocovar', '.rds' ) )
		}
		return(
			sg
		)
	}
	.cdata <- cdata 
	.cdata$time = decimal_date( cdata$date + lag )
	.cdata <- .cdata[, c('time', 'sctime', 'national_ContainmentHealthIndex', 'smooth_national_ContainmentHealthIndex', 'intercept') ]
	sg = mlskygrid( tre
		, tau = TAU
		, res = RES
		, sampleTimes = sts
		, NeStartTimeBeforePresent = max( sts ) - decimal_date( as.Date( '2020-01-30' ))# min( cdata$time )
		, model = MODEL
		, formula = fo
		, data = .cdata
	)
	print( cbind( sg$time, sg$ne ))
	print( sg$beta )
	sg$data = .cdata 
	if ( is.null( path )){
		if ( !is.null( subn )){
			saveRDS(sg, file = paste0('g3.', subn, '/', ktre, '.', lag, '.rds' ) )
		} else{
			saveRDS(sg, file = paste0('g3' , '/', ktre, '.', lag, '.rds' ) )
		}
	} else{
		if ( !is.null( subn )){
			saveRDS(sg, file = paste0(path, '/', ktre, '.', lag, '.rds' ) )
		} else{
			saveRDS(sg, file = paste0(path, '/', ktre, '.', lag, '.rds' ) )
		}
	}
	sg 
}


PATH = 'g3' 
PATH150 = 'g3.150' 
PATH.nointercept = 'g3.nointercept'

dir.create( PATH150 )
dir.create( PATH )
dir.create( PATH.nointercept )



# whole tree, no covar & default fo, multiple lags 
for (lag in LAGS)
{
	parallel::mclapply( seq_along(tres), function(k){
		sg = run_mlesky( tres[[k]], subn = NULL, lag = lag  , ktre = k  )
		0
	}, mc.cores = NCPU )
}

# subsample, no covar & default fo, multiple lags 
for (lag in LAGS)
{
	parallel::mclapply( seq_along(tres), function(k){
		sg = run_mlesky( tres[[k]], subn = 150, lag = lag, ktre = k  )
		0
	}, mc.cores = NCPU )
}


# whole tree, no covar & no intercept fo, multiple lags 
for (lag in LAGS)
{
	parallel::mclapply( seq_along(tres), function(k){
		sg = run_mlesky( tres[[k]], subn = NULL, lag = lag, fo = diffLogNe ~ national_ContainmentHealthIndex-1  , ktre = k  )
		0
	}, mc.cores = NCPU )
}



# figures

DAXIS <-  seq( as.Date('2020-01-30'), as.Date('2020-05-01'), by = 1)
TAXIS <- decimal_date( DAXIS )

quantiles_95 <- function(x) {
  r <- quantile(x , probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
