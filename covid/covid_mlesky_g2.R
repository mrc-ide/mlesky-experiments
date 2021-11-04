library( ape )
library( lubridate )
library( glue )
library( mlesky )
library( treedater ) 
require(ggplot2)
require(lubridate)
library(cowplot)
library(grid)
library(gridExtra)


LAGS = seq(-28, 28, by=7)
NCPU = 6

# load data 
tres = read.tree( 'Sample_England_2021-02-03_n_5000_nreps_10_n_tree_dating_10_dated.nwk' ) 
cdata = readRDS( 'covars_transformed_oxcgrt_loess_4pc.rds' )
cdata$time <- decimal_date( cdata$date )
cdata$sctime <- scale( cdata$time )
cdata$national_ContainmentHealthIndex <- scale( cdata$national_ContainmentHealthIndex )

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
	, model = 1
)
sg = mlskygrid( tre
	, tau = NULL, tau_lower = .001, tau_upper = 10
	, sampleTimes = sts
	, res = RES
	, ncpu = 12
	, NeStartTimeBeforePresent = max( sts ) -  decimal_date( as.Date( '2020-01-30' ))
	, model = 1
	, ncross = 5
)
TAU = sg$tau 

# function to run on all trees 
run_mlesky <- function(tre, subn = NULL, lag = 0, ktre = 0)
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
			, model = 1
		)
		if ( !is.null( subn )){
			saveRDS(sg, file = paste0('g2.', subn, '/', ktre, '.', 'nocovar', '.rds' ) )
		} else{
			saveRDS(sg, file = paste0('g2' , '/', ktre, '.', 'nocovar', '.rds' ) )
		}
		return(
			sg
		)
	}
	.cdata <- cdata 
	.cdata$time = decimal_date( cdata$date + lag )
	.cdata <- .cdata[, c('time', 'sctime', 'national_ContainmentHealthIndex') ]
	sg = mlskygrid( tre
		, tau = TAU
		, res = RES
		, sampleTimes = sts
		, NeStartTimeBeforePresent = max( sts ) -  decimal_date( as.Date( '2020-01-30' ))# min( cdata$time )
		, model = 1
		, formula_order = 1
		, formula = ~ national_ContainmentHealthIndex
		, data = .cdata
	)
	print( cbind( sg$time, sg$ne ))
	print( sg$beta )
	sg$data = .cdata 
	if ( !is.null( subn )){
		saveRDS(sg, file = paste0('g2.', subn, '/', ktre, '.', lag, '.rds' ) )
	} else{
		saveRDS(sg, file = paste0('g2' , '/', ktre, '.', lag, '.rds' ) )
	}
	sg 
}


PATH = 'g2' 
PATH150 = 'g2.150' 
PATH1000 = 'g2.1000' 

dir.create( PATH1000 )
dir.create( PATH150 )
dir.create( PATH )

# whole tree, no covariates
parallel::mclapply( seq_along(tres), function(k){
	sg = run_mlesky( tres[[k]], subn = NULL, lag = NA  , ktre = k  )
	0
}, mc.cores = NCPU )


# find max cor with lag 
fns = list.files( path= 'g2', patt = paste0('^.*\\.', 'nocovar', '\\.rds$'), full.name=TRUE )
rmat = sapply( fns, function(fn){
	f = readRDS( fn )
	ne = approx( f$time, f$ne , rule=1, xout = cdata$time)$y
	r = c( diff(log(ne)), NA)
	r
})
shift_var = sapply( LAGS, function(lag){
	.cdata = cdata 
	.cdata$time2 = decimal_date( cdata$date + lag )
	approx( .cdata$time2, .cdata$national_ContainmentHealthIndex, xout = .cdata$time )$y
})
keeprows = apply( rmat, 1, function(x) !all(is.na(x)))
rmat <- rmat[ keeprows, ] 
shift_var = shift_var[ keeprows, ] 
cormat  = cor( rmat, shift_var  )
medcor = apply( cormat, 2, function(x) median( na.omit( x )))
bestlag <- LAGS[ which.min( medcor ) ]


# subn 
parallel::mclapply( seq_along(tres), function(k){
	sg = run_mlesky( tres[[k]], subn = 150, lag = NA  , ktre = k  )
	0
}, mc.cores = NCPU )
## with lag 
parallel::mclapply( seq_along(tres), function(k){
	sg = run_mlesky( tres[[k]], subn = 150, lag = bestlag  , ktre = k  )
	0
}, mc.cores = NCPU )
###1000
parallel::mclapply( seq_along(tres), function(k){
	sg = run_mlesky( tres[[k]], subn = 1000, lag = NA  , ktre = k  )
	0
}, mc.cores = NCPU )
## with lag 
parallel::mclapply( seq_along(tres), function(k){
	sg = run_mlesky( tres[[k]], subn = 1000, lag = bestlag  , ktre = k  )
	0
}, mc.cores = NCPU )

## whole tree with lag 
parallel::mclapply( seq_along(tres), function(k){
	sg = run_mlesky( tres[[k]], subn = NULL, lag = bestlag  , ktre = k  )
	0
}, mc.cores = NCPU )






# figures

DAXIS <-  seq( as.Date('2020-01-30'), as.Date('2020-05-01'), by = 1)
TAXIS <- decimal_date( DAXIS )
LAGS = seq(-28, 28, by=7)


# cor plot 
cordf <- do.call( rbind, lapply( 1:ncol(cormat), function(i){
	data.frame( value = cormat[, i ] , lag = LAGS[i]  )
}))
quantiles_95 <- function(x) {
  r <- quantile(x , probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
corpl = ggplot(cordf, aes(x = as.factor(lag), y = value))  +
  theme_bw() + labs(y = "Cross-correlation containment health index") + geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text=element_text(size=16),       axis.title=element_text(size=24), strip.text = element_blank())+
  stat_summary(fun.data = quantiles_95, geom="boxplot", fill = "red", alpha = 0.4) + labs(x = "Time lag (days)")
corpl


# ne plots 
tree_info <- readRDS("sample_dates.rds")
datesdf <- data.frame(dates = as.Date(tree_info[[1]]$dates_all_trees))

fns = list.files( path= 'g2', patt = paste0('^.*\\.', bestlag, '\\.rds$'), full.name=TRUE )	
nemat = sapply( fns, function(fn){
	f = readRDS( fn )
	approx( f$time, f$ne , rule=1, xout = TAXIS)$y
})
qq = apply( nemat, 1, function(x) quantile( na.omit(x), c( .025, .5, .975 )))
dg2 <- data.frame( time = TAXIS, ymin = qq[1, ], ymax = qq[3, ], y  = qq[2, ], analysis = 'With covariate' )

fns = list.files( path= 'g2', patt = paste0('^.*\\.', 'nocovar', '\\.rds$'), full.name=TRUE )	
nemat = sapply( fns, function(fn){
	f = readRDS( fn )
	approx( f$time, f$ne , rule=1, xout = TAXIS)$y
})
qq = apply( nemat, 1, function(x) quantile( na.omit(x), c( .025, .5, .975 )))
dg2_nocovar <- data.frame( time = TAXIS, ymin = qq[1,], ymax = qq[3,], y  = qq[2,], analysis = 'Without' )
dne = rbind( dg2, dg2_nocovar )

pne = ggplot( aes(x = as.Date( date_decimal( time)), y = y,  ymin = ymin, ymax = ymax
	, col = as.factor(analysis), fill = as.factor(analysis)
	) 
	, data = dne
	) + 
  geom_path(size=1.5) +  
  labs(x='', col = "Analysis", fill = "Analysis") + 
  ylab('Effective population size' ) + 
  geom_path(data = dne , aes(x = as.Date( date_decimal( time)), y = ymax), alpha = 0.6, size = 1, linetype = "dashed")+
  geom_path(data = dne , aes(x = as.Date( date_decimal( time)), y = ymin), alpha = 0.6, size = 1, linetype = "dashed")+
  geom_ribbon( alpha = .25) +
  theme_minimal() + 
  theme(legend.position='top',panel.grid.minor = element_blank())+
  scale_x_date(date_breaks = "1 month", date_labels = '%b', 
               limits = c(min(datesdf$dates) - 7, max(datesdf$dates) + 7))+theme(axis.text=element_text(size=12),
                                                                                 axis.title=element_text(size=14))  + #scale_y_log10() +
  scale_color_manual(values = c("#1b9e77", "#7570b3"))+
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))+
  annotate("rect", xmin=as.Date("2020-03-23"), xmax=min(as.Date("2020-05-12"), max(datesdf$dates) + 7), ymin=0, ymax=Inf, alpha = 0.02512, fill = "red", col = "red")+
  theme(plot.title = element_text(hjust = 0.5, size=18), legend.position='top', legend.title =element_text(size=16),  legend.text = element_text(size=12), legend.key.size = unit(1, "cm")) 



## with n150 
fns = list.files( path= 'g2.150', patt = paste0('^.*\\.', bestlag, '\\.rds$'), full.name=TRUE )	
nemat = sapply( fns, function(fn){
	f = readRDS( fn )
	approx( f$time, f$ne , rule=1, xout = TAXIS)$y
})
qq = apply( nemat, 1, function(x) quantile( na.omit(x), c( .025, .5, .975 )))
dg2.150 <- data.frame( time = TAXIS, ymin = qq[1, ], ymax = qq[3, ], y  = qq[2, ], analysis = 'With covariate' )

fns = list.files( path= 'g2.150', patt = paste0('^.*\\.', 'nocovar', '\\.rds$'), full.name=TRUE )	
nemat = sapply( fns, function(fn){
	f = readRDS( fn )
	approx( f$time, f$ne , rule=1, xout = TAXIS)$y
})
qq = apply( nemat, 1, function(x) quantile( na.omit(x), c( .025, .5, .975 )))
dg2_nocovar.150 <- data.frame( time = TAXIS, ymin = qq[1,], ymax = qq[3,], y  = qq[2,], analysis = 'Without' )
dne150 = rbind( dg2.150, dg2_nocovar.150 )

pne150 = ggplot( aes(x = as.Date( date_decimal( time)), y = y,  ymin = ymin, ymax = ymax
	, col = as.factor(analysis), fill = as.factor(analysis)
	) 
	, data = dne150
	) + 
  geom_path(size=1.5) +  
  labs(x='', col = "Analysis", fill = "Analysis") + 
  ylab('Effective population size' ) + 
  geom_path(data = dne150 , aes(x = as.Date( date_decimal( time)), y = ymax), alpha = 0.6, size = 1, linetype = "dashed")+
  geom_path(data = dne150 , aes(x = as.Date( date_decimal( time)), y = ymin), alpha = 0.6, size = 1, linetype = "dashed")+
  geom_ribbon( alpha = .25) +
  theme_minimal() + 
  theme(legend.position='top',panel.grid.minor = element_blank())+
  scale_x_date(date_breaks = "1 month", date_labels = '%b', 
               limits = c(min(datesdf$dates) - 7, max(datesdf$dates) + 7))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))  +
  scale_color_manual(values = c("#1b9e77", "#7570b3"))+
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))+
  annotate("rect", xmin=as.Date("2020-03-23"), xmax=min(as.Date("2020-05-12"), max(datesdf$dates) + 7), ymin=0, ymax=Inf, alpha = 0.02512, fill = "red", col = "red")+
  theme(plot.title = element_text(hjust = 0.5, size=18), legend.position='top', legend.title =element_text(size=16),  legend.text = element_text(size=12), legend.key.size = unit(1, "cm"))  + 
  scale_y_log10() 


## with 1000 
fns = list.files( path= PATH1000, patt = paste0('^.*\\.', bestlag, '\\.rds$'), full.name=TRUE )	
nemat = sapply( fns, function(fn){
	f = readRDS( fn )
	approx( f$time, f$ne , rule=1, xout = TAXIS)$y
})
qq = apply( nemat, 1, function(x) quantile( na.omit(x), c( .025, .5, .975 )))
dg2.1000 <- data.frame( time = TAXIS, ymin = qq[1, ], ymax = qq[3, ], y  = qq[2, ], analysis = 'With covariate' )

fns = list.files( path= PATH1000, patt = paste0('^.*\\.', 'nocovar', '\\.rds$'), full.name=TRUE )	
nemat = sapply( fns, function(fn){
	f = readRDS( fn )
	approx( f$time, f$ne , rule=1, xout = TAXIS)$y
})
qq = apply( nemat, 1, function(x) quantile( na.omit(x), c( .025, .5, .975 )))
dg2_nocovar.1000 <- data.frame( time = TAXIS, ymin = qq[1,], ymax = qq[3,], y  = qq[2,], analysis = 'Without' )
dne1000 = rbind( dg2.1000, dg2_nocovar.1000 )

pne1000 = ggplot( aes(x = as.Date( date_decimal( time)), y = y,  ymin = ymin, ymax = ymax
	, col = as.factor(analysis), fill = as.factor(analysis)
	) 
	, data = dne1000
	) + 
  geom_path(size=1.5) +  
  labs(x='', col = "Analysis", fill = "Analysis") + 
  ylab('Effective population size' ) + 
  geom_path(data = dne1000 , aes(x = as.Date( date_decimal( time)), y = ymax), alpha = 0.6, size = 1, linetype = "dashed")+
  geom_path(data = dne1000 , aes(x = as.Date( date_decimal( time)), y = ymin), alpha = 0.6, size = 1, linetype = "dashed")+
  geom_ribbon( alpha = .25) +
  theme_minimal() + 
  theme(legend.position='top',panel.grid.minor = element_blank())+
  scale_x_date(date_breaks = "1 month", date_labels = '%b', 
               limits = c(min(datesdf$dates) - 7, max(datesdf$dates) + 7))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14))  +
  scale_color_manual(values = c("#1b9e77", "#7570b3"))+
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))+
  annotate("rect", xmin=as.Date("2020-03-23"), xmax=min(as.Date("2020-05-12"), max(datesdf$dates) + 7), ymin=0, ymax=Inf, alpha = 0.02512, fill = "red", col = "red")+
  theme(plot.title = element_text(hjust = 0.5, size=18), legend.position='top', legend.title =element_text(size=16),  legend.text = element_text(size=12), legend.key.size = unit(1, "cm"))  + 
  scale_y_log10() 
pne1000


# growth rate 
fns = list.files( path= 'g2', patt = paste0('^.*\\.', bestlag, '\\.rds$'), full.name=TRUE )	
nemat = sapply( fns, function(fn){
	f = readRDS( fn )
	ne = approx( f$time, f$ne , rule=1, xout = TAXIS)$y
	c( diff( log(ne))*7, NA ) 
})
qq = apply( nemat, 1, function(x) quantile( na.omit(x), c( .025, .5, .975 )))
dgr2 <- data.frame( time = TAXIS, ymin = qq[1, ], ymax = qq[3, ], y  = qq[2, ], analysis = 'With covariate' )

fns = list.files( path= 'g2', patt = paste0('^.*\\.', 'nocovar', '\\.rds$'), full.name=TRUE )	
nemat = sapply( fns, function(fn){
	f = readRDS( fn )
	ne = approx( f$time, f$ne , rule=1, xout = TAXIS)$y
	c( diff( log(ne))*7, NA ) 
})
qq = apply( nemat, 1, function(x) quantile( na.omit(x), c( .025, .5, .975 )))
dgr2_nocovar <- data.frame( time = TAXIS, ymin = qq[1,], ymax = qq[3,], y  = qq[2,], analysis = 'Without' )
dgr = rbind( dgr2, dgr2_nocovar )

pgr = ggplot( aes(x = as.Date( date_decimal( time)), y = y,  ymin = ymin, ymax = ymax, colour = as.factor(analysis), fill = as.factor(analysis)) 
	, data = dgr ) + 
  geom_path(size=1.5) +  labs(x='', col = "Analysis", fill = "Analysis") + ylab('Growth rate (1/week)' ) + 
  geom_ribbon( alpha = .38) +
  geom_path(data = dgr, aes(x = as.Date( date_decimal( time)), y = ymax), alpha = 0.6, size = .35)+
  geom_path(data = dgr, aes(x = as.Date( date_decimal( time)), y = ymin), alpha = 0.6, size = .35)+
  scale_x_date(date_breaks = "1 month", date_labels = '%b', 
	, limits = c(min(DAXIS) - 7, max(DAXIS) + 7))+
  theme_minimal() + #theme(legend.position='top',panel.grid.minor = element_blank())+
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14)) + 
  theme(plot.title = element_text(hjust = 0.5, size=18), legend.position='top', legend.title =element_text(size=16),  legend.text = element_text(size=12), legend.key.size = unit(1, "cm"))  + 
  geom_hline(yintercept = 0)+
  scale_color_manual(values = c("#1b9e77", "#7570b3"))+
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))+
  annotate("rect", xmin=as.Date("2020-03-23"), xmax=min(as.Date("2020-05-12"), max(datesdf$dates) + 7), ymin=-Inf, ymax=Inf, alpha = 0.02512, fill = "red", col = "red")


