library( ape ) 
library( mlesky )
library( lubridate )
library( glue ) 
ncpu = 8
rhcolour = '#DD6200' 

md <- read.csv( 'Google_mobility-UK12June2022.csv', stringsAs=FALSE )
md$date <- dmy( md$date ) 

trees <- readRDS( 'b117-dated_trees.rds' )


.fitmlesky <- function( lag = 0  )
{
	print ( Sys.time() )
	print( lag )
	print( '-----------') 
	
	md$time <- decimal_date( md$date + lag )
	# cubic smoothing spline with 5x cv 
	md$y <- smooth.spline( md$time, md$transit_stations_percent_change_from_baseline,  cv = 5)$y
	md$y <- scale( md$y )
	
	
	fit = mlskygrid( trees[[1]], sampleTimes = trees[[1]]$sts 
		, res = NULL, tau = NULL, tau_lower = 1e-5, tau_upper = 10
		, NeStartTimeBeforePresent =  max( trees[[1]]$sts ) - decimal_date(as.Date('2020-11-01'))
		, ncpu = ncpu
		, model = 2 , formula = diffLogNe ~ y , data = md)
	taxis <- fit$time 
	res = pbmcapply::pbmclapply( 2:length(trees), function(irep){
		sts <- trees[[irep]]$sts 
		f1 <- mlskygrid( trees[[irep]], sampleTimes = sts, res = fit$res,
					 tau = fit$tau, tau_tol = fit$tau_tol , ncross = fit$ncross, quiet = fit$quiet
					 , NeStartTimeBeforePresent =  max( trees[[irep]]$sts ) - decimal_date(as.Date('2020-11-01'))
					 , ne0 = median( fit$ne ), adapt_time_axis = FALSE, formula = fit$formula,
					 data = fit$data, ncpu = ncpu, model = fit$model )
		af <- approxfun( f1$time, f1$ne, rule = 2)
		afgr <- approxfun( f1$time, f1$growthrate , rule =2)
		list(ne = af(taxis), beta = f1$beta, growthrate = afgr(taxis) )
	}, mc.cores = ncpu )
	
	nemat <- do.call( cbind, lapply( res, '[[', 'ne' ) )
	grmat <- do.call( cbind, lapply( res, '[[', 'growthrate' ))
	
	fit$ne_ci <- cbind( 
		nelb = unname(  apply( nemat, 1, FUN = function(x) quantile(x, c(.025))) )
		, ne = unname(  apply( nemat, 1, FUN = function(x) quantile(x, c(.50))) )
		, neub = unname(  apply( nemat, 1, FUN = function(x) quantile(x, c(.975))) )
	)
	
	fit$growthrate_ci <- cbind( 
		grlb = unname(  apply( grmat, 1, FUN = function(x) quantile(x, c(.025))) )
		, gr = unname(  apply( grmat, 1, FUN = function(x) quantile(x, c(.50))) )
		, grub = unname(  apply( grmat, 1, FUN = function(x) quantile(x, c(.975))) )
	)
		
	if ( !is.null( fit$beta ))
	{
		betamat <- do.call( cbind, lapply( res, '[[', 'beta' ) )
		fit$beta_ci <- cbind( 
			betalb= apply( betamat, MARGIN = 1 , FUN=function(x) quantile(x, prob=.025 ) )
			 , beta = apply( betamat, MARGIN = 1 , FUN=function(x) quantile(x, prob=.50 ) )
			 , betaub = apply( betamat, MARGIN = 1 , FUN=function(x) quantile(x, prob=.975 ) )
			)
	}

	saveRDS( fit, file = glue('b117.1-fit-lag{lag}.rds' ) )
	with( fit$covar.df, matplot( time, cbind(scale(y),scale(diffLogNe)) ))
	print( fit$beta_ci )
}


lags <- seq(-15, 36, by = 3 )
if (TRUE)
{
	for ( lag in lags )
		.fitmlesky( lag )
}

betas <- sapply( lags , function( lag ){
	f = readRDS( glue('b117.1-fit-lag{lag}.rds' ) ) 
	f$beta_ci['y', ] 
})

library( ggplot2 ) 
lagdf = as.data.frame( cbind( lag = lags, t( betas )) )
pcoef = ggplot(data = lagdf) + geom_errorbar( aes(x = lag, ymin = betalb, ymax = betaub) ) + 
	geom_point( aes( x=lag, y = beta ) ) + 
	xlab('Lag') + 
	ylab(' \u0394 log(Ne) / mobility ') +
	geom_hline( yintercept = 0, lty = 2) + 
	theme_classic() 

ggsave( pcoef, file = glue('b117.1-pcoef.png') , width = 3.85, height = 2.8)
ggsave( pcoef, file = glue('b117.1-pcoef.svg') , width = 3.85, height = 2.8)
ggsave( pcoef, file = glue('b117.1-pcoef.pdf') , width = 3.85, height = 2.8)



f = readRDS( 'b117.1-fit-lag21.rds'  ) 
pne = plot(f, logy=F, ggplot=T) + theme_classic() + xlab('')
pne$data$date <- as.Date( date_decimal( pne$data$t ))
(pne1 <- ggplot( data = pne$data , aes(x = date) ) + 
	geom_ribbon(aes(ymin = nelb, ymax = neub), alpha = .50, fill = rhcolour, colour='white') + 
	geom_line( aes(y =  nemed), lwd = .8) + 
	theme_classic() + xlab('') + ylab('Effective population size' )
)
ggsave( pne1, file = glue('b117.1-pne.png') , width = 3.5, height = 2.8)
ggsave( pne1, file = glue('b117.1-pne.pdf') , width = 3.5, height = 2.8)



# growth rate figure 

md1 <- md[ md$date >= as.Date( '2020-10-31' ) & md$date < as.Date( '2021-02-01' ) , ]
md1$time <- decimal_date( md1$date )
# cubic smoothing spline with 5x cv 
md1$sctransit <- scale( md1$transit_stations_percent_change_from_baseline )
md1$y <- smooth.spline( md1$time, md1$sctransit,  cv = 5)$y
md1$grlb <- approx( f$time, f$growthrate_ci[,'grlb'] , xout = md1$time, rule = 2)$y
md1$grub <- approx( f$time, f$growthrate_ci[,'grub'] , xout = md1$time, rule = 2)$y
md1$gr <- approx( f$time, f$growthrate_ci[,'gr'] , xout = md1$time, rule = 2)$y


meantransit <- mean( md1$transit_stations_percent_change_from_baseline )
sdtransit <- sd( md1$transit_stations_percent_change_from_baseline )
meangr <- mean( md1$gr )
sdgr <- sd( md1$gr )


pgr <- ggplot(data=md1, aes(x = date ) ) + 
  geom_line( aes(y = (y)*sdgr+meangr ), lty =2 , colour = rhcolour) + 
  geom_point( aes(y = (sctransit)*sdgr+meangr ) , colour = rhcolour) + 
  geom_line( aes(y = gr ), lty =1 ) + 
  geom_line( aes(y = grub ), lty = 3) + 
  geom_line( aes(y = grlb ), lty = 3) + 
  scale_y_continuous( 
	name = 'Growth rate Ne (1/year)' 
	, sec.axis = sec_axis(trans=~(.+meangr)*sdgr, name = 'Transit mobility score') 
  ) + 
  theme_classic() +
  theme( axis.title.y.right = element_text(color = rhcolour) )+  xlab('' )

ggsave( pgr, file = glue('b117.1-pgr.png') , width = 3.5, height = 2.8)
ggsave( pgr, file = glue('b117.1-pgr.pdf') , width = 3.5, height = 2.8)
