'
mlesky applied to NC dater trees

findings Ne 
* signif correlated with v = i / p^2 , coef ~ 1 
* signif correlated with p, __but in the wrong direction__ 
* not correlated with incidence 

'

library( ape ) 
library( mlesky ) 
library( lubridate )


dtre <- read.tree( 'nc_refs_b.dater' )

md <- read.csv( 'nc_refs_b.dates.csv', stringsAs=FALSE )
sts <- setNames(  md$date[ match( dtre$tip.label, md$name )],  dtre$tip.label )



dsg2 <- mlskygrid( dtre, sampleTimes = sts,  tau = NULL, tau_lower = 1e-3, tau_upper = 1e2 , model=2 , ncpu = 8) 

st0 = Sys.time() 
dsg2 <- parboot( dsg2, nrep = 100, ncpu = 4 )
st1 = Sys.time() 



plot( dsg2 ) 

# covariates 

i <- read.csv( 'incidence.csv', stringsAs=FALSE , skip = 4)
p <- read.csv( 'prevalence.csv', stringsAs=FALSE , skip = 4)

i$i <- readr::parse_number(i$Cases ) 
i$p <- readr::parse_number( p$Cases ) 
i$v <- with( i, i / p^2 )
i$time <- i$Year + 0.5

i$logv <- log(i$v )
i$logp <- log(i$p )
i$logi <- log(i$i )


dsg2v <- mlskygrid( dtre
 , tau = NULL
 , tau_lower = 1e-3
 , tau_upper = 1e2  
 , model=2 
 , ncpu = 8
 , sampleTimes = sts
 , data = i
 , formula = logNe ~ logv
 , res = NULL 
) 
dsg2v <- parboot( dsg2v , nrep = 100 , ncpu = 4 )
print( dsg2v$beta_ci )
'
> dsg2c$beta_ci
                 betalb      beta     betaub
(Intercept) -90.3187629 17.953544 -23.649749
logv          0.6411651  1.497021   2.410505
'

# incidence 
dsg2i <- mlskygrid( dtre
 , tau = 13.36#NULL
 , tau_lower = 1e-3
 , tau_upper = 1e2  
 , model=2 
 , ncpu = 8
 , sampleTimes = sts
 , data = i
 , formula = logNe ~ logi
 , res = NULL 
) 
dsg2i <- parboot( dsg2i , ncpu = 12 )
print( dsg2i$beta_ci )
'
                 betalb      beta     betaub
(Intercept) 45.55197806 -6.195265 246.820865
logi        -0.08771629  1.801982   8.780722
'

# prevalence 
dsg2p <- mlskygrid( dtre
 , tau = 13.36#NULL
 , tau_lower = 1e-3
 , tau_upper = 1e2  
 , model=2 
 , ncpu = 8
 , sampleTimes = sts
 , data = i
 , formula = logNe ~ logp
 , res = NULL 
) 
dsg2p <- parboot( dsg2p , ncpu = 12 )
print( dsg2p$beta_ci )
'
                betalb      beta     betaub
(Intercept) -188.65880 31.550449 -64.378371
logp          -7.19952 -5.864419  -1.460562
'

library( ggplot2 ) 
pldf = data.frame( 
	ne = dsg2v$ne_ci[,2] 
	, nelb = dsg2v$ne_ci[,1] 
	, neub = dsg2v$ne_ci[,3] 
	, time = dsg2v$time 
) 
pldf$date <- as.Date( date_decimal( pldf$time ) )
pldf <- pldf[ pldf$time > 2005 , ] 
pldf$`Cases per year` <- exp( approx ( i$time, i$logi, xout = pldf$time, rule=1 )$y)
pldf$PLWHIV <- exp( approx ( i$time, i$logp, xout = pldf$time, rule=1 )$y)
pldf$nu <- 1/exp( approx ( i$time, i$logv, xout = pldf$time, rule=1 )$y )
#~ U+1D708
#~ `\u1d708 (t)`

p0 <- ggplot( data = pldf ) + 
  geom_path( aes(date, ne ), colour = 'black' , lwd = 1.4) + 
  geom_ribbon( aes( date, ymin = nelb, ymax =  neub ) , alpha = .25, fill = 'black', colour = 'white' ) + 
  xlab('') + ylab('Effective population size') + theme_classic()
  
p1 <- ggplot( data = pldf )  + geom_path(aes(date, scale(nu,scale=T)) , colour = '#DC9917')  + 
  geom_path( aes(date, scale( `Cases per year`, scale=T)), colour = '#008000' ) + 
  geom_path( aes(date, scale(PLWHIV,scale=T)), colour = '#A52A2A') + 
  xlab('') + ylab('') + theme_classic()
  
