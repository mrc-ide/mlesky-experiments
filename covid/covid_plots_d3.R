# Plotting mlesky analyses
# Required: 
# - RDS files mlesky analyses from covid_mlesky_d2.R

require(ggplot2)
require(lubridate)
library(cowplot)
library(grid)
library(gridExtra)



#' Function to extract key parameters from mlesky dfs
#' @param ofn Name of mlesky output RDS file to extract
extract_outcomes_mlesky <- function(ofn) {
  tN = readRDS(ofn)
  taxis_global = tN$time
  covar = "ContainmentHealthIndex"
  
  q_ne = t(apply( tN$ne, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ))
  gr =  apply( tN$ne, 2, function(x) c(exp( diff(log(x)) )^6.5, NA)  ) 
  q_gr = t( apply( gr, 1, function(x) quantile( na.omit(x), c(.5, .025, .975 )) ) )
  colnames( q_ne ) =  c( 'y', 'ylb', 'yub' )
  pldf0 = as.data.frame( q_ne ) ;  pldf0$time = tN$time ; pldf0$covar = covar; pldf0$group = group
  colnames( q_gr ) = c( 'y', 'ylb', 'yub' )
  gpldf0 = as.data.frame( q_gr ) ; gpldf0$time = tN$time ; gpldf0$covar = covar; gpldf0$group = group
  return(list(
    pldf0 = pldf0,
    gpldf0 = gpldf0,
    beta = data.frame(value = unname(tN$beta), covar = covar)
  ))
}


# Figure 5

# Read mlesky files
analysis_no_covar <- "covid/analysis_no_covar.rds"
analysis_covar_OxCGRT <- "covid/analysis_covar_OxCGRT.rds"

# Read sample dates
tree_info <- readRDS("covid/sample_dates.rds")
datesdf <- data.frame(dates = as.Date(tree_info[[1]]$dates_all_trees))


# Extract analysis no covar
out_no_covar <- lapply(analysis_no_covar, extract_outcomes_mlesky)

# Effective population size
Ne = lapply(out_no_covar, function(x) x$pldf0)
Ne = do.call(rbind, Ne)

# Reproduction number
Rt = lapply(out_no_covar, function(x) x$gpldf0)
Rt = do.call(rbind, Rt)

# Labelling for comparison
Ne_mlesky_england_no_covar =  Ne
Rt_mlesky_england_no_covar = Rt
Ne_mlesky_england_no_covar$analysis = "No covariate"
Rt_mlesky_england_no_covar$analysis = "No covariate"

# Extract analysis with covar
out_covar_OxCGRT <- lapply(analysis_covar_OxCGRT, extract_outcomes_mlesky)

# Effective population size
Ne = lapply(out_covar_OxCGRT, function(x) x$pldf0)
Ne = do.call(rbind, Ne)

# Reproduction number
Rt = lapply(out_covar_OxCGRT, function(x) x$gpldf0)
Rt = do.call(rbind, Rt)

# Labelling for comparison
Ne$analysis = "With covariate"
Rt$analysis = "With covariate"

# merging two dfs
Ne = rbind(Ne, Ne_mlesky_england_no_covar)
Rt = rbind(Rt, Rt_mlesky_england_no_covar)


# Ne through time
p0 = ggplot( aes(x = as.Date( date_decimal( time)), y = y,  ymin = ylb, ymax = yub, col = as.factor(analysis), fill = as.factor(analysis)) , data = Ne ) + 
  geom_path(size=1.5) +  labs(x='', col = "Analysis", fill = "Analysis") + ylab('Effective population size' ) + 
  geom_path(data = Ne , aes(x = as.Date( date_decimal( time)), y = yub), alpha = 0.6, size = 1, linetype = "dashed")+
  geom_path(data = Ne , aes(x = as.Date( date_decimal( time)), y = ylb), alpha = 0.6, size = 1, linetype = "dashed")+
  geom_ribbon( alpha = .25 ) +
  
  theme_minimal() + theme(legend.position='top',panel.grid.minor = element_blank())+
  scale_x_date(date_breaks = "1 month", date_labels = '%b', 
               limits = c(min(datesdf$dates) - 7, max(datesdf$dates) + 7))+theme(axis.text=element_text(size=12),
                                                                                 axis.title=element_text(size=14))  + #scale_y_log10() +
  scale_color_manual(values = c("#1b9e77", "#7570b3"))+
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))+
  annotate("rect", xmin=as.Date("2020-03-23"), xmax=min(as.Date("2020-05-12"), max(datesdf$dates) + 7), ymin=0, ymax=Inf, alpha = 0.02512, fill = "red", col = "red")+
  ggtitle("Epidemic trajectory in England (Spring 2020)")+
  theme(plot.title = element_text(hjust = 0.5, size=18), legend.position='top', legend.title =element_text(size=16),  legend.text = element_text(size=12), legend.key.size = unit(1, "cm"))

p0


# Rt through time
p1 = ggplot( aes(x = as.Date( date_decimal( time)), y = y,  ymin = ylb, ymax = yub, col = as.factor(analysis), fill = as.factor(analysis) ) , data = Rt ) + 
  geom_path(size=1.5) +  labs(x='', col = "Analysis", fill = "Analysis") + ylab('Reproduction number' ) + 
  geom_path(data = Rt , aes(x = as.Date( date_decimal( time)), y = yub), alpha = 0.6, size = 1, linetype = "dashed")+
  geom_path(data = Rt , aes(x = as.Date( date_decimal( time)), y = ylb ), alpha = 0.6, size = 1, linetype = "dashed")+
  geom_ribbon( alpha = .25 ) +
  
  theme_minimal() + theme(legend.position='',panel.grid.minor = element_blank())+
  scale_x_date(date_breaks = "1 month", date_labels = '%b', 
               limits = c(min(datesdf$dates) - 7, max(datesdf$dates) + 7))+theme(axis.text=element_text(size=12),
                                                                                 axis.title=element_text(size=14)) + scale_y_log10()+
  scale_color_manual(values = c("#1b9e77", "#7570b3"))+
  scale_fill_manual(values = c("#1b9e77", "#7570b3"))+
  
  # geom_vline(xintercept = as.Date("2020-06-30"))+ geom_vline(xintercept = as.Date("2020-08-01")) + 
  geom_hline(yintercept = 1)+
  annotate("rect", xmin=as.Date("2020-03-23"), xmax=min(as.Date("2020-05-12"), max(datesdf$dates) + 7), ymin=0, ymax=Inf, alpha = 0.02512, fill = "red", col = "red")
p1
gc()




# Reading in confirmed case data
confirmed_cases <- read.csv("cases_england.csv")
confirmed_cases = confirmed_cases[confirmed_cases$areaName == "England",]
confirmed_cases$time = decimal_date(as.Date(as.character(confirmed_cases$date)))
confirmed_cases = confirmed_cases[confirmed_cases$time > min(Ne_skygrowth_england$time) & confirmed_cases$time <( max(Ne_skygrowth_england$time)+0.1), ]
confirmed_cases$y = scales:::rescale(loess(confirmed_cases$newCasesBySpecimenDate ~ confirmed_cases$time,degree=1,  span = 0.2)$fitted, to = c(min(Ne_skygrowth_england$ylb), max(Ne_skygrowth_england$yub)))
confirmed_cases$ylb = scales:::rescale(loess(confirmed_cases$newCasesBySpecimenDate ~ confirmed_cases$time,degree=1,  span = 0.2)$fitted, to = c(min(Ne_skygrowth_england$ylb), max(Ne_skygrowth_england$yub)))
confirmed_cases$yub = scales:::rescale(loess(confirmed_cases$newCasesBySpecimenDate ~ confirmed_cases$time,degree=1,  span = 0.2)$fitted, to = c(min(Ne_skygrowth_england$ylb), max(Ne_skygrowth_england$yub)))
confirmed_cases$analysis = NA


# adding confirmed cases to Ne plot
p0 = p0 +  geom_line(data = confirmed_cases, aes(x = as.Date( date_decimal( time)), y = y), col = "black", size = 0.8,alpha = 0.6,linetype = "dashed")

# plotting all plots together
P0 = cowplot::plot_grid(p0, p1, p2, align = "h", nrow =3, rel_heights = c(1, 0.9, 0.6 ))

P0



# saving plot
ggsave( plot = P0, file = paste0('covid/England_national_lockdown.pdf'), width = 12, height = 12 )




# Figure 6

# Extract mlesky file
mlesky_files = list.files()[grep("Containment", list.files())]
out <- lapply(mlesky_files, extract_outcomes_mlesky, compare = "lineage", lag = T)


# extract beta
beta = lapply(out, function(x) x$beta)
beta = lapply(beta, function(x) {
  y = t(x)
  rownames(y) <- c()
  return(data.frame(y, lineage = y[nrow(y)]))
})

for(i in 1:length(mlesky_files)) {
  if(length(grep("lag", mlesky_files[i])) == 0) {
    lag = 0
  } else {
    lag =  paste0(strsplit(strsplit(strsplit(strsplit(mlesky_files[i], "/")[[1]][length(strsplit(mlesky_files[i], "/")[[1]])], "2021")[[1]][2], "lag_")[[1]][2], "number")[[1]][1]  )
  }
  beta[[i]]$lag = as.numeric(lag)
}

# merge analyses and make into long format
betadf_Containment = as.data.frame(do.call(rbind, beta))
betadf_Containment$national_ContainmentHealthIndex = as.numeric(as.character(betadf_Containment$national_ContainmentHealthIndex ))
betadf_Containment$index = "Containment"
names(betadf_Containment) = c("value", "lineage", "lag", "index")
betadf = betadf_Containment
betadf_melted = reshape2::melt(betadf, id.vars = c("lineage", "index", "value", "lag"))


# plot
bpl = ggplot(betadf_melted, aes(x = as.factor(lag), y = value)) + facet_wrap(~lineage, nrow = 1) + #geom_boxplot() + 
  theme_bw() + labs(y = "beta") + geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text=element_text(size=16),       axis.title=element_text(size=24), strip.text = element_blank())+
  stat_summary(fun.data = quantiles_95, geom="boxplot", fill = "grey", alpha = 0.4) + labs(x = "Time delay (days)")
bpl

# saving plot
ggsave(bpl, filename = paste0("covid_covar.pdf"), width = 12)
