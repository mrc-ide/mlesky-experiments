library(devtools)
# Install mlesky updated version
#devtools::install_github("emvolz-phylodynamics/mlesky")
library(mlesky)
library(ape)
#install.packages("pbmcapply")
library(pbmcapply)
library(dplyr)
library(ggplot2)
library(ggpubr)
#library(data.table)
library(reshape2)
#install.packages("Metrics")
#library(Metrics)
library(stringr)
library(ggforce)
library(glue)

ncpu = 7
nsim = 500
scaleFUN = function(x) sprintf("%.2f", x) # format to 2 decimal places
scale4FUN = function(x) sprintf("%.4f", x) # format to 4 decimal places

set.seed(3)
#set.seed(425672)
sampleDates_200 =seq(-20,0,0.1) #2000,2020,0.1
sampleDates_100 =seq(-20,0,0.2) #2000,2020,0.2

# Simulate coalescent trees using the simCoal method implemented in mlesky
simulate_ntrees <- function(sampDates, alphaFun) {
  tres = pbmclapply(1:nsim, function(i) {
    simCoal(sampDates, alphaFun, alphaMin=0.1)
  }, mc.cores = ncpu)
  return(tres)
}

# Get Ne estimates on common time axis and plot coverage probabilities over time
# Return a list contaning 8 items:
# (i) pboot_ne_ind: Ne estimates for each tree without getting common time axis
# (ii) pboot_ne_adj_df: Ne estimates using common time axis
# (iii) cov_prob: all Ne estimates and true value and its coverages on dataframe
# (iv) cov_prob_matrix: coverage probability matrix (rows = time axis elements; columns = simulation indexes)
# (v) cov_prob_avg: coverage probability over entire time axis
# (vi) cov_prob_over_time: coverage probability for each common time axis point over different simulations (rowMeans output)
# (vii) cov_prob_over_time_df: same as (vi) but on data.frame format to allow easy plotting using ggplot2
# (viii) common_time_ax: common time axis to use for plots
# (ix) mean_abs_error_df: Mean absolute error (MAE) for each simulation index
# (x) rmse_df: Root mean squared error (RMSE) for each simulation index
# Also plot the coverage probability over time to `coverage_plots` folder that is appended to the supplied `out_path`
# and print most of these returned values on the screen in the end of execution
get_nsim_estimates <- function(sim_trees, alphaFun, model, res_choice, out_path_cov_prob, out_path_error, out_path_rmse) { #sampDates
  print("Fit")
  #sts = sampDates; names(sts) = sim_trees[[1]]$tip.label
  fit = pbmclapply(X=sim_trees, FUN=mlskygrid, res=res_choice,tau=NULL,tau_lower = 0.001,tau_upper = 1000,model = model,ncpu=ncpu) #sampleTimes=sts
  pboot = pboot_ne = list()
  for(i in 1:length(sim_trees)){
    print(paste("===",i,"==="))
    if(length(fit[[i]]$time) <= 2) {
      next
    }
    pboot[[i]] = parboot(fit[[i]], nrep=200, ncpu=ncpu)
    pboot_ne[[i]] = pboot[[i]]["ne_ci"]
    pboot_ne[[i]]$nelb = pboot[[i]]$ne_ci[,1]
    pboot_ne[[i]]$ne = pboot[[i]]$ne_ci[,2]
    pboot_ne[[i]]$neub = pboot[[i]]$ne_ci[,3]
    pboot_ne[[i]]$time = pboot[[i]]$time
    pboot_ne[[i]]$time = as.data.frame(pboot_ne[[i]]$time)
    #pboot_ne[[i]]$true_ne = alphaFun(pboot_ne[[i]]$time); #pboot_ne[[i]]$true_ne = as.data.frame(pboot_ne[[i]]$true_ne)
    pboot_ne[[i]] = cbind(pboot_ne[[i]]$time, pboot_ne[[i]]$ne_ci) #pboot_ne[[i]]$true_ne
    colnames(pboot_ne[[i]]) = c("time","nelb","est_ne","neub") #"true_ne"
  }
  # Get estimates for each tree without getting common time axis
  pboot_ne_ind = bind_rows(pboot_ne, .id = 'n_sim')
  print("Estimates without common time axis:")
  print(pboot_ne_ind)

  # Get first most recent coalescent time and last less recent coalescent time among all simulated trees to plot out Ne values on common time axis
  first_coal_sims = pboot_ne_ind %>% group_by(n_sim) %>% filter(row_number()==1)
  mr_first_coal_time = max(first_coal_sims$time)
  
  last_coal_sims = pboot_ne_ind %>% group_by(n_sim) %>% filter(row_number()==n())
  oldest_last_coal_time = min(last_coal_sims$time)
  
  # Get common time axis for all simulations
  common_time_ax = approx(pboot_ne_ind$time, pboot_ne_ind$est_ne, xout=seq(mr_first_coal_time, oldest_last_coal_time, length.out=nrow(pboot_ne_ind)/length(sim_trees)), rule=2)$x
  print("Common time axis:")
  print(common_time_ax)
  
  # Split df by simulation index again and compute approximate Ne estimates based on common time axis
  pboot_ne_ind_common = split(pboot_ne_ind, pboot_ne_ind$n_sim)
  time_adj = true_ne_adj = ne_est_adj = nelb_adj = neub_adj = pboot_ne_adj = list()
  for(i in 1:length(pboot_ne_ind_common)) { 
    time_adj[[i]] = common_time_ax; #time_adj[[i]] = as.data.frame(time_adj[[i]])
    true_ne_adj[[i]] = alphaFun(time_adj[[i]]); #true_ne_adj[[i]] = as.data.frame(true_ne_adj[[i]])
    ne_est_adj[[i]] = approx(pboot_ne_ind_common[[i]]$time, pboot_ne_ind_common[[i]]$est_ne, xout=common_time_ax, rule=2)$y; ne_est_adj[[i]] = as.data.frame(ne_est_adj[[i]])
    nelb_adj[[i]] = approx(pboot_ne_ind_common[[i]]$time, pboot_ne_ind_common[[i]]$nelb, xout=common_time_ax, rule=2)$y; nelb_adj[[i]] = as.data.frame(nelb_adj[[i]])
    neub_adj[[i]] = approx(pboot_ne_ind_common[[i]]$time, pboot_ne_ind_common[[i]]$neub, xout=common_time_ax, rule=2)$y; neub_adj[[i]] = as.data.frame(neub_adj[[i]])
    pboot_ne_adj[[i]] = cbind(time_adj[[i]], true_ne_adj[[i]], nelb_adj[[i]], ne_est_adj[[i]], neub_adj[[i]])
    colnames(pboot_ne_adj[[i]]) = c("time","true_ne", "nelb", "est_ne", "neub")
  }
  # Join again now with common time axis
  pboot_ne_adj_df = bind_rows(pboot_ne_adj, .id = 'n_sim')
  print("Estimates using common time axis:")
  print(pboot_ne_adj_df)
  
  # Calculate Mean Absolute Error (MEA) and Root Mean Square Error (RMSE)
  mean_abs_error_df = pboot_ne_adj_df %>% group_by(n_sim) %>% summarise(mean_abs_error = sum(abs(true_ne - est_ne))/length(common_time_ax))
  lvls_mae <- str_sort(unique(mean_abs_error_df$n_sim), numeric = TRUE)
  mean_abs_error_df$n_sim <- factor(mean_abs_error_df$n_sim, levels = lvls_mae)
  #mean_abs_error_df$n_sim = factor(mean_abs_error_df$n_sim, levels=unique(mean_abs_error_df$n_sim))
  rmse_df = pboot_ne_adj_df %>% group_by(n_sim) %>% summarise(rmse_val = sqrt(mean((true_ne - est_ne)^2)))
  lvls_rmse <- str_sort(unique(rmse_df$n_sim), numeric = TRUE)
  rmse_df$n_sim <- factor(rmse_df$n_sim, levels = lvls_rmse)
  #rmse_df$n_sim = factor(rmse_df$n_sim, levels=unique(rmse_df$n_sim))
  # Got same results from Metrics package functions
  #mae_metrics = pboot_ne_adj_df %>% group_by(n_sim) %>% summarise(mae(actual = true_ne, predicted = est_ne))
  #rmse_metrics = pboot_ne_adj_df %>% group_by(n_sim) %>% summarise(rmse(actual = true_ne, predicted = est_ne))
  # Summary error and RMSE over all simulations
  mean_abs_error_avg = mean(mean_abs_error_df$mean_abs_error)
  rmse_avg = mean(rmse_df$rmse_val)
  
  # Calculate coverage probability (whether the true Ne value is covered [=1] or not [=0] by the estimated confidence intervals)
  cov_prob = pboot_ne_adj_df %>% mutate(cover = ifelse((true_ne >= nelb & true_ne <= neub), 1, 0))
  print("Coverage probability full df:")
  print(cov_prob)
  
  cov_prob_min = cov_prob %>% select(n_sim, time, cover)
  cov_prob_min = dcast(cov_prob_min,time~n_sim, value.var = "cover")
  print(cov_prob_min)
  cov_prob_matrix = as.matrix(cov_prob_min)
  rownames(cov_prob_matrix) = cov_prob_matrix[,1]
  cov_prob_matrix = cov_prob_matrix[,-1]
  print("Coverage probability matrix:")
  print(cov_prob_matrix)
  
  # Get coverage probability over entire time axis (Overall)
  cov_prob_avg = mean(cov_prob_matrix)
  print(paste("Overall coverage probability:",cov_prob_avg))
  # Get coverage probability for each common time point
  cov_prob_over_time = rowMeans(cov_prob_matrix)
  
  # Convert to df to plot using ggplot
  cov_prob_over_time_df = data.frame(cov_prob_over_time)
  cov_prob_over_time_df$time = rownames(cov_prob_over_time_df)
  cov_prob_over_time_df$time = as.numeric(cov_prob_over_time_df$time)
  rownames(cov_prob_over_time_df) = NULL
  cov_prob_over_time_df$cov_prob_over_time = as.numeric(cov_prob_over_time_df$cov_prob_over_time)
  # print(paste("Coverage probability for each time point in common time axis:",cov_prob_over_time_df))
  
  cp_out_pref = "coverage_plots/"
  error_out_pref = "error_plots/"
  system(paste0("mkdir -p ",cp_out_pref,dirname(out_path_cov_prob)))
  system(paste0("mkdir -p ",error_out_pref,dirname(out_path_error)))
  system(paste0("mkdir -p ",error_out_pref,dirname(out_path_rmse)))

  p1 = ggplot(data=cov_prob_over_time_df, aes(x=time, y=cov_prob_over_time)) + 
    geom_line(color="steelblue") + geom_point(color="steelblue") + coord_cartesian(ylim = c(0,1)) + 
    scale_x_continuous(labels=scaleFUN, breaks = common_time_ax) + #scale_y_continuous(expand = c(0, 0)) +
    labs(x = "Time since most recent sample", y = paste("Coverage probability for",length(sim_trees),"simulations"), caption = paste0("Overall coverage probability = ", scaleFUN(cov_prob_avg))) +
    theme_bw() + theme(plot.title = element_text(hjust=0.5), plot.caption = element_text(hjust=1),
                       axis.text.x = element_text(size=7, angle = 90, vjust = 0.5, hjust=1))

  print(paste0("Coverage probability for each time point plotted to: ",cp_out_pref,out_path_cov_prob))
  ggsave(paste0(cp_out_pref,out_path_cov_prob), plot=p1, units="in", width=7, height=7, dpi=600)
  
  p2 = ggplot(data=mean_abs_error_df, aes(x=as.numeric(as.character(n_sim)), y=mean_abs_error)) + 
    geom_line(color="steelblue") + geom_point(color="steelblue") +
    labs(x = "Simulation index", y = "Mean Absolute Error (MAE)", caption = paste0("Overall MAE = ", scaleFUN(mean_abs_error_avg))) +
    theme_bw() + theme(plot.title = element_text(hjust=0.5), plot.caption = element_text(hjust=1))
  
  print(paste0("Error (MAE) for each simulation index plotted to: ",error_out_pref,out_path_error))
  ggsave(paste0(error_out_pref,out_path_error), plot=p2, units="in", width=7, height=7, dpi=600)
  
  p3 = ggplot(data=rmse_df, aes(x=as.numeric(as.character(n_sim)), y=rmse_val)) + 
    geom_line(color="steelblue") + geom_point(color="steelblue") +
    labs(x = "Simulation index", y = "Root Mean Square Error (RMSE)", caption = paste0("Overall RMSE = ", scaleFUN(rmse_avg))) +
    theme_bw() + theme(plot.title = element_text(hjust=0.5), plot.caption = element_text(hjust=1))
  
  print(paste0("Root Mean Square Error (RMSE) for each simulation index plotted to: ",error_out_pref,out_path_rmse))
  ggsave(paste0(error_out_pref,out_path_rmse), plot=p3, units="in", width=7, height=7, dpi=600)
  
  return(list(pboot_ne_ind, pboot_ne_adj_df, cov_prob, cov_prob_matrix, cov_prob_avg, cov_prob_over_time, cov_prob_over_time_df, common_time_ax, mean_abs_error_df, rmse_df))
}

# Compare coverage probabilities for different models (skykappa, skygrid, and skygrowth) and the same Ne function
compare_cov_prob_models <- function(m1_cp, m2_cp, m3_cp, common_t_ax, out_path) {
  cp_out_pref = "coverage_plots/"
  system(paste0("mkdir -p ",cp_out_pref,dirname(out_path)))
  
  m1_cp$model = "Skykappa"; m2_cp$model = "Skygrid"; m3_cp$model = "Skygrowth"
  m_comb = rbind(m1_cp, m2_cp, m3_cp)
  m_comb$model = factor(m_comb$model, levels=c("Skykappa", "Skygrid", "Skygrowth"))
  
  p = ggplot(m_comb, aes(time, cov_prob_over_time, color = model)) + #shape=model
    geom_line() + geom_point() + scale_color_manual(name = "Model", values = c("steelblue", "darkred", "darkolivegreen")) +
    coord_cartesian(ylim = c(0,1)) +
    scale_x_continuous(labels=scaleFUN, breaks = common_t_ax) +
    labs(x = "Time since most recent sample", y = "Coverage probability") + #for 500 simulations
    theme_bw() + theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(size=7, angle = 90, vjust = 0.5, hjust=1))
  
  print(paste0("Coverage probability for each time point plotted to: ",cp_out_pref,out_path))
  ggsave(paste0(cp_out_pref,out_path), plot=p, units="in", width=7, height=7, dpi=600)
  
  return(p)
}

# Compare coverage probability estimates for the same model using different sample sizes to see precision loss when dropping tips
compare_same_model_diff_samp_size <- function(cp_more_tips, cp_intermed_tips, common_t_ax, out_path) {
  cp_out_pref = "coverage_plots/"
  system(paste0("mkdir -p ",cp_out_pref,dirname(out_path)))
  
  cp_more_tips$sample_size = 200; cp_intermed_tips$sample_size = 100
  cp_comb = rbind(cp_more_tips, cp_intermed_tips)
  cp_comb$sample_size = factor(cp_comb$sample_size, levels=c(200, 100))
  
  p = ggplot(cp_comb, aes(time, cov_prob_over_time, color = sample_size)) + 
    geom_line() + geom_point() + scale_color_manual(name = "Sample size", values = c("steelblue", "darkred", "darkolivegreen")) +
    coord_cartesian(ylim = c(0,1)) +
    scale_x_continuous(labels=scaleFUN, breaks = common_t_ax) +
    labs(x = "Time since most recent sample", y = "Coverage probability") + #for 500 simulations
    theme_bw() + theme(plot.title = element_text(hjust=0.5),axis.text.x = element_text(size=7, angle = 90, vjust = 0.5, hjust=1))
  
  print(paste0("Coverage probability for each time point plotted to: ",cp_out_pref,out_path))
  ggsave(paste0(cp_out_pref,out_path), plot=p, units="in", width=7, height=7, dpi=600)
  
  return(p)
}

config_rem_ax1 <- theme(axis.ticks.x = element_blank(),axis.text.x = element_blank(),axis.title.x = element_blank(),plot.margin = margin(r=0.5,t=0.5,unit="cm"))
config_rem_ax2 <- theme(plot.margin = margin(r=0.5,t=0.5,unit="cm"))

# Compare errors (MAE and RMSE) for different models (skykappa, skygrid, and skygrowth) and the same Ne function
compare_err_models <- function(m1_mae, m2_mae, m3_mae, m1_rmse, m2_rmse, m3_rmse, common_t_ax, out_path_mae, out_path_rmse) {
  err_out_pref = "error_plots/"
  system(paste0("mkdir -p ",err_out_pref,dirname(out_path_mae)))
  system(paste0("mkdir -p ",err_out_pref,dirname(out_path_rmse)))
  
  m1_mae$model = "Skykappa"; m2_mae$model = "Skygrid"; m3_mae$model = "Skygrowth"
  mae_comb = rbind(m1_mae, m2_mae, m3_mae)
  mae_comb$model = factor(mae_comb$model, levels=c("Skykappa", "Skygrid", "Skygrowth"))
  
  p1 = ggplot(mae_comb, aes(x=model, y=mean_abs_error, color=model)) +
  	geom_violin() + geom_sina(size=0.4) +
  	scale_color_manual(name = "Model", values = c("steelblue", "darkred", "darkolivegreen")) +
  	labs(x="Model", y="MAE across simulations") + coord_cartesian(ylim=c(0,20)) +
  	theme_bw() + theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7), legend.position="none")
  
  print(paste0("Mean absolute error (MAE) comparison for different models plotted to: ",err_out_pref,out_path_mae))
  ggsave(paste0(err_out_pref,out_path_mae), plot=p1, units="in", width=3.5, height=2, dpi=600)
  
  m1_rmse$model = "Skykappa"; m2_rmse$model = "Skygrid"; m3_rmse$model = "Skygrowth"
  rmse_comb = rbind(m1_rmse, m2_rmse, m3_rmse)
  rmse_comb$model = factor(rmse_comb$model, levels=c("Skykappa", "Skygrid", "Skygrowth"))
  
  p2 = ggplot(rmse_comb, aes(x=model, y=rmse_val, color=model)) + 
  	geom_violin() + geom_sina(size=0.4) +
  	scale_color_manual(name = "Model", values = c("steelblue", "darkred", "darkolivegreen")) +
  	labs(x="Model", y="RMSE across simulations") + coord_cartesian(ylim=c(0,20)) +
  	theme_bw() + theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7), legend.position="none")
  
  print(paste0("Root Mean Square Error (RMSE) comparison for different models plotted to: ",err_out_pref,out_path_rmse))
  ggsave(paste0(err_out_pref,out_path_rmse), plot=p2, units="in", width=3.5, height=2, dpi=600)
  
  return(list(p1, p2, mae_comb, rmse_comb))
}

# Combine all 9 estimates horizontally into a sina/violin plot (x-axis = model)
compare_err_models_combine_all <- function(mae_comb_m1, mae_comb_m2, rmse_comb_m1, rmse_comb_m2, out_path_mae, out_path_rmse) {
	
	err_out_pref = "error_plots/"
	system(paste0("mkdir -p ",err_out_pref,dirname(out_path_mae)))
	system(paste0("mkdir -p ",err_out_pref,dirname(out_path_rmse)))
	
	mae_comb_m1$sample_size <- 200; mae_comb_m2$sample_size <- 100
	rmse_comb_m1$sample_size <- 200; rmse_comb_m2$sample_size <- 100
	mae_comb_all <- rbind(mae_comb_m1, mae_comb_m2)
	rmse_comb_all <- rbind(rmse_comb_m1, rmse_comb_m2)
	
	mae_comb_all$sample_size <- factor(mae_comb_all$sample_size, levels=c(200, 100))
	rmse_comb_all$sample_size <- factor(rmse_comb_all$sample_size, levels=c(200, 100))
	
	p1 = ggplot(mae_comb_all, aes(x=model, y=mean_abs_error, color=sample_size)) + 
		geom_violin() + geom_sina(size=0.05) +
		scale_color_manual(name = "Sample size", values = c("steelblue", "darkred", "darkolivegreen")) +
		labs(x="Model", y="MAE across simulations") + coord_cartesian(ylim=c(0,30)) +
		theme_bw() + config_combined + theme(legend.position=c(0.08,0.80)) #=c(0.16,0.80)
	
	p1_f = p1 + facet_grid(. ~ model, scales="free_x")
	
	print(paste0("Mean Absolute Error (MAE) comparison for same model and different sample sizes (combined) plotted to: ",err_out_pref,out_path_mae))
	ggsave(paste0(err_out_pref,out_path_mae), plot=p1_f, units="in", width=6, height=3, dpi=600)
	
	p2 = ggplot(rmse_comb_all, aes(x=model, y=rmse_val, color=sample_size)) + 
		geom_violin() + geom_sina(size=0.05) +
		scale_color_manual(name = "Sample size", values = c("steelblue", "darkred", "darkolivegreen")) +
		labs(x="Model", y="RMSE across simulations") + coord_cartesian(ylim=c(0,30)) +
		theme_bw() + config_combined + theme(legend.position=c(0.08,0.80))
	
	p2_f = p2 + facet_grid(. ~ model, scales="free_x")
	
	print(paste0("Root mean squared error (RMSE) comparison for same model and different sample sizes (combined) plotted to: ",err_out_pref,out_path_rmse))
	ggsave(paste0(err_out_pref,out_path_rmse), plot=p2_f, units="in", width=6, height=3, dpi=600)
	
	return(list(mae_comb_all, mae_comb_all))
}

compare_err_same_model_diff_samp_size <- function(mae_more_tips, mae_intermed_tips, rmse_more_tips, rmse_intermed_tips, common_t_ax, out_path_mae, out_path_rmse) {
	err_out_pref = "error_plots/"
	system(paste0("mkdir -p ",err_out_pref,dirname(out_path_mae)))
	system(paste0("mkdir -p ",err_out_pref,dirname(out_path_rmse)))
	
	mae_more_tips$sample_size = 200; mae_intermed_tips$sample_size = 100
	mae_comb = rbind(mae_more_tips, mae_intermed_tips)
	mae_comb$sample_size = factor(mae_comb$sample_size, levels=c(200, 100))
	
	fct_lbls <- as_labeller(c(`200`="n=200",`100`="n=100"))
	
	p1 = ggplot(mae_comb, aes(x=sample_size, y=mean_abs_error, color=sample_size)) +
		geom_violin() + geom_sina(size=0.4) + geom_point(size=0.7) +
		scale_color_manual(name = "Sample size", values = c("steelblue", "darkred", "darkolivegreen")) +
		labs(x="Sample size", y="MAE across simulations") + coord_cartesian(ylim=c(0,20)) +
		theme_bw() + theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7), legend.position="none")
	
	print(paste0("Mean absolute error (MAE) comparison for same model and different sample sizes plotted to: ",err_out_pref,out_path_mae))
	ggsave(paste0(err_out_pref,out_path_mae), plot=p1, units="in", width=3.5, height=2, dpi=600)
	
	rmse_more_tips$sample_size = 200; rmse_intermed_tips$sample_size = 100
	rmse_comb = rbind(rmse_more_tips, rmse_intermed_tips)
	rmse_comb$sample_size = factor(rmse_comb$sample_size, levels=c(200, 100))
	
	p2 = ggplot(rmse_comb, aes(x=sample_size, y=rmse_val, color=sample_size)) + 
		geom_violin() + geom_sina(size=0.4) + geom_point(size=0.7) +
		scale_color_manual(name = "Sample size", values = c("steelblue", "darkred", "darkolivegreen")) +
		labs(x="Sample size", y="RMSE across simulations") + coord_cartesian(ylim=c(0,20)) +
		theme_bw() + theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7), legend.position="none")
	
	print(paste0("Root Mean Square Error (RMSE) comparison for same model and different sample sizes plotted to: ",err_out_pref,out_path_rmse))
	ggsave(paste0(err_out_pref,out_path_rmse), plot=p2, units="in", width=3.5, height=2, dpi=600)
	
	return(list(p1, p2, mae_comb, rmse_comb))
}

config_combined <- theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7), strip.text.x = element_blank(),
																									legend.key.size = unit(0.35, 'cm'),legend.key.height = unit(0.35, 'cm'),legend.key.width = unit(0.35, 'cm'),legend.title = element_text(size=6),legend.text = element_text(size=4))

# Combine all 9 estimates horizontally into a sina/violin plot (x-axis = sample size)
compare_err_same_model_diff_samp_size_combine_all <- function(mae_comb_m1, mae_comb_m2, mae_comb_m3, rmse_comb_m1, rmse_comb_m2, rmse_comb_m3, out_path_mae, out_path_rmse) {
	
	err_out_pref = "error_plots/"
	system(paste0("mkdir -p ",err_out_pref,dirname(out_path_mae)))
	system(paste0("mkdir -p ",err_out_pref,dirname(out_path_rmse)))
	
	mae_comb_m1$model <- "Skykappa"; mae_comb_m2$model <- "Skygrid"; mae_comb_m3$model <- "Skygrowth"
	rmse_comb_m1$model <- "Skykappa"; rmse_comb_m2$model <- "Skygrid"; rmse_comb_m3$model <- "Skygrowth"
	mae_comb_all <- rbind(mae_comb_m1, mae_comb_m2, mae_comb_m3)
	rmse_comb_all <- rbind(rmse_comb_m1, rmse_comb_m2, rmse_comb_m3)

	mae_comb_all$model = factor(mae_comb_all$model, levels=c("Skykappa", "Skygrid", "Skygrowth"))
	rmse_comb_all$model = factor(rmse_comb_all$model, levels=c("Skykappa", "Skygrid", "Skygrowth"))
	
	p1 = ggplot(mae_comb_all, aes(x=sample_size, y=mean_abs_error, color=model)) + 
		geom_violin() + geom_sina(size=0.05) +
		scale_color_manual(name = "Model", values = c("steelblue", "darkred", "darkolivegreen")) +
		labs(x="Sample size", y="MAE across simulations") + coord_cartesian(ylim=c(0,30)) +
		theme_bw() + config_combined + theme(legend.position=c(0.08,0.80)) #=c(0.16,0.80)
	
	p1_f = p1 + facet_grid(. ~ sample_size, scales="free_x")
	
	print(paste0("Mean Absolute Error (MAE) comparison for same model and different sample sizes (combined) plotted to: ",err_out_pref,out_path_mae))
	ggsave(paste0(err_out_pref,out_path_mae), plot=p1_f, units="in", width=6, height=3, dpi=600)
	
	p2 = ggplot(rmse_comb_all, aes(x=sample_size, y=rmse_val, color=model)) + 
		geom_violin() + geom_sina(size=0.05) +
		scale_color_manual(name = "Model", values = c("steelblue", "darkred", "darkolivegreen")) +
		labs(x="Sample size", y="RMSE across simulations") + coord_cartesian(ylim=c(0,30)) +
		theme_bw() + config_combined + theme(legend.position=c(0.08,0.80))
	
	p2_f = p2 + facet_grid(. ~ sample_size, scales="free_x")
	
	print(paste0("Root mean squared error (RMSE) comparison for same model and different sample sizes (combined) plotted to: ",err_out_pref,out_path_rmse))
	ggsave(paste0(err_out_pref,out_path_rmse), plot=p2_f, units="in", width=6, height=3, dpi=600)
	
	return(list(mae_comb_all, mae_comb_all))
}

# Since estimates for different sample sizes have different time axis, get ranges and approx values to get same time axis for all to fill table
# This way we can have only 1 table per function instead of 3
generate_summary_table_cov_prob <- function(alphaFun, t200_m1, t200_m2, t200_m3, t100_m1, t100_m2, t100_m3, out_path_cp) {
  # 1 = n200+skykappa; 2 = n200+skygrid; 3 = n200+skygrowth; 4 = n100+skykappa; 5 = n100+skygrid; 6 = n100+skygrowth
  t200_m1$id = 1; t200_m2$id = 2; t200_m3$id = 3;
  t100_m1$id = 4; t100_m2$id = 5; t100_m3$id = 6;
  
  all_estims = rbind(t200_m1, t200_m2, t200_m3, t100_m1, t100_m2, t100_m3)
  
  # Get boundaries for common time axis across different sample sizes
  first_coal_sims = all_estims %>% group_by(id) %>% filter(row_number()==1)
  mr_first_coal_time = max(first_coal_sims$time)
  last_coal_sims = all_estims %>% group_by(id) %>% filter(row_number()==n())
  oldest_last_coal_time = min(last_coal_sims$time)
  
  common_time_ax = approx(all_estims$time, all_estims$est_ne, xout=seq(mr_first_coal_time, oldest_last_coal_time, length.out=nrow(all_estims)/(nsim*9)), rule=2)$x
  print("Common time axis:")
  print(common_time_ax)
  
  all_estims_nsim_id = split(all_estims, list(all_estims$n_sim, all_estims$id), drop = TRUE)
  pboot_ne_adj = list()
  for(i in seq_along(all_estims_nsim_id)) {
    time_adj = common_time_ax; #time_adj = as.data.frame(time_adj)
    true_ne_adj = alphaFun(time_adj); #true_ne_adj = as.data.frame(time_adj)
    ne_est_adj = approx(all_estims_nsim_id[[i]]$time, all_estims_nsim_id[[i]]$est_ne, xout=common_time_ax, rule=2)$y; ne_est_adj = as.data.frame(ne_est_adj)
    nelb_adj = approx(all_estims_nsim_id[[i]]$time, all_estims_nsim_id[[i]]$nelb, xout=common_time_ax, rule=2)$y; nelb_adj = as.data.frame(nelb_adj)
    neub_adj = approx(all_estims_nsim_id[[i]]$time, all_estims_nsim_id[[i]]$neub, xout=common_time_ax, rule=2)$y; neub_adj = as.data.frame(neub_adj)
    nsim = all_estims_nsim_id[[i]]$n_sim[1:length(common_time_ax)]
    id = all_estims_nsim_id[[i]]$id[1:length(common_time_ax)]
    pboot_ne_adj[[i]] = cbind(time_adj, true_ne_adj, nelb_adj, ne_est_adj, neub_adj, nsim, id)
    colnames(pboot_ne_adj[[i]]) = c("time","true_ne", "nelb", "est_ne", "neub", "n_sim", "id")
  }
  # Join again now with common time axis
  pboot_ne_nsim_id_adj_df = bind_rows(pboot_ne_adj)
  # Calculate coverage probability (whether the true Ne value is covered [=1] or not [=0] by the estimated confidence intervals)
  cov_prob = pboot_ne_nsim_id_adj_df %>% mutate(cover_table = ifelse((true_ne >= nelb & true_ne <= neub), 1, 0))
  # View(cov_prob)
  cov_prob_min = cov_prob %>% select(time, cover_table, n_sim, id)
  cov_prob_min_splt = split(cov_prob_min, cov_prob_min$id)
  print(length(cov_prob_min_splt))
  cov_prob_min_id = cov_prob_avg = cov_prob_over_time = cov_prob_over_time_df = list()
  for(i in 1:6) {
    cov_prob_min_id[[i]] = dcast(cov_prob_min_splt[[i]],time~n_sim, value.var = "cover_table")
    cov_prob_matrix = as.matrix(cov_prob_min_id[[i]])
    rownames(cov_prob_matrix) = cov_prob_matrix[,1]
    cov_prob_matrix = cov_prob_matrix[,-1]
    # Get coverage probability over entire time axis (Overall)
    cov_prob_avg[[i]] = mean(cov_prob_matrix)
    # Get coverage probability for each common time point
    cov_prob_over_time[[i]] = rowMeans(cov_prob_matrix)
    cov_prob_over_time_df[[i]] = data.frame(cov_prob_over_time[[i]])
    cov_prob_over_time_df[[i]]$time = rownames(cov_prob_over_time_df[[i]])
    cov_prob_over_time_df[[i]]$time = as.numeric(cov_prob_over_time_df[[i]]$time)
    rownames(cov_prob_over_time_df[[i]]) = NULL
    cov_prob_over_time_df[[i]]$cov_prob_over_time = as.numeric(cov_prob_over_time_df[[i]]$cov_prob_over_time)
    cov_prob_over_time_df[[i]] = cov_prob_over_time_df[[i]] %>% select(time, cov_prob_over_time)
  }
  axis_len = lapply(cov_prob_over_time_df, nrow)
  min_axis_len = min(unlist(axis_len))
  for(i in 1:length(cov_prob_over_time_df)) { 
    cov_prob_over_time_df[[i]] = cov_prob_over_time_df[[i]][1:min_axis_len,]
  }
  cov_prob_all = Reduce(cbind, cov_prob_over_time_df)
  to_del <- seq(3, ncol(cov_prob_all), 2)
  cov_prob_all_f <-  cov_prob_all[,-to_del]
  colnames(cov_prob_all_f) = c("time","n200+skykappa","n200+skygrid","n200+skygrowth","n100+skykappa","n100+skygrid","n100+skygrowth")
  
  summ_cp_pref = "tables/"
  system(paste0("mkdir -p ",summ_cp_pref))
  write.table(format(cov_prob_all_f,digits=4), file=paste0(summ_cp_pref,out_path_cp), quote=FALSE, row.names=TRUE, col.names=NA, sep=",")
  return(list(cov_prob_all_f, cov_prob_avg))
}

generate_summary_table_rmse <- function(alphaFun, t200_m1, t200_m2, t200_m3, t100_m1, t100_m2, t100_m3, out_path_rmse) {
  t200_m1$id = 1; t200_m2$id = 2; t200_m3$id = 3;
  t100_m1$id = 4; t100_m2$id = 5; t100_m3$id = 6;
  all_estims = rbind(t200_m1, t200_m2, t200_m3, t100_m1, t100_m2, t100_m3)
  rmse_stat_df = all_estims %>% group_by(id) %>% summarise(avg_rmse = mean(rmse_val), median_rmse = median(rmse_val), iqr_rmse = IQR(rmse_val))
  rmse_stat_df = rmse_stat_df %>% select(avg_rmse, median_rmse, iqr_rmse)
  rmse_stat_df = t(rmse_stat_df)
  colnames(rmse_stat_df) = c("n200+skykappa","n200+skygrid","n200+skygrowth","n100+skykappa","n100+skygrid","n100+skygrowth")
  summ_err_pref = "tables/"
  system(paste0("mkdir -p ",summ_err_pref))
  write.table(format(rmse_stat_df,digits=4), file=paste0(summ_err_pref,out_path_rmse), quote=FALSE, row.names=TRUE, col.names=NA, sep=",")
  return(rmse_stat_df)
}

plot_ne_trajectories <- function(sim_out, sample_size, model, out_path) {
	sim_out_indiv_sims <- split(sim_out[[1]], as.integer(sim_out[[1]]$n_sim))
	sim_out_p1_ne <- list()
	#sim_out_indiv_sims_time <- sim_out_indiv_sims_ne_ci <- 
	sim_out_indiv_sims_time <- sim_out_indiv_sims_mlsky_obj <- list()
	for(l in 1:length(sim_out_indiv_sims)) {
		sim_out_indiv_sims[[l]] <- sim_out_indiv_sims[[l]] %>% select(time, nelb, est_ne, neub)
		sim_out_indiv_sims_time[[l]] <- sim_out_indiv_sims[[l]] %>% select(time)
		sim_out_indiv_sims_time[[l]]$time <- sim_out_indiv_sims_time[[l]][["time"]]
		sim_out_indiv_sims_mlsky_obj[[l]] <- sim_out_indiv_sims_time[[l]]
		sim_out_indiv_sims[[l]] <- sim_out_indiv_sims[[l]] %>% select(nelb, est_ne, neub)
		colnames(sim_out_indiv_sims[[l]]) <- c("nelb","ne","neub")
		#sim_out_indiv_sims_mlsky_obj[[l]] <- sim_out_indiv_sims[[l]]
		sim_out_indiv_sims_mlsky_obj[[l]]$ne_ci <- sim_out_indiv_sims[[l]]
		class(sim_out_indiv_sims_mlsky_obj[[l]]) <- "mlskygrid"
		# sim_out_indiv_sims_mlsky_obj[[l]] <- c(sim_out_indiv_sims_time[[l]], sim_out_indiv_sims_ne_ci[[l]]$ne_ci)
		#sim_out_indiv_sims_mlsky_obj[[l]] <- 
		sim_out_p1_ne[[l]] <- plot(sim_out_indiv_sims_mlsky_obj[[l]], ggplot=TRUE, logy=FALSE) + ggplot2::ylab("Ne") + ggplot2::xlab("Time") + ggplot2::coord_cartesian(ylim=c(0,100)) + ggplot2::ggtitle(glue("{sample_size},{model}")) + theme_classic()
		system(glue("mkdir -p {out_path}/{model}/{sample_size}"))
		ggsave(glue("{out_path}/{model}/{sample_size}/sim_id{l}.png"), plot=sim_out_p1_ne[[l]], units="in", width=5, height=3.5, dpi=300) #6,4,600
	}
}

plot_ne_trajectories <- function(sim_out, sample_size, model, out_path, max_v) {
	sim_out_indiv_sims <- split(sim_out[[1]], as.integer(sim_out[[1]]$n_sim))
	sim_out_p1_ne <- list()
	#sim_out_indiv_sims_time <- sim_out_indiv_sims_ne_ci <- 
	sim_out_indiv_sims_time <- sim_out_indiv_sims_mlsky_obj <- list()
	for(l in 1:length(sim_out_indiv_sims)) {
		sim_out_indiv_sims[[l]] <- sim_out_indiv_sims[[l]] %>% select(time, nelb, est_ne, neub)
		sim_out_indiv_sims_time[[l]] <- sim_out_indiv_sims[[l]] %>% select(time)
		sim_out_indiv_sims_time[[l]]$time <- sim_out_indiv_sims_time[[l]][["time"]]
		sim_out_indiv_sims_mlsky_obj[[l]] <- sim_out_indiv_sims_time[[l]]
		sim_out_indiv_sims[[l]] <- sim_out_indiv_sims[[l]] %>% select(nelb, est_ne, neub)
		colnames(sim_out_indiv_sims[[l]]) <- c("nelb","ne","neub")
		#sim_out_indiv_sims_mlsky_obj[[l]] <- sim_out_indiv_sims[[l]]
		sim_out_indiv_sims_mlsky_obj[[l]]$ne_ci <- sim_out_indiv_sims[[l]]
		class(sim_out_indiv_sims_mlsky_obj[[l]]) <- "mlskygrid"
		# sim_out_indiv_sims_mlsky_obj[[l]] <- c(sim_out_indiv_sims_time[[l]], sim_out_indiv_sims_ne_ci[[l]]$ne_ci)
		#sim_out_indiv_sims_mlsky_obj[[l]] <- 
		sim_out_p1_ne[[l]] <- plot(sim_out_indiv_sims_mlsky_obj[[l]], ggplot=TRUE, logy=FALSE) + ggplot2::ylab("Ne") + ggplot2::xlab("Time") + ggplot2::coord_cartesian(ylim=c(0,max_v)) + ggplot2::ggtitle(glue("{sample_size},{model}")) + theme_classic()
		system(glue("mkdir -p {out_path}/{model}/{sample_size}"))
		ggsave(glue("{out_path}/{model}/{sample_size}/sim_id{l}.png"), plot=sim_out_p1_ne[[l]], units="in", width=5, height=3.5, dpi=300) #6,4,600
	}
}

plot_ne_trajectories_true_ne_without_ci <- function(sim_out, sample_size, model, out_path) {
	sim_out_indiv_sims <- split(sim_out[[2]], as.integer(sim_out[[2]]$n_sim))
	#View(sim_out_indiv_sims[[1]])
	sim_out_p1_ne <- list()
	sim_out_indiv_sims_time <- sim_out_indiv_sims_mlsky_obj <- list()
	for(l in 1:length(sim_out_indiv_sims)) {
		sim_out_indiv_sims[[l]] <- sim_out_indiv_sims[[l]] %>% select(time, true_ne, est_ne)
		colnames(sim_out_indiv_sims[[l]]) <- c("time","true_ne","est_ne")
		sim_out_p1_ne[[l]] <- ggplot(sim_out_indiv_sims[[l]] , aes(x=time)) +
			geom_line(aes(y = true_ne), color = "darkred") + 
			geom_line(aes(y = est_ne), color = "steelblue")
		system(glue("mkdir -p {out_path}/{model}/{sample_size}"))
		ggsave(glue("{out_path}/{model}/{sample_size}/sim_id{l}.png"), plot=sim_out_p1_ne[[l]], units="in", width=5, height=3.5, dpi=300) #6,4,600
	}
}
