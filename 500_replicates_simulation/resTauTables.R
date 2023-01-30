# Clean all variables from environment
rm(list=ls())

library(glue)
library(data.table)

options(scipen=1, digits=8)

source("500_replicates_simulation/commonFunctions.R")

sinusFun=function(x){sin(x)*10+12}

# This simulations will fix `res`, `tau` or both
get_nsim_estimates2 <- function(sim_trees, sampDates, alphaFun, model, which_fixed, res_val, tau_val) {
	print("Fit")
	sts = sampDates; names(sts) = sim_trees[[1]]$tip.label
	if(which_fixed=="res") {
		fit = pbmclapply(X=sim_trees, FUN=mlskygrid, sampleTimes=sts, res=res_val, tau=NULL, tau_lower = 0.001, tau_upper = 20, model = model, ncpu=ncpu)
	}else if(which_fixed=="tau") {
		fit = pbmclapply(X=sim_trees, FUN=mlskygrid, sampleTimes=sts, res=NULL, tau=tau_val, model = model, ncpu=ncpu)
	} else if(which_fixed=="both") { # both fixed
		fit = pbmclapply(X=sim_trees, FUN=mlskygrid, sampleTimes=sts, res=res_val, tau=tau_val, model = model, ncpu=ncpu)
	} else {
		stop("Valid options for `which_fixed` are: `res`, `tau` and `both`")
	}
	pboot = pboot_ne = list()
	for(i in 1:length(sim_trees)){
		print(paste("===",i,"==="))
		if(length(fit[[i]]$time) <= 2) {
			next
		}
		pboot[[i]] = parboot(fit[[i]], nrep=200, ncpu=ncpu) #dd=T (use ddSimCoal or not?)
		pboot_ne[[i]] = pboot[[i]]["ne_ci"]
		pboot_ne[[i]]$nelb = pboot[[i]]$ne_ci[,1]
		pboot_ne[[i]]$ne = pboot[[i]]$ne_ci[,2]
		pboot_ne[[i]]$neub = pboot[[i]]$ne_ci[,3]
		pboot_ne[[i]]$time = pboot[[i]]$time; pboot_ne[[i]]$time = as.data.frame(pboot_ne[[i]]$time)
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
	
	# Calculate Root Mean Square Error (RMSE)
	# mean_abs_error_df = pboot_ne_adj_df %>% group_by(n_sim) %>% summarise(mean_abs_error = sum(abs(true_ne - est_ne))/length(common_time_ax))
	# lvls_mae <- str_sort(unique(mean_abs_error_df$n_sim), numeric = TRUE)
	# mean_abs_error_df$n_sim <- factor(mean_abs_error_df$n_sim, levels = lvls_mae)
	#mean_abs_error_df$n_sim = factor(mean_abs_error_df$n_sim, levels=unique(mean_abs_error_df$n_sim))
	rmse_df = pboot_ne_adj_df %>% group_by(n_sim) %>% summarise(rmse_val = sqrt(mean((true_ne - est_ne)^2)))
	lvls_rmse <- str_sort(unique(rmse_df$n_sim), numeric = TRUE)
	rmse_df$n_sim <- factor(rmse_df$n_sim, levels = lvls_rmse)
	#rmse_df$n_sim = factor(rmse_df$n_sim, levels=unique(rmse_df$n_sim))
	# Summary RMSE over all simulations
	#mean_abs_error_avg = mean(mean_abs_error_df$mean_abs_error)
	rmse_avg = mean(rmse_df$rmse_val)
	
	return(list(pboot_ne_adj_df, rmse_df, common_time_ax)) #pboot_ne_ind, mean_abs_error_df
}

# To use same simulated trees as in other analysis, load `t50_sinus_trees` (using n=50 samples from sinus function)
t50_sinus_trees <- readRDS(file="500_replicates_simulation/env/t50_sinus_trees.rds")

res_choices <- c(5,10,15,20,50) #5
res_choices_both <- c(5,10,20)
tau_choices <- c(0.001, 1, 20) #3

# `res` fixed, `tau estimated`
t50_sinus_m1_res_fix <- list()
for(i in 1:length(res_choices)) {
	print(glue("===res={res_choices[i]}==="))
	t50_sinus_m1_res_fix[[i]] = get_nsim_estimates2(t50_sinus_trees, sampleDates_50, sinusFun, 1, which_fixed="res", res_val=res_choices[i], tau_val=NULL)
}

#saveRDS(t50_sinus_m1_res_fix, file = "500_replicates_simulation/env/t50_sinus_m1_res_fix.rds")

# `tau` fixed, `res estimated`
t50_sinus_m1_tau_fix <- list()
for(j in 1:length(tau_choices)) {
	print(glue("===tau={tau_choices[j]}==="))
	t50_sinus_m1_tau_fix[[j]] = get_nsim_estimates2(t50_sinus_trees, sampleDates_50, sinusFun, 1, which_fixed="tau", res_val=NULL, tau_val=tau_choices[j])
}

#saveRDS(t50_sinus_m1_tau_fix, file = "500_replicates_simulation/env/t50_sinus_m1_tau_fix.rds")

t50_sinus_m1_both_fix <- matrix(data=list(), nrow=length(res_choices_both), ncol=length(tau_choices)) #5x3
for(k in 1:length(res_choices_both)) {
	for(l in 1:length(tau_choices)) {
		print(glue("===res={res_choices_both[k]}; tau={tau_choices[l]}==="))
		t50_sinus_m1_both_fix[[k,l]] = get_nsim_estimates2(t50_sinus_trees, sampleDates_50, sinusFun, 1, which_fixed="both", res_val=res_choices_both[k], tau_val=tau_choices[l])
	}
}

#saveRDS(t50_sinus_m1_both_fix, file = "500_replicates_simulation/env/t50_sinus_m1_both_fix.rds")

generate_summary_table_rmse2 <- function(alphaFun, obj_type, rmse_obj, coln, out_path_rmse) { #t200_m1, t200_m2, t200_m3, t50_m1, t50_m2, t50_m3, t20_m1, t20_m2, t20_m3
	res_list <- list()
	res_df <- data.frame()
	if(obj_type=="list") {
		for(i in 1:length(rmse_obj)) {
			res_list[[i]] <- rmse_obj[[i]][[2]]
		}
		res_df <- rbindlist(res_list, idcol="id")
	} else if(obj_type=="matrix") {
		cnt <- 0
		for(j in 1:length(res_choices_both)) {
			for(k in 1:length(tau_choices)) {
				#print(glue("{j},{k}"))
				cnt <- cnt+1
				#print(cnt)
				res_list[[cnt]] <- rmse_obj[[j,k]][[2]]
			}
		}
		res_df <- rbindlist(res_list, idcol="id")
	} else {
		stop("Valid options for `obj_type` are: `list` and `matrix`")
	}
	rmse_stat_df = res_df %>% group_by(id) %>% summarise(avg_rmse = mean(rmse_val), median_rmse = median(rmse_val), iqr_rmse = IQR(rmse_val))
	rmse_stat_df = rmse_stat_df %>% select(avg_rmse, median_rmse, iqr_rmse)
	rmse_stat_df = t(rmse_stat_df)
	colnames(rmse_stat_df) = coln
	summ_err_pref = "tables/"
	system(paste0("mkdir -p ",summ_err_pref))
	write.table(format(rmse_stat_df,digits=4), file=paste0(summ_err_pref,out_path_rmse), quote=FALSE, row.names=TRUE, col.names=NA, sep=",")
	return(rmse_stat_df)
}

generate_summary_table_rmse2(sinusFun, obj_type="list", rmse_obj=t50_sinus_m1_res_fix, coln=c("res=5","res=10","res=15","res=20","res=50"), out_path_rmse="res_fixed_rmse.csv")
generate_summary_table_rmse2(sinusFun, obj_type="list", rmse_obj=t50_sinus_m1_tau_fix, coln=c("tau=0.001","tau=1","tau=20"), out_path_rmse="tau_fixed_rmse.csv")
generate_summary_table_rmse2(sinusFun, obj_type="matrix", rmse_obj=t50_sinus_m1_both_fix, coln=c("res=5;tau=0.001","res=5;tau=1","res=5;tau=20","res=10;tau=0.001","res=10;tau=1","res=10;tau=20","res=20;tau=0.001","res=20;tau=1","res=20;tau=20"), out_path_rmse="both_fixed_rmse.csv")
