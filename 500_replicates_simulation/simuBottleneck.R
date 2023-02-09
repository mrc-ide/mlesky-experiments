# Clean all variables from environment
rm(list=ls())

# IMPORTANT: This will take a while and require significant CPU resources to run, so please load the RData file as in the line below to have the results directly
load(file="500_replicates_simulation/env/simuBottle.RData")
source("500_replicates_simulation/commonFunctions.R")
start <- Sys.time()

bottleFun=function(x){
	res = vector()
	for(i in 1:length(x)) {
		if(x[i] < -15 || x[i] > -10) {
			res[i]=10
		}
		else {
			res[i]=1
		}
	}
	return(res)
}

# IMPORTANT: Begin skipping the lines below if you have loaded the environment (RData file) previously
t200_bottle_trees = simulate_ntrees(sampleDates_200, bottleFun) #to use same trees in `resTauTablesBottle.R`: saveRDS(t200_bottle_trees, file = "500_replicates_simulation/env/t200_bottle_trees.rds")
t100_bottle_trees = simulate_ntrees(sampleDates_100, bottleFun)

#sampleDates_200
t200_bottle_m1 = get_nsim_estimates(t200_bottle_trees, bottleFun, 1, res_choice=NULL, "bottle/bottle_t200_m1.png", "bottle/bottle_t200_m1_mae.png", "bottle/bottle_t200_m1_rmse.png")
t200_bottle_m2 = get_nsim_estimates(t200_bottle_trees, bottleFun, 2, res_choice=NULL, "bottle/bottle_t200_m2.png", "bottle/bottle_t200_m2_mae.png", "bottle/bottle_t200_m2_rmse.png")
t200_bottle_m3 = get_nsim_estimates(t200_bottle_trees, bottleFun, 3, res_choice=NULL, "bottle/bottle_t200_m3.png", "bottle/bottle_t200_m3_mae.png", "bottle/bottle_t200_m3_rmse.png")

#sampleDates_100
t100_bottle_m1 = get_nsim_estimates(t100_bottle_trees, bottleFun, 1, res_choice=NULL, "bottle/bottle_t100_m1.png", "bottle/bottle_t100_m1_mae.png", "bottle/bottle_t100_m1_rmse.png")
t100_bottle_m2 = get_nsim_estimates(t100_bottle_trees, bottleFun, 2, res_choice=NULL, "bottle/bottle_t100_m2.png", "bottle/bottle_t100_m2_mae.png", "bottle/bottle_t100_m2_rmse.png")
t100_bottle_m3 = get_nsim_estimates(t100_bottle_trees, bottleFun, 3, res_choice=NULL, "bottle/bottle_t100_m3.png", "bottle/bottle_t100_m3_mae.png", "bottle/bottle_t100_m3_rmse.png")
# IMPORTANT: Finish skipping here. The lines below will reproduce the plots!

# compare different models for the same function and sample size
t200_bottle_models = compare_cov_prob_models(t200_bottle_m1[[7]], t200_bottle_m2[[7]], t200_bottle_m3[[7]], t200_bottle_m1[[8]], "comp_models/cov_prob/bottle_t200.png")
t100_bottle_models = compare_cov_prob_models(t100_bottle_m1[[7]], t100_bottle_m2[[7]], t100_bottle_m3[[7]], t100_bottle_m1[[8]], "comp_models/cov_prob/bottle_t100.png")

ggarrange(t200_bottle_models, t100_bottle_models, nrow=2, ncol=1, labels=c("n=200","n=100"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("coverage_plots/comp_models/bottle_models.pdf", units="in", width=9, height=8, dpi=600)

# compare different sample sizes for same function and model
tx_bottle_skykappa = compare_same_model_diff_samp_size(t200_bottle_m1[[7]], t100_bottle_m1[[7]], t200_bottle_m1[[8]], "comp_models/cov_prob/bottle_diff_samp_skykappa.png")
tx_bottle_skygrid = compare_same_model_diff_samp_size(t200_bottle_m2[[7]], t100_bottle_m2[[7]], t200_bottle_m2[[8]], "comp_models/cov_prob/bottle_diff_samp_skygrid.png")
tx_bottle_skygrowth = compare_same_model_diff_samp_size(t200_bottle_m3[[7]], t100_bottle_m3[[7]], t200_bottle_m3[[8]], "comp_models/cov_prob/bottle_diff_samp_skygrowth.png")

ggarrange(tx_bottle_skykappa, tx_bottle_skygrid, tx_bottle_skygrowth, nrow=3, ncol=1, labels=c("Skykappa","Skygrid","Skygrowth"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("coverage_plots/comp_models/bottle_diff_samp_sizes.pdf", units="in", width=9, height=12, dpi=600)

# compare error considering different models for the same function and sample size
t200_bottle_models_err = compare_err_models(t200_bottle_m1[[9]], t200_bottle_m2[[9]], t200_bottle_m3[[9]], t200_bottle_m1[[10]], t200_bottle_m2[[10]], t200_bottle_m3[[10]], t200_bottle_m1[[8]], "comp_models/mae/bottle_t200.png", "comp_models/rmse/bottle_t200.png")
t100_bottle_models_err = compare_err_models(t100_bottle_m1[[9]], t100_bottle_m2[[9]], t100_bottle_m3[[9]], t100_bottle_m1[[10]], t100_bottle_m2[[10]], t100_bottle_m3[[10]], t100_bottle_m1[[8]], "comp_models/mae/bottle_t100.png", "comp_models/rmse/bottle_t100.png")

# Sina/violin combined plot
tx_bottle_err_models_all <- compare_err_models_combine_all(t200_bottle_models_err[[3]], t100_bottle_models_err[[3]], t200_bottle_models_err[[4]], t100_bottle_models_err[[4]], "comp_models/mae/bottle_diff_models_COMBINED.pdf", "comp_models/rmse/bottle_diff_models_COMBINED.pdf")

# compare error considering different sample sizes for same function and model
tx_bottle_skykappa_err = compare_err_same_model_diff_samp_size(t200_bottle_m1[[9]], t100_bottle_m1[[9]], t200_bottle_m1[[10]], t100_bottle_m1[[10]], t200_bottle_m1[[8]], "comp_models/mae/bottle_diff_samp_skykappa.png", "comp_models/rmse/bottle_diff_samp_skykappa.png")
tx_bottle_skygrid_err = compare_err_same_model_diff_samp_size(t200_bottle_m2[[9]], t100_bottle_m2[[9]], t200_bottle_m2[[10]], t100_bottle_m2[[10]], t200_bottle_m2[[8]], "comp_models/mae/bottle_diff_samp_skygrid.png", "comp_models/rmse/bottle_diff_samp_skygrid.png")
tx_bottle_skygrowth_err = compare_err_same_model_diff_samp_size(t200_bottle_m3[[9]], t100_bottle_m3[[9]], t200_bottle_m3[[10]], t100_bottle_m3[[10]], t200_bottle_m3[[8]], "comp_models/mae/bottle_diff_samp_skygrowth.png", "comp_models/rmse/bottle_diff_samp_skygrowth.png")

# Sina/violin combined plot
tx_bottle_err_all <- compare_err_same_model_diff_samp_size_combine_all(tx_bottle_skykappa_err[[3]], tx_bottle_skygrid_err[[3]], tx_bottle_skygrowth_err[[3]], tx_bottle_skykappa_err[[4]], tx_bottle_skygrid_err[[4]], tx_bottle_skygrowth_err[[4]], "comp_models/mae/bottle_diff_samp_COMBINED.pdf", "comp_models/rmse/bottle_diff_samp_COMBINED.pdf")

end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
print(paste("Total time elapsed: ",total_time,"mins"))

summary_bottle_cp = generate_summary_table_cov_prob(bottleFun, t200_bottle_m1[[3]], t200_bottle_m2[[3]], t200_bottle_m3[[3]],
                                                    t100_bottle_m1[[3]],t100_bottle_m2[[3]], t100_bottle_m3[[3]], "bottle_cp.csv")

summary_bottle_rmse = generate_summary_table_rmse(bottleFun, t200_bottle_m1[[10]], t200_bottle_m2[[10]], t200_bottle_m3[[10]],
                                                    t100_bottle_m1[[10]],t100_bottle_m2[[10]], t100_bottle_m3[[10]], "bottle_rmse.csv")

# If running from scratch, the line below will have all the variables generated
system("mkdir -p 500_replicates_simulation/env/")
save.image(file="500_replicates_simulation/env/simuBottle.RData")

plot_ne_trajectories(t200_bottle_m1, "200", "skykappa", "trajectories_bottle", 20)
plot_ne_trajectories(t200_bottle_m2, "200", "skygrid", "trajectories_bottle", 20)
plot_ne_trajectories(t200_bottle_m3, "200", "skygrowth", "trajectories_bottle", 20)

plot_ne_trajectories_true_ne_without_ci(t200_bottle_m1, "200", "skykappa", "trajectories_bottle_true_ne")
plot_ne_trajectories_true_ne_without_ci(t200_bottle_m2, "200", "skygrid", "trajectories_bottle_true_ne")
plot_ne_trajectories_true_ne_without_ci(t200_bottle_m3, "200", "skygrowth", "trajectories_bottle_true_ne")