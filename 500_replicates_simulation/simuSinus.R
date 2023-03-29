# Clean all variables from environment
rm(list=ls())

# IMPORTANT: This will take a while and require significant CPU resources to run, so please load the RData file as in the line below to have the results directly
load(file="500_replicates_simulation/env/simuSinus.RData")
source("500_replicates_simulation/commonFunctions.R")
start <- Sys.time()

sinusFun=function(x){sin(x)*10+12}

# IMPORTANT: Begin skipping the lines below if you have loaded the environment (RData file) previously
t200_sinus_trees = simulate_ntrees(sampleDates_200, sinusFun)
t100_sinus_trees = simulate_ntrees(sampleDates_100, sinusFun) 

#sampleDates_200
t200_sinus_m1 = get_nsim_estimates(t200_sinus_trees, sinusFun, 1, res_choice=30, "sinus/sinus_t200_m1.png", "sinus/sinus_t200_m1_mae.png", "sinus/sinus_t200_m1_rmse.png")
t200_sinus_m2 = get_nsim_estimates(t200_sinus_trees, sinusFun, 2, res_choice=30, "sinus/sinus_t200_m2.png", "sinus/sinus_t200_m2_mae.png", "sinus/sinus_t200_m2_rmse.png")
t200_sinus_m3 = get_nsim_estimates(t200_sinus_trees, sinusFun, 3, res_choice=30, "sinus/sinus_t200_m3.png", "sinus/sinus_t200_m3_mae.png", "sinus/sinus_t200_m3_rmse.png")

#sampleDates_100
t100_sinus_m1 = get_nsim_estimates(t100_sinus_trees, sinusFun, 1, res_choice=30, "sinus/sinus_t100_m1.png", "sinus/sinus_t100_m1_mae.png", "sinus/sinus_t100_m1_rmse.png")
t100_sinus_m2 = get_nsim_estimates(t100_sinus_trees, sinusFun, 2, res_choice=30, "sinus/sinus_t100_m2.png", "sinus/sinus_t100_m2_mae.png", "sinus/sinus_t100_m2_rmse.png")
t100_sinus_m3 = get_nsim_estimates(t100_sinus_trees, sinusFun, 3, res_choice=30, "sinus/sinus_t100_m3.png", "sinus/sinus_t100_m3_mae.png", "sinus/sinus_t100_m3_rmse.png")
# IMPORTANT: Finish skipping here. The lines below will reproduce the plots!

# compare different models for the same function and sample size
t200_sinus_models = compare_cov_prob_models(t200_sinus_m1[[7]], t200_sinus_m2[[7]], t200_sinus_m3[[7]], t200_sinus_m1[[8]], "comp_models/cov_prob/sinus_t200.png")
t100_sinus_models = compare_cov_prob_models(t100_sinus_m1[[7]], t100_sinus_m2[[7]], t100_sinus_m3[[7]], t100_sinus_m1[[8]], "comp_models/cov_prob/sinus_t100.png")

ggarrange(t200_sinus_models, t100_sinus_models, nrow=2, ncol=1, labels=c("n=200","n=100"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
	theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("coverage_plots/comp_models/sinus_models.pdf", units="in", width=7, height=8, dpi=600)

# compare different sample sizes for same function and model
tx_sinus_skykappa = compare_same_model_diff_samp_size(t200_sinus_m1[[7]], t100_sinus_m1[[7]], t200_sinus_m1[[8]], "comp_models/cov_prob/sinus_diff_samp_skykappa.png")
tx_sinus_skygrid = compare_same_model_diff_samp_size(t200_sinus_m2[[7]], t100_sinus_m2[[7]], t200_sinus_m2[[8]], "comp_models/cov_prob/sinus_diff_samp_skygrid.png")
tx_sinus_skygrowth = compare_same_model_diff_samp_size(t200_sinus_m3[[7]], t100_sinus_m3[[7]], t200_sinus_m3[[8]], "comp_models/cov_prob/sinus_diff_samp_skygrowth.png")

ggarrange(tx_sinus_skykappa, tx_sinus_skygrid, tx_sinus_skygrowth, nrow=3, ncol=1, labels=c("Skykappa","Skygrid","Skygrowth"), label.x = 0.05, label.y = 1.05, font.label=list(size=12), legend="right", common.legend = TRUE) +
	theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"))
ggsave("coverage_plots/comp_models/sinus_diff_samp_sizes.pdf", units="in", width=9, height=12, dpi=600)

# compare error considering different models for the same function and sample size
t200_sinus_models_err = compare_err_models(t200_sinus_m1[[9]], t200_sinus_m2[[9]], t200_sinus_m3[[9]], t200_sinus_m1[[10]], t200_sinus_m2[[10]], t200_sinus_m3[[10]], t200_sinus_m1[[8]], "comp_models/mae/sinus_t200.png", "comp_models/rmse/sinus_t200.png")
t100_sinus_models_err = compare_err_models(t100_sinus_m1[[9]], t100_sinus_m2[[9]], t100_sinus_m3[[9]], t100_sinus_m1[[10]], t100_sinus_m2[[10]], t100_sinus_m3[[10]], t100_sinus_m1[[8]], "comp_models/mae/sinus_t100.png", "comp_models/rmse/sinus_t100.png")

# Sina/violin combined plot
tx_sinus_err_models_all <- compare_err_models_combine_all(t200_sinus_models_err[[3]], t100_sinus_models_err[[3]], t200_sinus_models_err[[4]], t100_sinus_models_err[[4]], "comp_models/mae/sinus_diff_models_COMBINED.pdf", "comp_models/rmse/sinus_diff_models_COMBINED.pdf")

# compare error considering different sample sizes for same function and model
tx_sinus_skykappa_err = compare_err_same_model_diff_samp_size(t200_sinus_m1[[9]], t100_sinus_m1[[9]], t200_sinus_m1[[10]], t100_sinus_m1[[10]], t200_sinus_m1[[8]], "comp_models/mae/sinus_diff_samp_skykappa.png", "comp_models/rmse/sinus_diff_samp_skykappa.png")
tx_sinus_skygrid_err = compare_err_same_model_diff_samp_size(t200_sinus_m2[[9]], t100_sinus_m2[[9]], t200_sinus_m2[[10]], t100_sinus_m2[[10]], t200_sinus_m2[[8]], "comp_models/mae/sinus_diff_samp_skygrid.png", "comp_models/rmse/sinus_diff_samp_skygrid.png")
tx_sinus_skygrowth_err = compare_err_same_model_diff_samp_size(t200_sinus_m3[[9]], t100_sinus_m3[[9]], t200_sinus_m3[[10]], t100_sinus_m3[[10]], t200_sinus_m3[[8]], "comp_models/mae/sinus_diff_samp_skygrowth.png", "comp_models/rmse/sinus_diff_samp_skygrowth.png")

# Sina/violin combined plot
tx_sinus_err_all <- compare_err_same_model_diff_samp_size_combine_all(tx_sinus_skykappa_err[[3]], tx_sinus_skygrid_err[[3]], tx_sinus_skygrowth_err[[3]], tx_sinus_skykappa_err[[4]], tx_sinus_skygrid_err[[4]], tx_sinus_skygrowth_err[[4]], "comp_models/mae/sinus_diff_samp_COMBINED.pdf", "comp_models/rmse/sinus_diff_samp_COMBINED.pdf")

end <- Sys.time()
total_time <- as.numeric (end - start, units = "mins")
print(paste("Total time elapsed: ",total_time,"mins"))

summary_sinus_cp = generate_summary_table_cov_prob(sinusFun, t200_sinus_m1[[3]], t200_sinus_m2[[3]], t200_sinus_m3[[3]],
																																																				t100_sinus_m1[[3]],t100_sinus_m2[[3]], t100_sinus_m3[[3]], "sinus_cp.csv")

summary_sinus_rmse = generate_summary_table_rmse(sinusFun, t200_sinus_m1[[10]], t200_sinus_m2[[10]], t200_sinus_m3[[10]],
																																																		t100_sinus_m1[[10]],t100_sinus_m2[[10]], t100_sinus_m3[[10]], "sinus_rmse.csv")

# If running from scratch, the line below will have all the variables generated
system("mkdir -p 500_replicates_simulation/env/")
save.image(file="500_replicates_simulation/env/simuSinus.RData")
saveRDS(t200_sinus_m1, "500_replicates_simulation/env/t200_sinus_m1.rds")

plot_ne_trajectories(t200_sinus_m1, "200", "skykappa", "trajectories_sinus", 100)
plot_ne_trajectories(t200_sinus_m2, "200", "skygrid", "trajectories_sinus", 100)
plot_ne_trajectories(t200_sinus_m3, "200", "skygrowth", "trajectories_sinus", 100)

plot_ne_trajectories_true_ne_without_ci(t200_sinus_m1, "200", "skykappa", "trajectories_sinus_true_ne")
plot_ne_trajectories_true_ne_without_ci(t200_sinus_m2, "200", "skygrid", "trajectories_sinus_true_ne")
plot_ne_trajectories_true_ne_without_ci(t200_sinus_m3, "200", "skygrowth", "trajectories_sinus_true_ne")
