# Clean all variables from environment
rm(list=ls())

library(glue)
library(data.table)
library(ggplot2)
library(ggforce)

options(scipen=1, digits=8)

#load(file="500_replicates_simulation/env/simuResTauTablesSinus.RData")
source("500_replicates_simulation/resTauCommon.R")

sinusFun=function(x){sin(x)*10+12}

# To use same simulated trees as in other analysis, load `t200_sinus_trees` (using n=200 samples from sinus function)
t200_sinus_trees <- readRDS(file="500_replicates_simulation/env/t200_sinus_trees.rds")
t200_sinus_m1 <- readRDS(file="500_replicates_simulation/env/t200_sinus_m1.rds") # tau-optimised simulations

tau_choices <- c(0.001, 1, 20, 1000)

t200_sinus_m1_both_fix <- list()
for(l in 1:length(tau_choices)) {
	print(glue("tau={tau_choices[l]}==="))
	t200_sinus_m1_both_fix[[l]] = get_nsim_estimates2(t200_sinus_trees, sampleDates_200, sinusFun, 1, which_fixed="both", res_val=30, tau_val=tau_choices[l])
}

#saveRDS(t200_sinus_m1_both_fix, file = "500_replicates_simulation/env/t200_sinus_m1_both_fix.rds")

# plot RMSE varying tau = {0.001, 1, 20, 1000}
t200_sinus_m1_both_fix[[1]][[2]]$tau <- "0.001"
t200_sinus_m1_both_fix[[2]][[2]]$tau <- "1"
t200_sinus_m1_both_fix[[3]][[2]]$tau <- "20"
t200_sinus_m1_both_fix[[4]][[2]]$tau <- "1000"
t200_sinus_m1[[10]]$tau <- "Optimized"
rmse_t200_sinus_m1_both_fix <- rbind(t200_sinus_m1[[10]], t200_sinus_m1_both_fix[[1]][[2]], t200_sinus_m1_both_fix[[2]][[2]], t200_sinus_m1_both_fix[[3]][[2]], t200_sinus_m1_both_fix[[4]][[2]])
rmse_t200_sinus_m1_both_fix$tau <- factor(rmse_t200_sinus_m1_both_fix$tau, levels = unique(rmse_t200_sinus_m1_both_fix$tau))

# p1 <- ggplot(data=rmse_t200_sinus_m1_both_fix, aes(x=as.numeric(as.character(n_sim)), y=rmse_val, color=tau)) + 
# 	geom_line() + geom_point() + #color="steelblue" 
# 	scale_color_manual(name = "tau", values = c("#a64d79ff", "steelblue", "darkred", "darkolivegreen", "#e69138ff")) +
# 	labs(x = "Simulation index", y = "Root Mean Square Error (RMSE)", caption = paste0("Overall RMSE = ", scaleFUN(rmse_t200_sinus_m1_both_fix$rmse_val))) +
# 	coord_cartesian(ylim=c(0,20)) +
# 	theme_bw() + theme(plot.title = element_text(hjust=0.5), plot.caption = element_text(hjust=1))

p2 <- ggplot(rmse_t200_sinus_m1_both_fix, aes(x=tau, y=rmse_val, color=tau)) + 
	geom_violin() + geom_sina(size=0.3) +
	scale_color_manual(name = "tau", values = c("#a64d79ff", "steelblue", "darkred", "darkolivegreen", "#e69138ff")) +
	labs(x="tau", y="RMSE across simulations") + coord_cartesian(ylim=c(0,20)) +
	theme_bw() + theme(plot.title = element_text(hjust=0.5), axis.text=element_text(size=7), axis.title=element_text(size=7), legend.position="none")
ggsave("error_plots/comp_models/rmse/sinus_optTau_vs_fixed.pdf", plot=p2, units="in", width=6, height=4.5, dpi=600)

system("mkdir -p 500_replicates_simulation/env/")
save.image(file="500_replicates_simulation/env/simuResTauSinus.RData")