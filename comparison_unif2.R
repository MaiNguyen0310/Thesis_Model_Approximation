# code to run the hybrid simulation of the markovian hitchhikers' model 
## with a homogeneous carrier population.
# the relative error over time is plotted and saved.

source("definitions_VJ.R")
source("calculate_D.R")
source("definitions_numerics.R")
source("experiment_uniform.R")

library("paletteer")

# times to save the density plots
save_times_phys <- c(0.001,0.5,1,
                     1.5,2,
                     2.5,3,
                     3.5,4,
                     4.5,5,
                     5.5,6,
                     6, 6.5,
                     7, 7.5)   # physical time


save_times <- save_times_phys * kappa    # convert to nd time

# number of hitchhikers
Nh_list <- c(50000,25000,10000)

results <- lapply(Nh_list, exp_u_uniform)
names(results) <- paste0("Nh=", Nh_list)

# save plot of relative error over time
png("errors_vs_time_unif2_alpha4.png", width=800, height=600)
plot(NULL, xlim=c(0,max(results[[1]]$times)), ylim=c(0, max(sapply(results, function(r) max(r$errors)))),
     xlab="Time", ylab="Relative error")
cols <- paletteer::paletteer_c("ggthemes::Sunset-Sunrise Diverging", length(results))
for(i in 1:length(results)){
  lines(results[[i]]$times, results[[i]]$errors, col=cols[i], lwd=0.8)
}
legend("bottomright", legend=names(results), col=cols, lwd=2, cex=0.85)
dev.off()

# save results for later plot adjustments
saveRDS(results, file="results_unif_alpha4.rds")
saveRDS(final_errors, file="final_errors_unif_alpha4.rds")