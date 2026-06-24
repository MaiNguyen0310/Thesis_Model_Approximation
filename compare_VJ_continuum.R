source("definitions_VJ.R")
source("calculate_D.R")
source("definitions_numerics.R")
source("run_experiment.R")

library("paletteer")

#D <- D_scalar          # effective diffusion coefficient for transporters
D <- 1
mu0 <- 1            # attachment rate scale 
 
Nh_list <- c(50000)

results <- lapply(Nh_list, run_experiment)
names(results) <- paste0("Nh=", Nh_list)

png("errors_vs_time.png", width=800, height=600)
plot(NULL, xlim=c(0,max(results[[1]]$times)), ylim=c(0, max(sapply(results, function(r) max(r$errors)))),
     xlab="Time", ylab="Relative error")
cols <- paletteer::paletteer_c("grDevices::Temps", length(results))
for(i in 1:length(results)){
  lines(results[[i]]$times, results[[i]]$errors, col=cols[i], lwd=0.8)
}
legend("bottomright", legend=names(results), col=cols, lwd=2, cex=0.65)
dev.off()

final_errors <- sapply(results, function(r) tail(r$errors,1))

png("final_errors.png", width=800, height=600)
plot(final_errors, type="b", log="x",
     xlab="Number of hitchhikers", ylab="Final relative L2 error")
dev.off()

saveRDS(results, file="results.rds")
saveRDS(final_errors, file="final_errors.rds")