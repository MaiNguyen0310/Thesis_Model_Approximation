# visusally compare the TAMSD with the ensemble MSD
params <- readRDS("2_Nh50000_/params.rds")
dt_sim <- params$dt_sim
library(ggplot2)
source("moments_system.R")
source("msd_single_particle.R")

single <- simulate_single()
cat("Fraction mobile =", mean(single$state), "\n")
cat("Expected alpha  =", alpha_the, "\n")

lag_time  <- 0.05
lag_index <- round(lag_time/dt_sim)

msd_df <- running_msd(
  single$traj,
  single$time,
  lag_index
)

outdir <- "msd_single"
dir.create(outdir, showWarnings = FALSE)

## Save the concentration fields
saveRDS(
  u_list,
  file = file.path(outdir, "u_list.rds")
)

## Save the trajectory
saveRDS(
  single,
  file = file.path(outdir, "single_particle.rds")
)

## Save the running MSD
saveRDS(
  msd_df,
  file = file.path(outdir, "running_msd.rds")
)


## Save useful simulation settings
saveRDS(
  list(
    dt_sim = dt_sim,
    Tmax = Tmax,
    lag_time = lag_time,
    lag_index = lag_index,
    alpha_the = alpha_the,
    D_phys = D_phys,
    Deff = Deff,
    L = L
  ),
  file = file.path(outdir, "settings.rds")
)



timeseries <- readRDS("2_Nh50000_/timeseries.rds")
idx <- which.min(abs(timeseries$times_all - lag_index*dt_sim))
ensemble_msd <- timeseries$msd_all[idx]



df_plot <- data.frame(
  time = msd_df$time,
  msd  = msd_df$msd
)

p_msd <- ggplot(df_plot, aes(x = time, y = msd)) +
  geom_line(size = 0.85, colour = "black") +
  geom_hline(yintercept = ensemble_msd,
             colour = "#B81840FF",
             linetype = 2,
             linewidth = 0.85) +
  labs(
    x = "Observation time",
    y = paste0("Estimated MSD (lag = ", lag_time, ")"),
    colour = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  scale_colour_manual(values = c("black", "#B81840FF"))

p_msd <- p_msd +
  geom_line(aes(y = NA, colour = "Single-particle estimate", linetype =2)) +
  geom_line(aes(y = NA, colour = "Ensemble MSD"))

p_msd


# plot trajectory
traj_df <- data.frame( x = single$traj[,1], 
                       y = single$traj[,2], 
                       state = factor(single$state) )

ggplot(traj_df, aes(x, y, colour = state)) +
  geom_path() + 
  coord_equal() + 
  theme_minimal()
