# visually compare the MSD for the simulations and MSD system
# panel figure
source("msd_experiment.R")
source("moments_system.R")

library(ggplot2)
library(grid)


results_msd_50000 <- readRDS("MSD_Nh50000_/msd_timeseries.rds")
df_sim_50000 <- data.frame(results_msd_50000$time, results_msd_50000$msd)

results_msd_25000 <- readRDS("MSD_Nh25000_/msd_timeseries.rds")
df_sim_25000 <- data.frame(results_msd_25000$time, results_msd_25000$msd)

results_msd_10000 <- readRDS("MSD_Nh10000_/msd_timeseries.rds")
df_sim_10000 <- data.frame(results_msd_10000$time, results_msd_10000$msd)


names(df_sim_50000) <- names(df_sim_25000) <- names(df_sim_10000) <- c("time", "msd")


df_50000 <- df_sim_50000
df_50000$type <- "Simulation (50000 Hitchhikers)"

df_25000 <- df_sim_25000
df_25000$type <- "Simulation (25000 Hitchhikers)"

df_10000 <- df_sim_10000
df_10000$type <- "Simulation (10000 Hitchhikers)"


df2 <- df
df2$type <- "Closed moments system"

df_all <- rbind(df_50000[, c("time","msd","type")],
                df_25000[, c("time","msd","type")],
                df_10000[, c("time","msd","type")],
                 df2[, c("time","msd","type")])


# [0, 0.5]
p_zoom1 <- ggplot(subset(df_all, time >= 0 & time <= 1),
                  aes(time, msd, colour = type)) +
  geom_line(size = 0.7) +
  scale_colour_manual(values = c("black", "#33608CFF", "#B81840FF", "#F6BA57FF")) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.background = element_rect(colour = "black", fill = "white", size = 0.6)
  )+
  labs(x = "Time", y = "MSD") 
g1 <- ggplotGrob(p_zoom1)

# [4,5]
p_zoom2 <- ggplot(subset(df_all, time >= 4 & time <= 5),
                  aes(time, msd, colour = type)) +
  geom_line(size = 0.7) +
  scale_colour_manual(values = c("black", "#33608CFF", "#B81840FF", "#F6BA57FF")) +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.background = element_rect(colour = "black", fill = "white", size = 0.6)
  )+
  labs(x = "Time", y = "MSD") 
g2 <- ggplotGrob(p_zoom2)

# Main plot
p_main <- ggplot(df_all, aes(time, msd, colour = type)) +
  geom_line(size = 0.85) +
  scale_colour_manual(values = c("black", "#33608CFF", "#B81840FF", "#F6BA57FF")) +
  xlim(0, Tmax) +
  labs(x = "Time", y = "Mean Squared Displacement", colour = "") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

# Top zoom box rectangle 
xmin1 <- 0.05 * Tmax
xmax1 <- 0.45 * Tmax
ymin1 <- 0.58 * ymax
ymax1 <- 0.99 * ymax  

# Bottom zoom box rectangle 
xmin2 <- 0.55 * Tmax
xmax2 <- 0.95 * Tmax
ymin2 <- 0.002 * ymax
ymax2 <- 0.4 * ymax

# Final plot
p_final <- p_main +
  
annotation_custom(
  grob = g1,
  xmin = xmin1, xmax = xmax1,
  ymin = ymin1, ymax = ymax1
) +
  
  # diagonal connectors 
  geom_segment(
    aes(
      x = 0.0,
      y = df_all$msd[which.min(abs(df_all$time - 0))],
      xend = xmin1,
      yend = ymin1
    ),
    linewidth = 0.4, linetype = 2
  ) +
  geom_segment(
    aes(
      x = 1.0,
      y = df_all$msd[which.min(abs(df_all$time - 0.5))],
      xend = xmax1,
      yend = ymin1
    ),
    linewidth = 0.4, linetype = 2
  ) +
  
  #  BOTTOM ZOOM BOX 
annotation_custom(
  grob = g2,
  xmin = xmin2, xmax = xmax2,
  ymin = ymin2, ymax = ymax2
) +
  
  # diagonal connectors 
  geom_segment(
    aes(
      x = 5.5,
      y = df_all$msd[which.min(abs(df_all$time - 4))],
      xend = xmin2,
      yend = ymax2    
    ),
    linewidth = 0.4, linetype = 2
  ) +
  geom_segment(
    aes(
      x = 6,5,
      y = df_all$msd[which.min(abs(df_all$time - 5))],
      xend = xmax2,
      yend = ymax2     #  upper-right corner
    ),
    linewidth = 0.4, linetype = 2
  )

p_final
