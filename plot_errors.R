library(ggplot2)
library(paletteer)

results <- readRDS("results_unif_alpha4.rds")
# Build a tidy data frame
df <- do.call(rbind, lapply(names(results), function(nm) {
  data.frame(
    time   = results[[nm]]$times,
    error  = results[[nm]]$errors,
    system = nm
  )
}))

# Choose a scientific diverging palette
cols <- paletteer_c("ggthemes::Sunset-Sunrise Diverging", length(unique(df$system)))

p <- ggplot(df, aes(x = time, y = error, colour = system)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols) +
  labs(
    x = "Time",
    y = "Relative error",
    colour = "Number of hitchhikers"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

p
