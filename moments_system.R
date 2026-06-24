 source("calculate_D.R")
 source("calculate_ST2.R")
source("estimate_u0.R")

library(ggplot2)
library(deSolve)

Tmax <- 7.5
L <- 10

ST2 <- 2
#kappa <- D_phys/L^2
kappa <- 0.03878075

# parameters
nu <- lambda_off
lambda <- lambda_tumble/eps^2

mu_c  <- mu0 * u0     
ST2 <- ST2      # from earlier computation (ST2 = 2 in this case)
N0 <- 50000

dt <- 0.01
times <- seq(0, Tmax, by = dt)

moment_system <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
  Dsi2 <- y[["Dsi2"]]
  Dsm2 <- y[["Dsm2"]]
  Bsi  <- y[["Bsi"]]
  Bsm  <- y[["Bsm"]]
  Vsi2 <- y[["Vsi2"]]
  Vsm2 <- y[["Vsm2"]]
  Nsm  <- y[["Nsm"]]
  
  # Nsi is determined by conservation
  Nsi <- N0 - Nsm
  
  # ODEs
  dDsi2 <- nu * Dsm2 - mu_c * Dsi2
  dDsm2 <- 2 * Bsm - dDsi2 
  
  dBsi  <- nu * Bsm - mu_c * Bsi
  dBsm  <- Vsm2 - dBsi - lambda * Bsm
  
  dVsi2 <- nu * Vsm2 - mu_c * Vsi2
  dVsm2 <- -dVsi2 + lambda * ST2 * Nsm -  lambda * Vsm2
  
  dNsm  <- -nu * Nsm + mu_c * Nsi
  
  list(c(dDsi2, dDsm2, 
         dBsi, dBsm, 
         dVsi2, dVsm2, 
         dNsm),
       dDsi2 = dDsi2,
       dDsm2 = dDsm2)
  })
}

# initial conditions 
y0 <- c(
  Dsi2 = 0,
  Dsm2 = 0,
  Bsi  = 0,
  Bsm  = 0,
  Vsi2 = 0,
  Vsm2 = 0,
  Nsm  = 0
)

out <- ode(
  y = y0,
  times = times,
  func = moment_system,
  parms = NULL
)


t_phys <- out[, "time"]

# extract solution
Dsi2 <- out[ , "Dsi2"]
Dsm2 <- out[ , "Dsm2"]

# convert time
t_phys <- out[, "time"] / kappa

# convert MSD
msd_phys <- L^2 * (out[, "Dsi2"] + out[, "Dsm2"]) / N0


df <- data.frame(
  time = out[, "time"],
  msd  = msd_phys,
  alpha_est = L^2 * (out[,"dDsi2"] + out[,"dDsm2"])/ (N0*4*D_phys)
)

p <- ggplot(df, aes(x = time, y = msd)) +
  geom_line(size = 0.9, colour = "black") +
  xlim(0, Tmax) +
  labs(
    x = "Time",
    y = expression("Mean Squared Displacement  " * (MSD))
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(size = 0.3, colour = "grey85"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

p

alpha_est <- L^2 * (out[,"dDsi2"] + out[,"dDsm2"])/ (N0*4*D_phys)


# plot(out[,"time"], alpha_est, type = 'l', lwd =2)
# abline(h=alpha_the)

alpha_the <- (mu_c)/(nu + mu_c)
p2 <- ggplot(df, aes(x = time, y = alpha_est)) +
  geom_line(size = 0.9, colour = "black") +
  labs(
    x = "Time",
    y = expression("Approximation of "~alpha~"based on the slope of the MSD")
  ) +
  xlim(0, Tmax) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(size = 0.3, colour = "grey85"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) + 
  geom_segment(x = 0, 
               xend = Tmax, 
               y=alpha_the,
               yend = alpha_the)+
  annotate(
    "text",
    x = 0, 
    y = alpha_the,
    label = expression(alpha),
    hjust = 1.2,   # move left of the line
    vjust = -0.2,  # slight vertical offset
    size = 5
  )

p2








