library(deSolve)
source("calculate_D.R")
source("definitions_numerics.R")
source("experiment_uniform.R")

set.seed(123)

Tmax = 7.5

run_msd_experiment <- function(N_h, dt_sim = 0.001) {
  run_id <- paste0("MSD_Nh", N_h, "_")
  dir.create(run_id, showWarnings = FALSE, recursive = TRUE)
  
  # --- parameters ---
  L <- 10
  lambda_off <- lambda_off_hat
  mu0 <- mu0_hat
  u0 <- 1              # normalized uniform transporter density
  
  # constant attachment rate in the uniform-u test
  rate_on <- mu0 * u0
  
  # save parameters
  params <- list(
    N_h = N_h,
    Tmax = Tmax,
    dt_sim = dt_sim,
    L = L,
    D_phys = D_phys,
    lambda_off = lambda_off,
    mu0 = mu0,
    u0 = u0,
    rate_on = rate_on
  )
  saveRDS(params, file = file.path(run_id, "params.rds"))
  
  # --- initial conditions ---
  # hitchhikers start near the center
  pos_h <- cbind(
    rnorm(N_h, L/2, 1/2),
    rnorm(N_h, L/2, 1/2)
  )
  
  state_h <- rep(0L, N_h)  # 0 = immotile, 1 = motile
  
  # exponential clocks for detachment
  jump_off <- rexp(N_h, rate = lambda_off)

  # MSD bookkeeping
  pos_h0 <- pos_h
  pos_h_unwrapped <- pos_h
  
  # time series
  n_steps <- ceiling(Tmax / dt_sim) + 1
  time_history <- numeric(n_steps)
  msd_history  <- numeric(n_steps)
  attached_history <- numeric(n_steps)
  
  time_history[1] <- 0
  msd_history[1] <- 0
  attached_history[1] <- mean(state_h)
  
  iter <- 2
  time_sim <- 0
  
  while (time_sim < Tmax) {
    
    # attached particles diffuse
    motile <- which(state_h == 1)
    if (length(motile) > 0) {
      step_res <- step_brownian_reflect(pos_h[motile, ], D_phys, dt_sim, L)
      pos_h[motile, ] <- step_res$pos
      pos_h_unwrapped[motile, ] <- pos_h_unwrapped[motile, ] + step_res$dX
    }
    
    # detachment
    jump_off <- jump_off - dt_sim
    off_now <- which(state_h == 1 & jump_off <= 0)
    if (length(off_now) > 0) {
      state_h[off_now] <- 0L
      jump_off[off_now] <- Inf
    }
    
    # attachment (uniform u)
    immotile <- which(state_h == 0L)
    if (length(immotile) > 0) {
      p_on <- 1 - exp(-rate_on * dt_sim)
      if (p_on > 0) {
        attach_now <- immotile[runif(length(immotile)) < p_on]
        if (length(attach_now) > 0) {
          state_h[attach_now] <- 1L
          jump_off[attach_now] <- rexp(length(attach_now), rate = lambda_off)
        }
      }
    }
    
    time_sim <- time_sim + dt_sim
    
    # MSD from unwrapped positions
    disp <- pos_h_unwrapped - pos_h0
    msd <- mean(rowSums(disp^2))
    
    # store
    if (iter <= n_steps) {
      time_history[iter] <- time_sim
      msd_history[iter] <- msd
      attached_history[iter] <- mean(state_h)
      iter <- iter + 1
    }
    message(paste("time", time_sim, "MSD:", msd))
  }
  
  # trim storage
  time_history <- time_history[1:(iter - 1)]
  msd_history  <- msd_history[1:(iter - 1)]
  attached_history <- attached_history[1:(iter - 1)]
  
  out <- list(
    params = params,
    time = time_history,
    msd = msd_history,
    attached_fraction = attached_history
  )
  
  saveRDS(out, file = file.path(run_id, "msd_timeseries.rds"))
  return(out)
}

run_msd_experiment(10000)
run_msd_experiment(25000)
run_msd_experiment(50000)

