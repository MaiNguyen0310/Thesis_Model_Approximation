# function and definitions for the single-particle MSD experiment
library(ggplot2)

source("calculate_D.R")
source("definitions_VJ.R")
source("definitions_numerics.R")
params <- readRDS("2_Nh50000_/params.rds")
dt_sim <- params$dt_sim
source("solve_pde_u.R")
source("moments_system.R")



Deff <- alpha_the*D_phys


lag_time <- 0.05      
lag_index <- round(lag_time/dt_sim)


simulate_single <- function(){
  
  
  times <- seq(0,Tmax,by=dt_sim)
  M <- length(times)
  

  
  traj <- matrix(0,M,2)
  
  state <- 0
  states <- numeric(M) 
  states[1] <- state
  
  jump_off <- Inf
  
  pos <- c(
    rnorm(1, L/2, 1/2),
    rnorm(1, L/2, 1/2)
  )
  
  traj[1,] <- pos
  
  for(i in 2:M){
    if(i %% 100 == 0){
      message(sprintf(
        "t = %.3f / %.3f",
        times[i], Tmax
      ))
    }
  
    ## switch on
    
    if(state==0){
      u_curr <- u_list[[i]]
      
      chi <- to_nd_pos(matrix(pos, nrow = 1))
      
      local_u <- interp2(
        u_curr,
        chi[1],
        chi[2],
        dx_nd,
        dy_nd
      )
      
      rate_on <- mu0 * local_u
      
      if(runif(1) < 1-exp(-rate_on*dt_sim)){
        state <- 1
        jump_off <- rexp(1, rate=lambda_off)
      }
      
    }
    
    ## diffuse
    if(state == 1){
      
      step <- step_brownian_reflect(
        matrix(pos, nrow = 1),
        D_phys,
        dt_sim,
        L
      )
      
      pos <- as.numeric(step$pos)
      
      jump_off <- jump_off - dt_sim
      
      if(jump_off <= 0){
        state <- 0
        jump_off <- Inf
      }
      
    }
    
    traj[i,] <- pos
    states[i] <- state
    
  }
  
  list(time=times,
       traj=traj,
       state=states)
  
}

compute_msd <- function(traj, dt_sim, max_frac = 0.1){
  
  M <- nrow(traj)
  
  maxlag <- floor(max_frac*M)
  
  lag <- (1:maxlag)*dt_sim
  msd <- numeric(maxlag)
  
  for(m in 1:maxlag){
    
    dx <- traj[(m+1):M,1] -
      traj[1:(M-m),1]
    
    dy <- traj[(m+1):M,2] -
      traj[1:(M-m),2]
    
    msd[m] <- mean(dx^2 + dy^2)
    
  }
  
  data.frame(
    lag = lag,
    msd = msd
  )
}


running_msd_curve <- function(traj, dt_sim,
                              every = 500,
                              max_frac = 0.1){
  
  M <- nrow(traj)
  
  out <- list()
  
  k <- 1
  
  for(i in seq(1000, M, by = every)){
    
    message(sprintf("MSD curve using trajectory up to t = %.2f",
                    (i-1)*dt))
    
    out[[k]] <- compute_msd(
      traj[1:i,],
      dt_sim,
      max_frac = max_frac
    )
    
    out[[k]]$time <- (i-1)*dt_sim
    
    k <- k+1
  }
  
  out
}




single <- simulate_single()
cat("Fraction mobile =", mean(single$state), "\n")
cat("Expected alpha  =", alpha_the, "\n")

running_msd <- function(traj, times, lag_index){
  
  M <- nrow(traj)
  
  msd <- rep(NA, M)
  
  for(i in (lag_index+2):M){
    
    if(i %% 500 == 0)
      message(sprintf("Running MSD: %.2f / %.2f",
                      times[i], Tmax))
    
    dx <- traj[(lag_index+1):i,1] -
      traj[1:(i-lag_index),1]
    
    dy <- traj[(lag_index+1):i,2] -
      traj[1:(i-lag_index),2]
    
    msd[i] <- mean(dx^2 + dy^2)
    
  }
  
  data.frame(
    time = times,
    msd = msd
  )
}

lag_time  <- 0.05
lag_index <- round(lag_time/dt_sim)

msd_df <- running_msd(
  single$traj,
  single$time,
  lag_index
)

timeseries <- readRDS("2_Nh50000_/timeseries.rds")
idx <- which.min(abs(timeseries$times_all - lag_index*dt_sim))
ensemble_msd <- timeseries$msd_all[idx]