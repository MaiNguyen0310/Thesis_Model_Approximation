# comparative experiment for u intially uniform, s start in the middle
library(fields)


exp_u_uniform <- function(N_h){
  run_id <- paste0("2_Nh", N_h, "_")
  dir.create(run_id, showWarnings = FALSE)
  D <- 1

  # choose a reference density for nondimensionalizing u
  u_ref <- N_t / L^2  
  
  # nondimensional total final time
  Tmax_nd <- kappa * Tmax
  dt_pde <- 0.125*min(dx_nd,dy_nd)^2

  Nt_pde <- ceiling(Tmax_nd / dt_pde)
  
  dt_sim <- dt_pde / kappa
  
  params <- list(
    N_h=N_h, N_t=N_t, L=L,
    lambda_off=lambda_off,
    mu0=mu0,
    D_phys=D_phys,
    dt_sim=dt_sim,
    dt_pde=dt_pde
  )
  saveRDS(params, file=file.path(run_id, "params.rds"))
  ## INITIAL CONDITIONS
  # Transporters: uniform distribution
  pos_t <- cbind(
     runif(N_t,0,L),
     runif(N_t,0,L)
   )
  
  # Hitchhikers simulation: start in middle
  pos_h <- cbind(rnorm(N_h, L/2, 1/2), 
                 rnorm(N_h, L/2, 1/2))
  
  # PDE initial conditions from simulation
  X <- outer(x_nd, rep(1, Ny_grid))
  Y <- outer(rep(1, Nx_grid), y_nd)
  
  # Approximate transporter density from simulation
  u <- bin_to_grid(to_nd_pos(pos_t), 
                    Nx_grid, Ny_grid, dx_nd, dy_nd) / (dx_nd * dy_nd * N_t)
  # hitchhikers PDE
  s_raw <- bin_to_grid(to_nd_pos(pos_h),
                       Nx_grid, Ny_grid, dx_nd, dy_nd)/(dx_nd * dy_nd)
  s <- smooth2d(s_raw)
  
  # store results
  errors <- numeric(length(save_times_phys))
  times  <- numeric(length(save_times_phys))
  alphas <- numeric(length(save_times_phys))
  
  # hitchhiker states
  state_h <- rep(0, N_h)
  
  # initialize rates
  jump_off <- rexp(N_h, rate = lambda_off)   
  
  time_sim <- 0       # physical time
  time_pde <- 0       # scaled time τ = κ t
  
  step <- 1           # index for saving/errors
  
  
  pos_h0 <- pos_h                  # initial positions
  pos_h_unwrapped <- pos_h         # will track true displacement
  pos_prev <- pos_h
  
  msd_vals <- numeric(length(save_times_phys))
  max_steps <- ceiling(Tmax/dt_sim) + 100 # extra room for rounding safety
  
  error_history <- numeric(max_steps)
  time_history  <- numeric(max_steps)
  msd_history   <- numeric(max_steps)
  
  iter <- 1
  
  while(time_sim < Tmax){
    
    ## SIMULATION STEP
    
    motile <- which(state_h == 1)
    # if(length(motile) > 0){
    #   pos_h[motile,] <- step_brownian_reflect(pos_h[motile,], D_phys, dt_sim, L)
    #   # treat reflective boundaries as periodic unwrap
    #   dx <- unwrap_increment(pos_h[motile,1], pos_prev[motile,1], L)
    #   dy <- unwrap_increment(pos_h[motile,2], pos_prev[motile,2], L)
    #   
    #   pos_h_unwrapped[motile,1] <- pos_h_unwrapped[motile,1] + dx
    #   pos_h_unwrapped[motile,2] <- pos_h_unwrapped[motile,2] + dy
    # }
    motile <- which(state_h == 1)
    if(length(motile) > 0){
      step_res <- step_brownian_reflect_track(pos_h[motile,], D_phys, dt_sim, L)
      
      pos_h[motile,] <- step_res$pos
      pos_h_unwrapped[motile,] <- pos_h_unwrapped[motile,] + step_res$dX
    }
   #  pos_prev <- pos_h

    # update status hitchhikers that jump off
    jump_off <- jump_off - dt_sim
    off_now <- which(state_h == 1 & jump_off <= 0)
    if(length(off_now) > 0){
      state_h[off_now] <- 0
      # attached_to[off_now] <- NA
      jump_off[off_now] <- Inf 
    }
    
    alpha <- mu0 * u / (lambda_off + mu0 * u)
    

    for(h in which(state_h == 0)){
      chi_h <- to_nd_pos(matrix(pos_h[h,], nrow = 1))

      local_u_nd <- interp2(u, chi_h[1], chi_h[2], dx_nd, dy_nd)
      local_u_phys <- local_u_nd

      rate_on <- mu0 * local_u_phys

      p <- 1-exp(-rate_on*dt_sim)
      if(runif(1) < p) {
           state_h[h] <- 1
           jump_off[h] <- rexp(1, rate = lambda_off)
      }
    }
    
    time_sim <- time_sim + dt_sim
    
    # ESTIMATE ALPHA
    # total density
    dens_total <- bin_to_grid(
      to_nd_pos(pos_h),
      Nx_grid, Ny_grid, dx_nd, dy_nd
    )
    
    # attached only
    attached_idx <- which(state_h == 1)
    
    dens_attached <- bin_to_grid(
      to_nd_pos(pos_h[attached_idx,,drop=FALSE]),
      Nx_grid, Ny_grid, dx_nd, dy_nd
    )
    
    alpha_emp <- matrix(0, Nx_grid, Ny_grid)
    
    mask <- dens_total > quantile(dens_total, 0.5)
    
    alpha_emp[mask] <- dens_attached[mask] / dens_total[mask]
    alpha_theory <- mu0 * u / (lambda_off + mu0 * u)
    ## ADVANCE PDE
    
    tau_target <- kappa * time_sim

    while(time_pde + dt_pde <= tau_target){
      
    u <- u + dt_pde * laplacian_neumann(u, dx_nd, dy_nd)
      u_phys <- u 
      alpha  <- mu0 * u_phys / (lambda_off + mu0 * u_phys)
      
      q <- alpha * s
      s <- s + dt_pde * laplacian_neumann(q, dx_nd, dy_nd)

      time_pde <- time_pde + dt_pde
    }
    # simulation density 
    s_sim_nd <- bin_to_grid(
      to_nd_pos(pos_h),
      Nx_grid, Ny_grid, dx_nd, dy_nd
    )
    s_sim_density <- smooth2d(s_sim_nd / (dx_nd * dy_nd))
    
    
    # error
    diff <- s_sim_density - s
    E2 <- sqrt(sum(diff^2) * dx_nd * dy_nd)
    E2_rel <- E2 / sqrt(sum(s^2) * dx_nd * dy_nd)
    
    mass_sim <- sum(s_sim_nd)
    mass_pde <- mass2d(s, dx_nd, dy_nd)
    
    # current MSD
    disp <- pos_h_unwrapped - pos_h0
    msd <- mean(rowSums(disp^2))
    
    # save every timestep
    error_history[iter] <- E2_rel
    time_history[iter]  <- time_sim
    msd_history[iter]   <- msd
    
    iter <- iter + 1
    
    ##  SAVE
    
    if(step <= length(save_times_phys) &&
       time_sim >= save_times_phys[step]){
      
      errors[step] <- E2_rel
      times[step]  <- time_sim
      alphas[step] <- mean(alpha)
     
      # calculate and save MSD
      disp <- pos_h_unwrapped - pos_h0
      msd <- mean(rowSums(disp^2))
      
      msd_vals[step] <- msd
      
     
      png(file.path(run_id, sprintf("u_t%.3f.png", time_sim)), width=1200, height=600)
      z_max_u <- quantile(u, 0.995)
      
      image(x_nd, y_nd, t(u),
            col=hcl.colors(200,"YlOrRd"),
            main=paste("PDE (t =", round(time_sim,2),")"))
      image.plot(zlim=c(0,z_max_u),col=hcl.colors(200,"YlOrRd"), legend.only=TRUE)
      dev.off()

      
      pde_plot <- s / N_h
      sim_plot <- s_sim_density / N_h
      
      vals  <- c(pde_plot, sim_plot)
      z_max <- quantile(vals, 0.995)
      
      fname <- sprintf("snapshot_Nh%d_t%.3f_alpha4.png", N_h, time_sim)
      
      png(file.path(run_id, sprintf("snapshot_Nh%d_t%.3f.png", N_h, time_sim)), width=1200, height=600)
      par(mfrow=c(1,2))
      
      image(x_nd, y_nd, t(pde_plot),
            col=hcl.colors(200,"YlOrRd"),
            zlim=c(0,z_max),
            main=paste("PDE density (t =", round(time_sim,2),")"),
            xlab = "[0,L]", xlab = "[0,L]" )
      
      image(x_nd, y_nd, t(sim_plot),
            col=hcl.colors(200,"YlOrRd"),
            zlim=c(0,z_max),
            main=paste("Simulation density (t =", round(time_sim,2),")"),
            xlab = "[0,L]", xlab = "[0,L]")
      
      image.plot(zlim=c(0,z_max), col=hcl.colors(200,"YlOrRd"), legend.only=TRUE)
      
      dev.off()
      
      step <- step + 1
      
      saveRDS(list(
        u=u,
        s=s,
        s_sim=s_sim_density,
        alpha=alpha,
        alpha_emp=alpha_emp
      ), file=file.path(run_id, sprintf("fields_t%.3f.rds", time_sim)))
    }
    cat(
      "t =", round(time_sim,3),
      "mass PDE =", round(mass_pde,2),
      "mass sim =", mass_sim,
      "error =", round(E2_rel,3),
      "mean alpha = ", round(mean(alpha),3),
      "mean attached = ", round(mean(state_h),3),
      "\n"
    )

  }
  saveRDS(list(
    snapshot_times  = times[1:(step-1)],
    snapshot_errors = errors[1:(step-1)],
    snapshot_msd    = msd_vals[1:(step-1)],
    
    times_all  = time_history[1:(iter-1)],
    errors_all = error_history[1:(iter-1)],
    msd_all    = msd_history[1:(iter-1)]
  ), file=file.path(run_id, "timeseries.rds"))
  
  return(list(
    times = time_history[1:(iter-1)],
    errors = error_history[1:(iter-1)],
    msd = msd_history[1:(iter-1)]
  ))
}
