source("definitions_numerics.R")
source("definitions_VJ.R")
source("calculate_D.R")

# PDE uses nondimensional diffusion coefficient 1
D <- 1

D_phys <- 15
# space scaling: chi = (eps/L) x
scale_x <- eps / L

# time scaling: tau = kappa * t
kappa <- D_phys * eps^2 / L^2

# choose a reference density for nondimensionalizing u
u_ref <- N_t / L^2   #me an physical transporter density

# nondimensional parameters
lambda_hat <- lambda_off / kappa
mu_hat     <- mu0 * u_ref / kappa

# nondimensional total final time
Tmax_nd <- kappa * Tmax


# function to run the comparison experiment
run_experiment <- function(N_h){

  # Nt_pde <- Tmax * 150
  dt_pde <- 0.1*min(dx_nd,dy_nd)^2
  # dt_pde <- Tmax_nd / Nt_pde
  Nt_pde <- ceiling(Tmax_nd / dt_pde)

  # dt_sim <- dt_pde * L^2 / (D * eps^2)
  dt_sim <- dt_pde / kappa
  
  ## INITIAL CONDITIONS
  # Transporters 
  sigma_init <- 1   # small radius in physical units
  pos_t <- cbind(
    rnorm(N_t, mean = L/2, sd = sigma_init),
    rnorm(N_t, mean = L/2, sd = sigma_init)
  )
  
  
  # Hitchhikers simulation: uniformly distributed
  pos_h <- cbind(runif(N_h, 0, L), runif(N_h, 0, L))
  
  # PDE initial conditions from simulation
  X <- outer(x_nd, rep(1, Ny_grid))
  Y <- outer(rep(1, Nx_grid), y_nd)
  
  # sigma <- 0.5 * eps 
  # transporters (continuous)
  # u0 <- exp(-((X-eps/2)^2 + (Y-eps/2)^2)/(2*sigma^2))
  # u  <- N_t * u0 / (sum(u0) * dx_nd * dy_nd)
  u <- bin_to_grid(to_nd_pos(pos_t), Nx_grid, Ny_grid, dx_nd, dy_nd) / (dx_nd * dy_nd)
  
  # hitchhikers PDE
  s <- matrix(N_h / eps^2, Nx_grid, Ny_grid)
  # store results
  times <- numeric(Nt_pde)
  errors <- numeric(Nt_pde)
  
  time_sim <- 0
  time_pde <- 0
  
  # initialize transporter directions, speeds, rates 
  dir_t <- t(replicate(N_t, random_direction()))
  speed_t <- runif(N_t, speed_min, speed_max)
  run_t <- rexp(N_t, rate = lambda_tumble)
  
  # hitchhiker states
  state_h <- rep(0, N_h)
  attached_to <- rep(NA, N_h)
  # initialize rates
  jump_on  <- rexp(N_h, rate = 1)
  jump_off <- rexp(N_h, rate = lambda_off)   
  
  for(n in 1:Nt_pde){
    
    
    # PDE step (transporters)
    # Lu <- laplacian(u, dx_nd, dy_nd)
    # u  <- u + dt_pde * D * Lu
    u <- u + dt_pde * D * laplacian_neumann(u, dx_nd, dy_nd)
    
    u_density_nd <- u
    # local multiplication of diffusion coefficent
    #alpha_raw <- mu0 * u_density / (lambda_off + mu0*u_density)
    # u_phys <- (eps / L)^2 * u
    # alpha_raw <- mu0 * u_phys / (lambda_off + mu0 * u_phys)
    # alpha <- alpha_raw
    
    # u_smooth <- smooth2d(u)
    # u_phys <- (eps / L)^2 * u_smooth
    u_phys <- (eps / L)^2 * u
    # alpha_raw <- mu0 * u_phys / (lambda_off + mu0 * u_phys)
    # alpha <- alpha_raw
    alpha <- mu0 * u_phys / (lambda_off + mu0 * u_phys)
  
    q <- alpha * s
    s <- s + dt_pde * D * laplacian_neumann(q, dx_nd, dy_nd)
    
    time_pde <- time_pde + dt_pde
    times[n] <- time_pde
    
    ## simulation during PDE time step
    while(time_sim < time_pde / kappa){
      
      pos_t <- step_brownian_reflect(pos_t, D_phys, dt_sim, L)
      pos_t[,1] <- reflect_1d(pos_t[,1], L)
      pos_t[,2] <- reflect_1d(pos_t[,2], L)
      
      run_t <- run_t - dt_sim
      ## transporter movement for simulations
      # change velocity if run time is over
      tumblers <- which(run_t <= 0)
      if(length(tumblers) > 0){
        dir_t[tumblers,] <- t(replicate(length(tumblers), random_direction()))
        speed_t[tumblers] <- runif(length(tumblers), speed_min, speed_max)
        run_t[tumblers] <- rexp(length(tumblers), rate = lambda_tumble)
      }
      
      ## hitchhiker movement
      motile <- which(state_h == 1)
      
      # movement for the attached hitchhikers
      if(length(motile) > 0){
        pos_h[motile,] <- step_brownian_reflect(
          pos_h[motile,,drop=FALSE],
          D_phys,
          dt_sim,
          L
        )
      }
      
      jump_off <- jump_off - dt_sim
      
      
      for(h in which(state_h == 0)){
        chi_h <- to_nd_pos(matrix(pos_h[h,], nrow = 1))

        local_u_nd <- interp2(u_density_nd, chi_h[1], chi_h[2], dx_nd, dy_nd)
        local_u_phys <- (eps / L)^2 * local_u_nd
        
        rate_on_local <- mu0 * local_u_phys / (lambda_off + mu0 * local_u_phys)
        
        jump_on[h] <- jump_on[h] - rate_on_local * dt_sim
        
        if(jump_on[h] <= 0){
          dists <- sqrt(rowSums((pos_t - pos_h[h,])^2))
          r_attach <- dx
          candidates <- which(dists < r_attach)
          
          if(length(candidates) == 0) next
          
          nearest <- sample(candidates, 1)
          state_h[h] <- 1
          attached_to[h] <- nearest
          pos_h[h,] <- pos_t[nearest,]
          jump_on[h] <- rexp(1, rate = 1)
        }
      }
      time_sim <- time_sim + dt_sim
    }
    
    ## compute E_2(t)
    s_sim_nd <- bin_to_grid(to_nd_pos(pos_h), Nx_grid, Ny_grid, dx_nd, dy_nd)
    
    s_sim_nd <- bin_to_grid(to_nd_pos(pos_h), Nx_grid, Ny_grid, dx_nd, dy_nd)
    s_sim_density <- s_sim_nd / (dx_nd * dy_nd)
    #mass_sim <- mass2d(s_sim_density, dx_nd, dy_nd)
    mass_sim <- sum(s_sim_nd)
    
    mass_pde <- mass2d(s, dx_nd, dy_nd)
    # mass_pde <- sum(s) * dx_nd * dy_nd
    #mass_u   <- mass2d(u, dx_nd, dy_nd)

    diff <- s_sim_density - s
    E2 <- sqrt(sum(diff^2) * dx_nd * dy_nd)
    E2_rel <- E2 / sqrt(sum(s^2) * dx_nd * dy_nd)
    
    # test: mass conv
    message(paste("mass PDE:", mass_pde, "mass sim:", mass_sim ,
                  "Nh:", N_h, "time:", time_sim,
                  "error:" , E2_rel, "\n"))
    errors[n] <- E2_rel
  }
  u_grid <- bin_to_grid(pos_t, Nx_grid, Ny_grid, dx, dy)
  
  # test: visualization of final state
  png("transporter_density.png", width = 600, height = 300)
  image(x_nd, y_nd, t(u / (dx_nd * dy_nd)), main = "Transporter PDE density")
  dev.off()
  
  png(paste0("final_states_",N_h,".png"), width=1200, height=600)
  par(mfrow=c(1,2))
  
  image(
    x = x_nd,
    y = y_nd,
    z = t(s / (dx_nd * dy_nd)),
    main = "PDE density (scaled domain)"
  )
  
  image(
    x = x_nd,
    y = y_nd,
    z = t(s_sim_nd/ (dx_nd * dy_nd)),
    main = "Simulation density (scaled domain)"
  )
  
  dev.off()
  
  #message(paste("finished for N_h =", N_h))
  return(list(times = times, errors = errors))
}
