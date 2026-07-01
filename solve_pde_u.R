# solve carrier PDE seperatelt
source("calculate_D.R")
source("definitions_VJ.R")
source("definitions_numerics2.R")
# Transporters
pos_t <- cbind(
  runif(N_t, 0, L),
  runif(N_t, 0, L)
)
# Initial transporter density
u <- bin_to_grid(
  to_nd_pos(pos_t),
  Nx_grid, Ny_grid,
  dx_nd, dy_nd
) / (dx_nd * dy_nd * N_t)

time_pde <- 0
dt_pde <- params$dt_pde

## Precompute transporter field at all particle times
dt <- 0.1
times <- seq(0, Tmax, by = dt)

u_list <- vector("list", length(times))
u_curr <- u
u_list[[1]] <- u_curr

time_pde <- 0

for(i in 2:length(times)){
  
  tau_target <- kappa * times[i]
  
  while(time_pde + dt_pde <= tau_target){
    
    u_curr <- u_curr +
      dt_pde * laplacian_neumann(u_curr, dx_nd, dy_nd)
    
    time_pde <- time_pde + dt_pde
  }
  
  u_list[[i]] <- u_curr
  
  if(i %% 100 == 0)
    cat(sprintf("Precomputing u: %d / %d\r", i, length(times)))
}