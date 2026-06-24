source("definitions_numerics.R")
set.seed(123)
pos_t <- cbind(
  runif(N_t,0,L),
  runif(N_t,0,L)
)
u <- bin_to_grid(to_nd_pos(pos_t), 
                 Nx_grid, Ny_grid, dx_nd, dy_nd) / (dx_nd * dy_nd * N_t)
u0 <- mean(u)
