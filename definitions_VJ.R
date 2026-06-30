set.seed(1234)
N_t <- 50000

lambda_off_hat <- 8 # 1/mean time to jump off a transporter

lambda_tumble <- 1/4 # 1/mean run length of a transporter


mu0_hat <- 4
speed_min <- 0
speed_max <- 2

eps <- 0.1

lambda_off <- lambda_off_hat/eps^2
mu0 <- mu0_hat/eps^2

Tmax <- 7.5
