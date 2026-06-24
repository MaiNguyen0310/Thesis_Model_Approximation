# code to estimate diffusion coefficient D
source("definitions_VJ.R")
library(MASS)   # for ginv()

K <- 2000               # number of discrete velocities

# Volume of velocity space 
vol_V <- pi * (speed_max^2 - speed_min^2)

angles <- seq(0, 2*pi, length.out = K+1)[-1]

# unit directions
dirs <- cbind(cos(angles), sin(angles))

# choose speeds uniformly in [speed_min, speed_max]
r <- sqrt(seq(speed_min^2, speed_max^2, length.out = K))

# actual velocities v_i
v <- dirs * r

#uniform turning kernel
Tmat <- matrix(1/K, nrow = K, ncol = K)

# L = -lambda I + lambda T
Lmat <- -lambda_tumble * diag(K) + lambda_tumble * Tmat

# COMPUTE THE PSEUDOINVERSE F
Fmat <- ginv(Lmat)   # Moore–Penrose pseudoinverse

Dmat <- matrix(0, 2, 2)

for(i in 1:K){
  for(j in 1:K){
    Dmat <- Dmat + Fmat[i,j] * (v[i,] %*% t(v[j,]))
  }
}

Dmat <- - Dmat / K

# If isotropic, scalar diffusion coefficient:
D_phys <- mean(diag(Dmat))