# estimate psi
source("definitions_VJ.R")


M <- 20000

# randomly choose velocities
vol_V <- pi * (speed_max^2 - speed_min^2)
angles <- seq(0, 2*pi, length.out = M+1)[-1]
dirs <- cbind(cos(angles), sin(angles))
u <- runif(M)
r <- sqrt(u*(speed_max^2 - speed_min^2) + speed_min^2)
v <- dirs * r



#uniform turning kernel
Tmat <- matrix(1/M, nrow = M, ncol = M)

Tvec <- rep(1/M, M)
vbar <- colSums(v * Tvec)



speeds_sq <- rowSums(v^2)      # ||v||^2 for each discrete velocity
S_T2     <- sum(speeds_sq * Tvec)

S_T2 <- mean(rowSums(v^2))


S_T2
