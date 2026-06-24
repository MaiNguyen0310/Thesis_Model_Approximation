# grid definition and numeric functions
# grid
L <- 10                # domain size [0,L] x [0,L]
Nx_grid <- 150
Ny_grid <- 150

dx <- L / (Nx_grid - 1)
dy <- L / (Ny_grid - 1)

x <- seq(0, L, length.out = Nx_grid)
y <- seq(0, L, length.out = Ny_grid)

# nondimensional domain is [0,1]^2
dx_nd <- 1 / (Nx_grid - 1)
dy_nd <- 1 / (Ny_grid - 1)

x_nd <- seq(0, 1, length.out = Nx_grid)
y_nd <- seq(0, 1, length.out = Ny_grid)

to_nd_pos <- function(pos_phys) {
  pos_phys /  L
}



# kappa <- eps^2 * D_phys / L^2
kappa <- D_phys / L^2

#to_nd_pos   <- function(pos_phys) pos_phys * scale_x
to_nd_time   <- function(t_phys)   t_phys * kappa
to_phys_time <- function(t_nd)     t_nd / kappa


# numeric steps
interp2 <- function(M, x, y, dx, dy){
  xi <- x / dx + 1
  yi <- y / dy + 1
  
  i <- floor(xi)
  j <- floor(yi)
  
  i <- max(min(i, nrow(M) - 1), 1)
  j <- max(min(j, ncol(M) - 1), 1)
  
  fx <- xi - i
  fy <- yi - j
  
  (1-fx)*(1-fy)*M[i, j] +
    fx*(1-fy)*M[i+1, j] +
    (1-fx)*fy*M[i, j+1] +
    fx*fy*M[i+1, j+1]
}

# reflect boundaries
reflect_1d <- function(z, L) {
  z <- abs(z)
  z <- z %% (2 * L)
  z[z > L] <- 2 * L - z[z > L]
  z
}

# diffusion - used for transporters and motile hitchhikers
# step_brownian_reflect <- function(pos, D, dt, L) {
#   pos <- as.matrix(pos)
#   if (nrow(pos) == 0){
#     return(pos)
#   }
#   sigma <- sqrt(2 * D * dt)
#   pos <- pos + matrix(rnorm(2 * nrow(pos), 0, sigma), ncol = 2)
#   pos[,1] <- reflect_1d(pos[,1], L)
#   pos[,2] <- reflect_1d(pos[,2], L)
#   pos
# }

step_brownian_reflect <- function(pos, D, dt, L) {
  pos <- as.matrix(pos)
  if (nrow(pos) == 0) {
    return(list(pos = pos, dX = matrix(0, 0, 2)))
  }
  
  sigma <- sqrt(2 * D * dt)
  dX <- matrix(rnorm(2 * nrow(pos), 0, sigma), ncol = 2)
  
  pos_new <- pos + dX
  pos_new[,1] <- reflect_1d(pos_new[,1], L)
  pos_new[,2] <- reflect_1d(pos_new[,2], L)
  
  list(pos = pos_new, dX = dX)
}


# discretesize density
bin_to_grid <- function(pos, Nx, Ny, dx, dy){
  dens <- matrix(0, Nx, Ny)
  if (nrow(pos) == 0) return(dens)
  
  # keep everything inside the last open cell
  x <- pmin(pmax(pos[,1], 0), (Nx - 1) * dx - .Machine$double.eps)
  y <- pmin(pmax(pos[,2], 0), (Ny - 1) * dy - .Machine$double.eps)

  xr <- x / dx
  yr <- y / dy
  
  i <- floor(xr) + 1L
  j <- floor(yr) + 1L
  
  fx <- xr - floor(xr)
  fy <- yr - floor(yr)
  
  for(k in seq_len(nrow(pos))){
    dens[i[k],   j[k]  ] <- dens[i[k],   j[k]  ] + (1-fx[k])*(1-fy[k])
    dens[i[k]+1, j[k]  ] <- dens[i[k]+1, j[k]  ] + fx[k]*(1-fy[k])
    dens[i[k],   j[k]+1] <- dens[i[k],   j[k]+1] + (1-fx[k])*fy[k]
    dens[i[k]+1, j[k]+1] <- dens[i[k]+1, j[k]+1] + fx[k]*fy[k]
  }
  
  dens
}

laplacian <- function(u, dx, dy){

  Nx <- nrow(u)
  Ny <- ncol(u)

  Lu <- matrix(0,Nx,Ny)

  for(i in 2:(Nx-1)){
    for(j in 2:(Ny-1)){

      Lu[i,j] <- (
        4*(u[i+1,j] + u[i-1,j] + u[i,j+1] + u[i,j-1]) +
          (u[i+1,j+1] + u[i+1,j-1] + u[i-1,j+1] + u[i-1,j-1]) -
          20*u[i,j]
      )/(6*dx^2)

    }
  }

  Lu
}

laplacian_neumann <- function(f, dx, dy){
  Nx <- nrow(f)
  Ny <- ncol(f)
  Lf <- matrix(0, Nx, Ny)
  
  for(i in 1:Nx){
    for(j in 1:Ny){
      
      fxm <- if(i == 1) f[i,j] else f[i-1,j]
      fxp <- if(i == Nx) f[i,j] else f[i+1,j]
      
      fym <- if(j == 1) f[i,j] else f[i,j-1]
      fyp <- if(j == Ny) f[i,j] else f[i,j+1]
      
      Lf[i,j] <- (fxp - 2*f[i,j] + fxm)/dx^2 +
        (fyp - 2*f[i,j] + fym)/dy^2
    }
  }
  
  Lf
}

grad <- function(f, dx, dy){
  Nx_grid <- nrow(f); Ny_grid <- ncol(f)
  gx <- matrix(0, Nx_grid, Ny_grid)
  gy <- matrix(0, Nx_grid, Ny_grid)

  for(i in 2:(Nx_grid-1)){
    for(j in 2:(Ny_grid-1)){
      gx[i,j] <- (f[i+1,j] - f[i-1,j]) / (2*dx)
      gy[i,j] <- (f[i,j+1] - f[i,j-1]) / (2*dy)
    }
  }
  return(list(gx = gx, gy = gy))
}

divergence_fv <- function(Fx, Fy, dx, dy){
  Nx <- nrow(Fx) - 1
  Ny <- ncol(Fy) - 1

  div <- matrix(0, Nx, Ny)

  for(i in 1:Nx){
    for(j in 1:Ny){
      div[i,j] <- (Fx[i+1,j] - Fx[i,j]) / dx +
        (Fy[i,j+1] - Fy[i,j]) / dy
    }
  }

  div
}
#
compute_flux <- function(s, alpha, D, dx, dy){
  Nx <- nrow(s)
  Ny <- ncol(s)
  Fx <- matrix(0, Nx+1, Ny) # flux in x-direction (vertical edges)
  Fy <- matrix(0, Nx, Ny+1) # flux in y-direction (horizontal edges)
  # x-fluxes
  for(i in 2:Nx){
    for(j in 1:Ny){
      s_avg <- 0.5 * (s[i,j] + s[i-1,j])
      alpha_avg <- 0.5 * (alpha[i,j] + alpha[i-1,j])
      grad_s <- (s[i,j] - s[i-1,j]) / dx
      grad_alpha <- (alpha[i,j] - alpha[i-1,j]) / dx
      Fx[i,j] <- -D * (alpha_avg * grad_s + grad_alpha * s_avg)
    }
  }
  # y-fluxes
  for(i in 1:Nx){
    for(j in 2:Ny){
      s_avg <- 0.5 * (s[i,j] + s[i,j-1])
      alpha_avg <- 0.5 * (alpha[i,j] + alpha[i,j-1])
      grad_s <- (s[i,j] - s[i,j-1]) / dy
      grad_alpha <- (alpha[i,j] - alpha[i,j-1]) / dy
      Fy[i,j] <- -D * (alpha_avg * grad_s + grad_alpha * s_avg)
    }
  }
  list(Fx=Fx, Fy=Fy)

}

# random direction choice
random_direction <- function() {
  theta <- runif(1, 0, 2*pi)
  c(cos(theta), sin(theta))
}

# smooth2d <- function(M){
#   K <- matrix(1,3,3)/9
#   out <- M
#   for(i in 2:(nrow(M)-1)){
#     for(j in 2:(ncol(M)-1)){
#       out[i,j] <- sum(M[(i-1):(i+1),(j-1):(j+1)] * K)
#     }
#   }
#   out
# }
smooth2d <- function(M){
  K <- matrix(1,3,3)/9
  out <- M
  Nx <- nrow(M)
  Ny <- ncol(M)

  for(i in 1:Nx){
    for(j in 1:Ny){
      ii <- pmax(pmin((i-1):(i+1), Nx), 1)
      jj <- pmax(pmin((j-1):(j+1), Ny), 1)
      out[i,j] <- mean(M[ii, jj])
    }
  }
  out
}


mass2d <- function(M, dx, dy){
  wx <- rep(1, nrow(M))
  wy <- rep(1, ncol(M))
  wx[c(1, nrow(M))] <- 0.5
  wy[c(1, ncol(M))] <- 0.5
  sum(M * outer(wx, wy)) * dx * dy
}


laplacian_neumann <- function(f, dx, dy){
  Nx <- nrow(f)
  Ny <- ncol(f)
  Lf <- matrix(0, Nx, Ny)
  
  for(i in 1:Nx){
    im1 <- if(i == 1) 2 else i - 1
    ip1 <- if(i == Nx) Nx - 1 else i + 1
    
    for(j in 1:Ny){
      jm1 <- if(j == 1) 2 else j - 1
      jp1 <- if(j == Ny) Ny - 1 else j + 1
      
      Lf[i,j] <-
        (f[ip1,j] - 2*f[i,j] + f[im1,j]) / dx^2 +
        (f[i,jp1] - 2*f[i,j] + f[i,jm1]) / dy^2
    }
  }
  
  Lf
}

unwrap_increment <- function(x_new, x_old, L){
  dx <- x_new - x_old
  
  dx[dx >  L/2] <- dx[dx >  L/2] - 2*L
  dx[dx < -L/2] <- dx[dx < -L/2] + 2*L
  
  dx
}
