# slightly modified self-avoiding walk model by Engbert et al. (2011) to simulate fixational drift. 
# implemented by Richard Schweitzer

# References:
# Richard Schweitzer and Martin Rolfs. Definition, modeling and detection of saccades in the face of post-saccadic oscillations.
# Ralf Engbert, Konstantin Mergenthaler, Petra Sinn, and Arkady Pikovsky. An integrated model of fixational eye movements and microsaccades. Proceedings of the National Academy of Sciences, 108(39):E765–E770, 2011.


self_avoiding_walk <- function(n_steps=50, # how many steps for the random walk?
                               relax_rate=0.001, # at what rate should the activation decay, i.e. (1-relax_rate)
                               start_ij=c(NaN,NaN), # where in the matrix should we start
                               L=51, # dimensions of the lattice
                               U_slope=1, # the slope of the 2D quadratic potential
                               U_m_sd=c(0.5, 0.2), # mean and SD of Gaussian-random starting values of the activation matrix
                               use_only_direct_neighbors=TRUE, # TRUE: 4 adjacent fields, FALSE: 8 adjacent fields
                               do_final_smooth=FALSE, # TRUE: runs a five-point smoother in the end
                               do_plot=FALSE # want a diagnostic plot?
) {
  require(assertthat)
  assert_that(length(start_ij)==2)
  assert_that(L>3)
  assert_that(relax_rate>0 & relax_rate<1)
  assert_that(length(U_m_sd)==2)
  # instantiate activation matrix
  h <- matrix(data = rnorm(n = L*L, mean = U_m_sd[1], sd = U_m_sd[2]), 
              nrow = L, ncol = L) 
  h[h<0] <- 0
  # instantiate 2D quadratic potential, which is constant
  i_0 <- (L-1)/2
  j_0 <- (L-1)/2
  u <- matrix(data = 0, nrow = L, ncol = L) 
  for (i in 1:L) {
    for (j in 1:L) {
      u[i, j] <- U_slope*L * (((i-i_0)/i_0)^2 + ((j-j_0)/j_0)^2)
    }
  }
  # pre-allocate results
  all_i <- vector(mode = "numeric", length = n_steps)
  all_j <- vector(mode = "numeric", length = n_steps)
  # where do we start?
  if (is.na(start_ij[1])) {
    i <- i_0
  } else {
    i <- start_ij[1]
  }
  if (is.na(start_ij[2])) {
    j <- j_0
  } else {
    j <- start_ij[2]
  }
  # now walk
  for (step_i in 1:n_steps) {
    # increase activation at current position i,j and decay activation on all other sides
    h_current <- h[i,j]
    h <- h * (1 - relax_rate)
    h[i,j] <- h_current + 1
    # determine next step, there are four options, pay attention to borders
    adjacent_activations <- c(NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
    try({adjacent_activations[1] <- h[i+1, j] + u[i+1, j]}, silent = TRUE) # downward
    try({adjacent_activations[2] <- h[i-1, j] + u[i-1, j]}, silent = TRUE) # upward
    try({adjacent_activations[3] <- h[i, j+1] + u[i, j+1]}, silent = TRUE) # rightward
    try({adjacent_activations[4] <- h[i, j-1] + u[i, j-1]}, silent = TRUE) # leftward
    if (use_only_direct_neighbors==FALSE) {
      try({adjacent_activations[5] <- h[i+1, j+1] + u[i+1, j+1]}, silent = TRUE) 
      try({adjacent_activations[6] <- h[i-1, j+1] + u[i-1, j+1]}, silent = TRUE)
      try({adjacent_activations[7] <- h[i-1, j-1] + u[i-1, j-1]}, silent = TRUE) 
      try({adjacent_activations[8] <- h[i+1, j-1] + u[i+1, j-1]}, silent = TRUE)
    }
    which_min <- which(!is.na(adjacent_activations) & 
                         adjacent_activations==min(adjacent_activations, na.rm = TRUE))
    if (is.null(which_min)) { stop(paste0("Could not determine which_min for i=", i, " and j=", j)) }
    if (length(which_min)>1) { # if many have the same value, then choose randomly
      which_min <- sample(which_min)[1]
    }
    # now choose where to go and save that info
    if (which_min==1) {
      i <- i+1
    } else if (which_min==2) {
      i <- i-1
    } else if (which_min==3) {
      j <- j+1
    } else if (which_min==4) {
      j <- j-1
    } else if (which_min==5) {
      i <- i+1
      j <- j+1
    } else if (which_min==6) {
      i <- i-1
      j <- j+1
    } else if (which_min==7) {
      i <- i-1
      j <- j-1
    } else if (which_min==8) {
      i <- i+1
      j <- j-1
    }
    all_i[step_i] <- i
    all_j[step_i] <- j
  }
  return_this <- data.frame(x = all_j, y = all_i)
  if (do_final_smooth) {
    five_point_smoother <- function(x) { # 5-point moving-window smoother
      N <- length(x)
      v <- vector(mode = "numeric", length = N)
      for (r in 1:N) {
        if (r==N || r == 1) { 
          v[r] = x[r]
        } else if (r==N-1 || r == 2) {
          v[r] = (x[r-1] + x[r] + x[r+1]) / 3
        } else {
          v[r] = (x[r-2] + x[r-1] + x[r] + x[r+1] + x[r+2]) / 5
        }
      }
      return(v)
    }
    return_this$x <- five_point_smoother(return_this$x)
    return_this$y <- five_point_smoother(return_this$y)
  }
  # plot
  if (do_plot) {
    require(ggplot2)
    all_ij <- return_this
    all_ij$step <- 1:nrow(all_ij)
    all_ij$x_r <- all_ij$x + rnorm(nrow(all_ij), sd = 0.05) # add some randomness for visualization
    all_ij$y_r <- all_ij$y + rnorm(nrow(all_ij), sd = 0.05)
    p_all_ij <- ggplot(data = all_ij, aes(x = x_r, y = y_r, color = step)) + 
      geom_hline(yintercept = i_0, linetype = "dotted") + 
      geom_vline(xintercept = j_0, linetype = "dotted") + 
      geom_point() + 
      geom_path() + 
      scale_color_viridis_c() + 
      theme_minimal() + 
      labs(x = "X [pix]", y = "Y [pix]", color = "Step")
    print(p_all_ij)
  }
  return(return_this)
}
