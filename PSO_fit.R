# Routines for fitting the PSO with the Model proposed by Bouzat et al. (2018) and Gasaneo et al. (2019)
# by Richard Schweitzer

# load the hypergeometric function, that is Eq. A1 from Del Punta, Rodriguez, Gasaneo, Bouzat (2019)
# We use one from the Characteristic Functions Toolbox (CharFun):
source('hypergeom1F1.R')


# Eq. 9 from Del Punta, Rodriguez, Gasaneo, Bouzat (2019)
x_t <- function(t, beta, mu, xm, A) { 
  tau <- (((beta+1)*xm) / (A*gamma((beta+mu+1)/mu)))^(1/(beta+1))
  res <- A * ((t^(beta+1)) / (beta+1)) * hypergeom1F1(z = -1*(t/tau)^mu, a = (beta+1)/mu, b = (beta+mu+1)/mu)
  return(res)
}


# Eq. 3 from Del Punta, Rodriguez, Gasaneo, Bouzat (2019)
F_t <- function(t, beta, mu, xm, A) { 
  tau <- (((beta+1)*xm) / (A*gamma((beta+mu+1)/mu)))^(1/(beta+1))
  res <- A * t^beta * exp(-((t^mu)/(tau^mu)))
  return(res)
}


# Eq. 11 from Del Punta et al. (2019), the derivative of F(t)
F_t_deriv <- function(t, beta, mu, xm, A) { 
  tau <- (((beta+1)*xm) / (A*gamma((beta+mu+1)/mu)))^(1/(beta+1))
  res <- A * exp(-(t/tau)^mu) * (beta*t^(beta-1) - (mu/(tau^mu)) * t^(beta+mu-1))
  return(res)
}


# Main model function
PSO_fit <- function(t_sac, xm, beta=1, mu=2, A, gamma0, k0, c=0, d=0, viscosity=FALSE, 
                    use_method="ode45") {
  require(deSolve)
  # define the function to solve the ordinary differential equation, Eq. 2 from Bouzat et al. (2018)
  solvr <- function(t_sac,y,params) {
    # all parameters
    beta <- params[1]
    mu <- params[2]
    tau <- params[3]
    A <- params[4]
    gamma0 <- params[5]
    k0 <- params[6]
    c <- params[7]
    d <- params[8]
    viscosity <- params[9]
    # The force felt by the iris, the derivative of F(t)
    Ft_dt <- A * exp(-(t_sac/tau)^mu) * (beta*t_sac^(beta-1) - (mu/(tau^mu)) * t_sac^(beta+mu-1))  # Eq. 11 from Del Punta et al. (2019)
    # which model: force dependence or inhomogenous viscosity? 
    if (viscosity==FALSE) { # ... the force-dependent model
      # compute gamma and k as force-dependent parameters, as described in Bouzat et al. (2018), Fig. 3
      Ft <- A * t_sac^beta * exp(-((t_sac^mu)/(tau^mu))) # Eq. 3 from Del Punta, Rodriguez, Gasaneo, Bouzat (2019)
      gamma <- gamma0 * exp(-c * Ft) # page 178101-4 in Bouzat et al. (2018)
      k <- k0 * exp(-d * Ft)
      # the forced harmonic oscillator equation, Eq. 10 from Del Punta et al. (2019)
      dy_1 <- y[2]
      dy_2 <- (-1) * gamma * y[2] - k * y[1] - Ft_dt
    } else { # ... the inhomogenous viscosity model
      # k is constant
      k <- k0
      # gamma depends in the distance to the equilibrium position zero
      dy_1 <- y[2]
      dy_2 <- (-1) * (gamma0 * (1 + c*exp(-d * y[1]))) * y[2] - k * y[1] - Ft_dt # page 032422-11 in Del Punta et al. (2019), a/b are now c/d
    }
    y_t <- list(c(dy_1, dy_2))
    return(y_t)
  }
  # assure that values are in sensible ranges
  #print(paste0("xm=", xm, ", beta=", beta, ", mu=", mu, ", A=", A, ", gamma0=", gamma0, ", k0=", k0, ", c=", c, ", d=", d))
  if (any(is.na(c(xm, beta, mu, A, gamma0, k0, c, d)))) {
    stop("At least one of the supplied parameters is NA!")
  }
  if (xm == 0) {
    print("xm was set to 0: Setting xm to 0.001!")
    xm <- 0.001
  } else if (xm < 0) {
    print("xm is negative: xm will now be abs(xm)!")
    xm <- abs(xm)
  }
  if (beta < 1) { 
    print(paste0("Beta is lower than 1: Resetting beta to 1!" ))
    beta <- 1
  }
  if (mu < 0) {
    print("mu is negative: mu will now be abs(mu)!")
    mu <- abs(mu)
  }
  if (A < 0) {
    print("A is negative: A will now be abs(A)!")
    A <- abs(A)
  }
  if (gamma0 < 0) {
    print("gamma0 is negative: gamma0 will now be abs(gamma0)!")
    gamma0 <- abs(gamma0)
  }
  if (k0 < 0) {
    print("k0 is negative: k0 will now be abs(k0)!")
    k0 <- abs(k0)
  }
  if (c < 0) {
    print("c is negative: c will now be abs(c)!")
    c <- abs(c)
  }
  if (d < 0) {
    print("d is negative: d will now be abs(d)!")
    d <- abs(d)
  }
  # compute tau, Eq. 6 from Del Punta, Rodriguez, Gasaneo, Bouzat (2019)
  tau <- (((beta+1)*xm) / (A*gamma((beta+mu+1)/mu)))^(1/(beta+1))
  # Then, get eyeball position over time, Eq. 9 from Del Punta, Rodriguez, Gasaneo, Bouzat (2019)
  x <- A * ((t_sac^(beta+1)) / (beta+1)) * hypergeom1F1(z = -1*(t_sac/tau)^mu, a = (beta+1)/mu, b = (beta+mu+1)/mu)
  # solving numerically for pupil motion
  res <- ode(y = c(0, 0), 
             method = use_method, 
             times = t_sac, 
             func = solvr, 
             parms = c(beta = beta, mu = mu, tau = tau, A = A, gamma0 = gamma0, k0 = k0, c = c, d = d, 
                       viscosity = viscosity)) 
  # pupil motion is eyeball motion plus wobble
  p <- x + res[, 2] # 1st col: time, 2nd col: 1, 3rd col: 2
  return(p)
}

