#' Confluent hypergeometric function
#' from the Characteristic Functions Toolbox (CharFun): 
#' https://www.rdocumentation.org/packages/CharFun/versions/0.1.0
#'
#' HYPERGEOM1F1 Computes the confluent hypergeometric function 1F1(a,b,z),
#' also known as the Kummer function M(a,b,z), for the real parameters a
#' and b (here assumed to be scalars), and the complex argument z
#' (could be scalar, vector or array).
#'
#' The algorithm is based on a Fortran program in S. Zhang & J. Jin
#' "Computation of Special Functions" (Wiley, 1996).
#' Converted to matlab by using f2matlab:
#' https://sourceforge.net/projects/f2matlab/
#' written by Ben Barrowes.
#'
#'
#' z numerical values (number, vector...)
#' a single value
#' b single value
#' hypergeometric function 1F1(a,b,z)

hypergeom1F1 <- function(z, a, b) {
  sz <- dim(z)
  z <- c(z)
  f <- unlist(lapply(z, function(z) cchg(z, a, b)))
  return(f)
}

cchg <- function (z, a, b) {
  chw <- 0
  pi <- 3.141592653589793
  ci <- 1i
  a0 <- a
  # print(paste("Running cchg with:", z, a, b))
  if (!is.na(a) & !is.na(b) & !is.na(z)) {
    if (b == 0 || b == -trunc(abs(b))) {
      chg <- 1e+300
    } else if (a == 0 || z == 0) {
      chg <- 1
    } else if (a == -1) {
      chg <- 1 - z / b
    } else if (a == b) {
      chg <- exp(z)
    } else if (a - b == 1) {
      chg <- (1 + z / b) * exp(z)
    } else if (a == 1 && b == 2) {
      chg <- (exp(z) - 1) / z
    } else if (a == trunc(a) && a < 0) {
      m   <- trunc(-a)
      cr  <- 1
      chg <- 1
      for (k in (1:m)) {
        cr  <- cr * (a + k - 1) / k / (b + k - 1) * z
        chg <- chg + cr
      }
    } else {
      x0 <- Re(z)
      if (x0 < 0) {
        a  <- b - a
        a0 <- a
        z  <- -z
      }
      if (a < 2) {nl <- 0}
      if (a >= 2) {
        nl <- 1
        la <- trunc(a)
        a  <- a - la - 1
      }
      for (n in (0:nl)) {
        if (a0 >= 2) {a <- a + 1}
        if (abs(z) < 20 + abs(b) || a < 0) {
          chg <- 1
          crg <- 1
          for (j in (1:500)) {
            crg <- crg * (a + j - 1) / (j * (b + j - 1)) * z
            chg <- chg + crg
            if (abs((chg - chw) / chg) < 1e-15) {break}
            chw = chg
          }
        } else {
          g1 <- gamma(a)
          g2 <- gamma(b)
          ba <- b - a
          g3 <- gamma(ba)
          cs1 <- 1
          cs2 <- 1
          cr1 <- 1
          cr2 <- 1
          for (i in (1:8)) {
            cr1 <- -cr1 * (a + i - 1) * (a - b + i) / (z * i)
            cr2 <-  cr2 * (b - a + i - 1) * (i - a) / (z * i)
            cs1 <-  cs1 + cr1
            cs2 <-  cs2 + cr2
          }
          x <- Re(z)
          y <- Im(z)
          if (x == 0 && y >= 0) {phi <- 0.5 * pi
          } else if (x == 0 && y <= 0) {phi <- -0.5 * pi
          } else {phi <- atan(y / x)}
          if (phi > -0.5 * pi && phi < 1.5 * pi) {ns <- 1}
          if (phi > -1.5 * pi && phi <= -0.5 * pi) {ns <- -1}
          cfac = exp(ns * ci * pi * a)
          if (y == 0) {cfac = cos(pi * a)}
          chg1 <- g2 / g3 * z ^ (-a) * cfac * cs1
          chg2 <- g2 / g1 * exp(z) * z ^ (a - b) * cs2
          chg <- chg1 + chg2
        }
        if (n == 0) {cy0 <- chg}
        if (n == 1) {cy1 <- chg}
      }
      if (a0 >= 2) {
        for (i in (1:(la - 1))){
          chg <- ((2 * a - b + z) * cy1 + (b - a) * cy0) / a
          cy0 <- cy1
          cy1 <- chg
          a   <- a + 1
        }
      }
      if (x0 < 0) {chg <- chg * exp(-z)}
    }
  } else {
    chg <- NaN
  }
  return(chg)
}

GAMMA <- function(x, z) {
  pi = 3.141592653589793
  if (x == trunc(x)) {
    if (x > 0) {
      ga = 1
      m1 = x - 1
      for  (k in seq(2, m1, max(0, length.out = m1-2))) {
        ga = ga * k
      }
    } else {ga = 1e+300}
  } else {
    if (abs(x) > 1) {
      z = abs(x)
      m = trunc(z)
      r = 1
      for  (k in (1:m)) {r = r * (z - k)}
      z = z - m
    } else
      z = x
  }
  g <- c(1.0e0,0.5772156649015329e0,-0.6558780715202538e0,
         -0.420026350340952e-1,0.1665386113822915e0,-0.421977345555443e-1,
         -0.96219715278770e-2,0.72189432466630e-2,-0.11651675918591e-2,
         -0.2152416741149e-3,0.1280502823882e-3,-0.201348547807e-4,
         -0.12504934821e-5,0.11330272320e-5,-0.2056338417e-6,0.61160950e-8,
         0.50020075e-8,-0.11812746e-8,0.1043427e-9,0.77823e-11,-0.36968e-11,
         0.51e-12,-0.206e-13,-0.54e-14,0.14e-14,0.1e-15)
  gr <- g[26]
  for (k in (25:1)) {gr = gr * z + g[k]}
  ga = 1 / (gr * z)
  if (abs(x)> 1) {ga = ga * r
  if(x < 0) {ga = -pi / (x * ga * sin(pi * x))}
  }
  return(ga)
}