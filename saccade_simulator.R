# Saccade simulator, presented in the book chapter "Definition, modeling and detection of saccades in the face of post-saccadic oscillations."
# implemented by Richard Schweitzer
saccade_simulator <- function(fixpos_xy, # a matrix or data.frame with 3 columns and >=2 rows, the latter each representing one fixation pos and dur
                              sac_params, # same as above, only containing 8 columns and >=1 rows (one less than fixations)
                              fix_params, # also a matrix or data.frame with X columns and >=2 rows, as many as fixations exist
                              pix_per_dva, # how many pixels per degree of visual angle?
                              detect_vel_criterion = 5, # ground-truth velocity criterion [dva/s]
                              sampling_rate = 1000, # what's the sampling rate we want to achieve?
                              noise_to_saccades_sd = 0, # that's the gaussian noise SD, that can be added to the saccade
                              smooth_to_fixation = TRUE, # smoothes the fixation period to reduce inter-sample noise
                              do_plot = FALSE # a plot of the results?
) {
  require(assertthat)
  require(data.table)
  # first, make sure of a few things
  assert_that((is.matrix(fixpos_xy) || is.data.frame(fixpos_xy)) & 
                nrow(fixpos_xy)>=2 & ncol(fixpos_xy)==3, 
              msg = "fixpos_xy must have >=2 rows and 3 columns: fixation X, Y [pix], and duration [ms]")
  assert_that((is.matrix(sac_params) || is.data.frame(sac_params)) & 
                nrow(sac_params)>=1 & ncol(sac_params)==8, 
              msg = "sac_params must have >=1 rows and 8 columns: beta, mu, A, gamma, k, c, d, viscosity")
  assert_that((is.matrix(fix_params) || is.data.frame(fix_params)) & 
                nrow(fix_params)>=2 & ncol(fix_params)==5, 
              msg = "fix_params must have >=2 rows and 5 columns: relax_rate, U_slope, U_mean, U_sd, only_direct")
  assert_that(nrow(fixpos_xy)==nrow(fix_params) & nrow(fix_params)==nrow(sac_params)+1, 
              msg = "fixpos_xy and fix_params must have the same number of rows, and sac_params one less!")
  # second, get important parameters for this run
  fix_lattice_size <- 2 * round(pix_per_dva) + 1 # two degrees of visual angle (from Engbert et al., 2011, PNAS)
  fix_lattice_center <- (fix_lattice_size-1)/2
  inter_sample_dur <- (1000/sampling_rate) # in milliseconds
  # new prepare loop
  current_time = 0
  resulting_df <- vector(mode = "list", length = nrow(fixpos_xy)) # resulting data.frame
  resulting_metrics <- vector(mode = "list", length = nrow(fixpos_xy)) # resulting saccade parameters
  # okay, now run through fixations
  current_xy = c(fixpos_xy[1, 1:2])
  for (fix_nr in 1:nrow(fixpos_xy)) {
    # 1. simulate fixation data
    fixation_steps <- fixpos_xy[fix_nr, 3] / inter_sample_dur + 1
    this_fixation <- self_avoiding_walk(n_steps = fixation_steps, L = fix_lattice_size, 
                                        start_ij = c(fix_lattice_center, fix_lattice_center),
                                        relax_rate = fix_params[fix_nr, 1], U_slope = fix_params[fix_nr, 2], 
                                        U_m_sd = c(fix_params[fix_nr, 3], fix_params[fix_nr, 4]), 
                                        use_only_direct_neighbors = as.logical(fix_params[fix_nr, 5]), 
                                        do_final_smooth = smooth_to_fixation) # NEW: option to smooth
    # enrich simulation results
    this_fixation$x <- this_fixation$x - fix_lattice_center + current_xy[1]
    this_fixation$y <- this_fixation$y - fix_lattice_center + current_xy[2]
    this_fixation$time <- seq(from = current_time, to = current_time+fixpos_xy[fix_nr, 3], 
                              by = inter_sample_dur)
    this_fixation$what <- "FIX"
    # update current time, x and y, so that our saccade can start from the latest fixation position
    time_fixon <- current_time
    x_fixon <- current_xy[1]
    y_fixon <- current_xy[2]
    current_time <- this_fixation$time[length(this_fixation$time)] + inter_sample_dur
    current_xy <- as.numeric(this_fixation[length(this_fixation$time), c("x", "y")]) 
    time_fixoff <- current_time - inter_sample_dur
    # 2. simulate saccade data (if present)
    if (fix_nr <= nrow(sac_params) & fix_nr < nrow(fix_params)) {
      # get saccade amplitude
      sac_amp <- xy_to_dist(vx = fixpos_xy[fix_nr+1, 1]-current_xy[1], 
                            vy = fixpos_xy[fix_nr+1, 2]-current_xy[2]) / pix_per_dva
      sac_dir <- xy_to_degree(vx = fixpos_xy[fix_nr+1, 1]-current_xy[1], 
                              vy = fixpos_xy[fix_nr+1, 2]-current_xy[2])
      # simulate saccade according to parameters
      # 1. use extremely high temporal res to be sure to be able to detect ground truth
      sac_time_full <- seq(from = 0, to = 300, by = 0.1)
      sac_prof_full <- PSO_fit(t_sac = sac_time_full, xm = sac_amp, 
                               beta = sac_params[fix_nr, 1], mu = sac_params[fix_nr, 2], A = sac_params[fix_nr, 3], 
                               gamma0 = sac_params[fix_nr, 4], k0 = sac_params[fix_nr, 5], 
                               c = sac_params[fix_nr, 6], d = sac_params[fix_nr, 7], 
                               viscosity = as.logical(sac_params[fix_nr, 8]) )
      # extract saccade metrics, like duration etc
      full_vel_vec <- get_vel_vec(sac_prof_full, sac_time_full)
      sac_vpeak <- max(abs(full_vel_vec)[!is.na(full_vel_vec) & !is.infinite(full_vel_vec) & 
                                           abs(full_vel_vec)>=detect_vel_criterion], na.rm = TRUE)
      sac_vpeak_when <- sac_time_full[!is.na(full_vel_vec) & !is.infinite(full_vel_vec) & 
                                        sac_vpeak==full_vel_vec]
      sac_vmean <- mean(abs(full_vel_vec)[!is.na(full_vel_vec) & !is.infinite(full_vel_vec) & 
                                            abs(full_vel_vec)>=detect_vel_criterion], na.rm = TRUE)
      full_dur_PSO <- sac_time_full[min(which(!is.na(full_vel_vec) & !is.infinite(full_vel_vec) & 
                                                full_vel_vec < detect_vel_criterion & sac_time_full > sac_vpeak_when), 
                                        na.rm = TRUE)]
      if (!is.na(full_dur_PSO) && 
          any(abs(full_vel_vec[!is.na(full_vel_vec) & !is.infinite(full_vel_vec) & 
                               sac_time_full>full_dur_PSO]) > detect_vel_criterion)) { # is the threshold passed once more?
        # get last time point where the velocity threshold was crossed ...
        time_last_above_threshold <- sac_time_full[max(which(!is.na(full_vel_vec) & !is.infinite(full_vel_vec) & 
                                                               abs(full_vel_vec) > detect_vel_criterion))] 
        # ... and get first time point after that below threshold
        full_dur_Full <- sac_time_full[min(which(!is.na(full_vel_vec) & !is.infinite(full_vel_vec) & 
                                                   full_vel_vec < detect_vel_criterion & sac_time_full > time_last_above_threshold))]
      } else { 
        full_dur_Full <- full_dur_PSO
      }
      time_sacon <- current_time
      x_sacon <- current_xy[1]
      y_sacon <- current_xy[2]
      time_sacoff <- full_dur_Full + current_time
      time_sacoff_PSO <- full_dur_PSO + current_time
      # 2. use specified temporal res
      sac_time <- seq(from = 0, to = ceiling(full_dur_Full), by = inter_sample_dur)
      sac_prof <- PSO_fit(t_sac = sac_time, xm = sac_amp, 
                          beta = sac_params[fix_nr, 1], mu = sac_params[fix_nr, 2], A = sac_params[fix_nr, 3], 
                          gamma0 = sac_params[fix_nr, 4], k0 = sac_params[fix_nr, 5], 
                          c = sac_params[fix_nr, 6], d = sac_params[fix_nr, 7], 
                          viscosity = as.logical(sac_params[fix_nr, 8]) )
      sac_x <- degree_to_xy(sac_prof, sac_dir)[[1]] * pix_per_dva + current_xy[1]
      sac_x <- sac_x + rnorm(n = length(sac_x), mean = 0, sd = noise_to_saccades_sd)
      sac_y <- degree_to_xy(sac_prof, sac_dir)[[2]] * pix_per_dva + current_xy[2]
      sac_y <- sac_y + rnorm(n = length(sac_y), mean = 0, sd = noise_to_saccades_sd)
      eyeball_prof <- x_t(t = sac_time, xm = sac_amp, 
                          beta = sac_params[fix_nr, 1], mu = sac_params[fix_nr, 2], A = sac_params[fix_nr, 3] )
      eyeball_x <- degree_to_xy(eyeball_prof, sac_dir)[[1]] * pix_per_dva + current_xy[1]
      eyeball_y <- degree_to_xy(eyeball_prof, sac_dir)[[2]] * pix_per_dva + current_xy[2]
      this_saccade <- data.frame(x = sac_x, y = sac_y, time = current_time + sac_time)
      this_saccade$what <- "SAC"
      this_saccade$what[this_saccade$time>=time_sacoff_PSO] <- "PSO"
      x_sacoff <- this_saccade$x[length(this_saccade$x)]
      y_sacoff <- this_saccade$y[length(this_saccade$y)]
      current_time <- this_saccade$time[length(this_saccade$time)] + inter_sample_dur # update to next time step
      current_xy <- as.numeric(this_saccade[length(this_saccade$time), c("x", "y")]) 
    } else { # no corresponding saccade found
      this_saccade <- NULL
      # no saccade, no saccade parameters
      sac_amp <- NaN
      sac_dir <- NaN
      sac_vpeak <- NaN
      sac_vmean <- NaN
      time_sacon <- NaN
      x_sacon <- NaN
      y_sacon <- NaN
      time_sacoff <- NaN
      time_sacoff_PSO <- NaN
      x_sacoff <- NaN
      y_sacoff <- NaN
    }
    # add to final data frame
    resulting_df[[fix_nr]] <- rbind(this_fixation, this_saccade)
    resulting_metrics[[fix_nr]] <- data.frame(t_fixon = time_fixon, # fixation
                                              x_fixon = x_fixon, y_fixon = y_fixon, 
                                              t_fixoff = time_fixoff, 
                                              t_sacon = time_sacon, # saccade 
                                              x_sacon = x_sacon, y_sacon = y_sacon, 
                                              t_sacoff = time_sacoff, t_sacoff_PSO = time_sacoff_PSO,
                                              x_sacoff = x_sacoff, y_sacoff = y_sacoff, 
                                              sac_vpeak = sac_vpeak, sac_vmean = sac_vmean )
  } # end of loop through fixations
  # convert to data frame
  resulting_df <- rbindlist(resulting_df)
  resulting_metrics <- rbindlist(resulting_metrics)
  # want a plot?
  if (do_plot) {
    par(mfrow=c(1,2)) 
    # x 
    plot(resulting_df$time, resulting_df$x, 
         xlab = "Time [ms]", ylab = "X position [pix]")
    abline(v = resulting_metrics$t_sacon, col = "black", lty = 2)
    abline(v = resulting_metrics$t_sacoff_PSO, col = "blue")
    abline(v = resulting_metrics$t_sacoff, col = "black")
    # legend
    legend("bottomright", 
           legend=c("Saccade onset", "PSO onset", "PSO offset"),
           col=c("black", "blue", "black"), 
           lty=c("dashed", "solid", "solid"), cex=0.8)
    # y
    plot(resulting_df$time, resulting_df$y, ylim = rev(range(resulting_df$y)), 
         xlab = "Time [ms]", ylab = "Y position [pix]")
    abline(v = resulting_metrics$t_sacon, col = "black", lty = 2)
    abline(v = resulting_metrics$t_sacoff_PSO, col = "blue")
    abline(v = resulting_metrics$t_sacoff, col = "black")
    par(mfrow=c(1,1)) 
  }
  # return data
  return(list(resulting_df, resulting_metrics))
}


# AUX functions for the simulator
xy_to_degree <- function(vx, vy) {
  # compute direction of movement 360deg. see: http://stackoverflow.com/questions/22421054/determine-movement-vectors-direction-from-velocity
  deg = atan2(vy, vx)*180/pi # larger y mean downward !
  deg[deg < 0] <- deg[deg < 0] + 360
  return(deg)
}

xy_to_dist <- function(vx, vy) {
  dist = sqrt(vx^2 + vy^2)
  return(dist)
}

degree_to_xy <- function(dist, deg) {
  x = dist * cos(deg * pi / 180)
  y = dist * sin(deg * pi / 180)
  return(list(x, y))
}

five_point_smoother <- function(x) {
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

get_vel_vec <- function(pos, t, do_smooth=FALSE) {
  if (do_smooth) {
    pos <- five_point_smoother(pos)
  }
  return(c(0, diff(pos)) / (c(1, diff(t))/1000))
}
