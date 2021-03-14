# This is the PSO detection method we can run on the saccades detected by the Engbert-Kliegl algorithm
# There are two methods at the moment: direction inversion and minimum velocity. 
# by Richard Schweitzer

# The algorithm adds 4 columns to the saccade table
# [,8] direction-inversion sample index
# [,9] overall saccade direction in degrees
# [,10] minimum-velocity sample index
# [,11] minimum velocity detected

check_for_PSO <- function(eye_x, eye_y, samp_rate, em_table, 
                          direction_range = 50) {
  # by Richard Schweitzer
  n_saccades_in_table <- nrow(em_table) # how many saccades were detected using the Engbert-kliegl algorithm?
  if (!is.null(em_table) && !is.null(n_saccades_in_table) && n_saccades_in_table>=1) {
    # new columns for (1) index direction inversion, (2) saccade direction, 
    # (3) index minimum velocity, (4) minimum velocity
    em_table <- cbind(em_table, 
                      rep(NaN, n_saccades_in_table), rep(NaN, n_saccades_in_table), 
                      rep(NaN, n_saccades_in_table), rep(NaN, n_saccades_in_table) )
    n_cols_in_table <- ncol(em_table)
    # check for every saccade
    for (current_saccade in 1:n_saccades_in_table) {
      # extract the offline detection results
      offline_onset <- em_table[current_saccade,1]
      offline_offset <- em_table[current_saccade,2]
      # what is the overall direction of the saccade?
      overall_sac_direction <- atan2(eye_y[offline_offset]-eye_y[offline_onset], 
                                     eye_x[offline_offset]-eye_x[offline_onset])*180/pi
      if (overall_sac_direction < 0) { overall_sac_direction <- overall_sac_direction + 360 }
      # what is the direction criterion?
      opposite_direction_limits <- c(max = overall_sac_direction+direction_range, # max
                                     min = overall_sac_direction-direction_range) # min
      opposite_direction_limits[opposite_direction_limits>360] <- opposite_direction_limits[opposite_direction_limits>360] - 360
      opposite_direction_limits[opposite_direction_limits<0] <- opposite_direction_limits[opposite_direction_limits<0] + 360
      opposite_direction_limits[opposite_direction_limits==0] <- opposite_direction_limits[opposite_direction_limits==0] + 0.00001
      # where should we start looking for the inversion point?
      start_looking <- floor(offline_onset+(offline_offset-offline_onset)/5)
      # compute direction vector
      sac_directions <- vecdir(x = as.matrix(cbind(eye_x, eye_y)), SAMPLING = samp_rate, TYPE = 2)
      sac_velocities <- vecvel(x = as.matrix(cbind(eye_x, eye_y)), SAMPLING = samp_rate, TYPE = 2)
      sac_velocities <- sqrt(sac_velocities[,1]^2 + sac_velocities[,2]^2)
      # extract the point, where direction deviates from the overall saccade direction, i.e., likely the PSO
      direction_inversion <- rep(NaN, times = length(sac_directions))
      minimum_velocity <- rep(NaN, times = length(sac_directions))
      for (it in start_looking:offline_offset) { # run through samples
        # check directions
        if (opposite_direction_limits[1] > opposite_direction_limits[2]) { # max > min, basically any saccade
          if (sac_directions[it] < opposite_direction_limits[1] && 
              sac_directions[it] > opposite_direction_limits[2]) {
            direction_inversion[it] = 0
          } else {
            direction_inversion[it] = 1
          }
        } else { # min > max, the case for a rightward saccade
          if ( (sac_directions[it] >= 0.0 && sac_directions[it] < opposite_direction_limits[1]) ||
               (sac_directions[it] <= 360.0 && sac_directions[it] > opposite_direction_limits[2]) ) {
            direction_inversion[it] = 0
          } else {
            direction_inversion[it] = 1
          }
        }
        # check velocities, is the current velocity the minimum velocity, but not the last above-threshold sample?
        minimum_velocity[it] <- as.numeric(sac_velocities[it]==min(sac_velocities[start_looking:offline_offset]) & 
                                             it!=offline_offset)
      }
      # DIRECTION: what's the first direction mismatch? (if there's any)
      if (any(!is.nan(direction_inversion) & direction_inversion==1)) {
        first_mismatch <- min(which(!is.nan(direction_inversion) & direction_inversion==1))
        em_table[current_saccade, n_cols_in_table-3] <- first_mismatch
        em_table[current_saccade, n_cols_in_table-2] <- overall_sac_direction
      } 
      # VELOCITY: where's our minimal velocity (if not the last sample)
      if (any(!is.nan(minimum_velocity) & minimum_velocity==1)) {
        first_mismatch <- min(which(!is.nan(minimum_velocity) & minimum_velocity==1))
        em_table[current_saccade, n_cols_in_table-1] <- first_mismatch
        em_table[current_saccade, n_cols_in_table] <- sac_velocities[first_mismatch]
      } 
    } # across saccades
  } else {
    print(paste("A valid em_table must be provided. check_for_PSO received:"))
    print(em_table)
  }
  return(em_table)
}
