# extracts features from eye movement data using an integration window moving through samples over time.
# This approach was first proposed by Raimondas Zemblys. 
# for more info, see: Zemblys, Niehorster, Komogortsev, and Holmqvist (2018) in Behavior research methods.
# implemented by Richard Schweitzer
# NOTE that the IM2C metric is not used at the moment, but you may uncomment it anytime.
extract_features <- function(x, y, time, pix_per_dva, 
                             window_sizes_ms=c(12, 25, 50, 100), 
                             p_bcea=0.68, use_sgolay=FALSE,
                             #i2mc_n_cluster=2, i2mc_noise_sd=0, 
                             do_scale=FALSE ) {
  require(circular) # for the rayleigh test
  require(assertthat)
  require(data.table)
  assert_that(length(x)==length(y) & length(x)==length(time), 
              msg = c("x, y, time must have the same length!"))
  # loop through time scales
  all_extracted_features <- NULL
  for (window_size_ms in window_sizes_ms) {
    if (use_sgolay) {
      require(signal)
      # produce velocity and acceleration vectors, filtered by a Savitzky–Golay filter with polynomial order 2 and a window size of 12 ms
      x_velvec <- get_vel_vec(pos = x / pix_per_dva, t = time, do_smooth = FALSE) 
      y_velvec <- get_vel_vec(pos = y / pix_per_dva, t = time, do_smooth = FALSE)
      x_velvec <- sgolayfilt(x = x_velvec, p = 2, n = 2*round((window_size_ms/2)/2)+1) # small window
      y_velvec <- sgolayfilt(x = y_velvec, p = 2, n = 2*round((window_size_ms/2)/2)+1)
      x_accvec <- get_vel_vec(x_velvec, time, do_smooth = FALSE)
      y_accvec <- get_vel_vec(y_velvec, time, do_smooth = FALSE)
      x_accvec <- sgolayfilt(x = x_accvec, p = 2, n = 2*round((window_size_ms/2)/2)+1)
      y_accvec <- sgolayfilt(x = y_accvec, p = 2, n = 2*round((window_size_ms/2)/2)+1)
    }
    # run through samples
    extracted_features <- vector(mode = "list", length = length(time))
    for (sample_nr in 1:length(time)) {
      # get window for sample
      window_index <- which(time>=(time[sample_nr]-window_size_ms/2) & time<=(time[sample_nr]+window_size_ms/2))
      small_window_index <- which(time>=(time[sample_nr]-window_size_ms/4) & time<=(time[sample_nr]+window_size_ms/4))
      pre_window_index <- which(time>(time[sample_nr]-window_size_ms) & time<=time[sample_nr])
      post_window_index <- which(time>=time[sample_nr] & time<(time[sample_nr]+window_size_ms))
      vel_acc_window_index <- which(time>=(time[sample_nr]-window_size_ms/8) & time<=(time[sample_nr]+window_size_ms/8))
      # feature 1: med-diff and mean-diff
      x_meddiff <- abs(median(x[post_window_index])-median(x[pre_window_index])) / pix_per_dva
      y_meddiff <- abs(median(y[post_window_index])-median(y[pre_window_index])) / pix_per_dva
      x_meandiff <- abs(mean(x[post_window_index])-mean(x[pre_window_index])) / pix_per_dva
      y_meandiff <- abs(mean(y[post_window_index])-mean(y[pre_window_index])) / pix_per_dva
      # feature 2: std and rms of sample-to-sample displacement, dispersion, and bivariate contour ellipse area
      x_std <- sd(x[window_index]) / pix_per_dva
      y_std <- sd(y[window_index]) / pix_per_dva
      x_rms <- sqrt(mean(diff(x[window_index])^2)) / pix_per_dva
      y_rms <- sqrt(mean(diff(y[window_index])^2)) / pix_per_dva
      dispersion <- (max(x[window_index])-min(x[window_index])) + (max(y[window_index])-min(y[window_index])) / pix_per_dva
      bcea <- sqrt(2*sd(x[window_index]/pix_per_dva)*sd(y[window_index]/pix_per_dva)*(1-p_bcea^2)^(1/2))
      # feature 3: Rayleightest for non-uniformity of directions, TEST: circular variance
      deg <- atan2(diff(y[small_window_index]), diff(x[small_window_index]))
      deg <- circular(x = deg, type = "directions", units = "radians")
      ray_test <- rayleigh.test(deg)
      ray_test_stat <- ray_test$statistic # large values shows non-uniformity
      if (is.null(ray_test_stat)) { ray_test_stat <- NaN }
      dirvar <- var.circular(deg) #range.circular(deg)
      # feature 4: velocity and acceleration (previously filtered)
      if (use_sgolay) {
        x_vel <- abs(x_velvec[sample_nr])
        y_vel <- abs(y_velvec[sample_nr])
        x_acc <- abs(x_accvec[sample_nr])
        y_acc <- abs(y_accvec[sample_nr])
      } else {
        x_velvec <- get_vel_vec(pos = x[vel_acc_window_index] / pix_per_dva, 
                                t = time[vel_acc_window_index], do_smooth = FALSE) 
        y_velvec <- get_vel_vec(pos = y[vel_acc_window_index] / pix_per_dva, 
                                t = time[vel_acc_window_index], do_smooth = FALSE)
        x_vel <- mean(abs(x_velvec))
        y_vel <- mean(abs(y_velvec))
        x_accvec <- get_vel_vec(x_velvec, time[vel_acc_window_index], do_smooth = FALSE)
        y_accvec <- get_vel_vec(y_velvec, time[vel_acc_window_index], do_smooth = FALSE)
        x_acc <- mean(abs(x_accvec))
        y_acc <- mean(abs(y_accvec))
      }
      # feature 5: rms-diff, std-diff, bcea-diff
      x_rms_diff <- abs(sqrt(mean(diff(x[post_window_index])^2)) - sqrt(mean(diff(x[pre_window_index])^2))) / pix_per_dva
      y_rms_diff <- abs(sqrt(mean(diff(x[post_window_index])^2)) - sqrt(mean(diff(y[pre_window_index])^2))) / pix_per_dva
      x_std_diff <- abs(sd(x[post_window_index]) - sd(x[pre_window_index])) / pix_per_dva
      y_std_diff <- abs(sd(y[post_window_index]) - sd(y[pre_window_index])) / pix_per_dva
      bcea_diff <- abs(sqrt(2*sd(x[post_window_index]/pix_per_dva)*sd(y[post_window_index]/pix_per_dva)*(1-p_bcea^2)^(1/2)) -
                         sqrt(2*sd(x[pre_window_index]/pix_per_dva)*sd(y[pre_window_index]/pix_per_dva)*(1-p_bcea^2)^(1/2)))
#      # feature 6: I2MC
#      clust_weight <- i2mc(x = x[small_window_index], y = y[small_window_index],
#                           n_clusters = i2mc_n_cluster, noise_sd = i2mc_noise_sd)
#      # potential feature 7: dirmean-diff and dirvar-diff
#      deg_post <- atan2(diff(y[post_window_index]), diff(x[post_window_index]))
#      deg_post <- circular(x = deg_post, type = "directions", units = "radians")
#      deg_pre <- atan2(diff(y[pre_window_index]), diff(x[pre_window_index]))
#      deg_pre <- circular(x = deg_pre, type = "directions", units = "radians")
#      suppressWarnings({dirmean_diff <- abs(mean.circular(deg_post) - mean.circular(deg_pre))})
#      suppressWarnings({dirvar_diff <- abs(var.circular(deg_post) - var.circular(deg_pre))})
#      if (is.na(dirmean_diff) | is.null(dirmean_diff)) { dirmean_diff <- NaN }
      # sum up to one data.frame
      extracted_features[[sample_nr]] <- data.frame(meddiff = xy_to_dist(x_meddiff, y_meddiff), 
                                                    meandiff = xy_to_dist(x_meandiff, y_meandiff),
                                                    std = xy_to_dist(x_std, y_std), 
                                                    rms = xy_to_dist(x_rms, y_rms), 
                                                    dispersion = dispersion, 
                                                    bcea = bcea,
                                                    ray_test = ray_test_stat, 
                                                    dirvar = as.numeric(dirvar),
                                                    vel = xy_to_dist(x_vel, y_vel),
                                                    acc = xy_to_dist(x_acc, y_acc),
                                                    rms_diff = xy_to_dist(x_rms_diff, y_rms_diff), 
                                                    std_diff = xy_to_dist(x_std_diff, y_std_diff), 
                                                    bcea_diff = bcea_diff#, 
                                                    #i2mc = clust_weight#, 
                                                    #dirmean_diff = as.numeric(dirmean_diff),
                                                    #dirvar_diff = as.numeric(dirvar_diff)
      )
    }
    extracted_features <- rbindlist(extracted_features)
    colnames(extracted_features) <- paste(colnames(extracted_features), window_size_ms, sep = "_")
    all_extracted_features <- cbind(all_extracted_features, extracted_features)
  }
  if (do_scale) {
    # https://stackoverflow.com/questions/24260891/scale-in-data-table-in-r
    all_extracted_features <- all_extracted_features[, lapply(.SD, function(x) as.vector(scale(x)))]
  }
  return(all_extracted_features)
}


# simple version of the I2MC algorithm of Hessels et al., 2017.
# It runs the kmeans algorithm to assign clusters to samples.
# implemented by Richard Schweitzer
i2mc <- function(x, y, n_clusters=2, noise_sd=0) {
  if (noise_sd!=0) {
    x <- x + rnorm(n = length(x), mean = 0, sd = noise_sd)
    y <- y + rnorm(n = length(y), mean = 0, sd = noise_sd)
  }
  clusters <- kmeans(cbind(x, y), n_clusters)
  which_clusters <- clusters$cluster
  # now detect where cluster transitions take place
  where_transitions <- c(FALSE, diff(which_clusters)!=0)
  # how many between-cluster transitions?
  n_transitions <- sum(as.numeric(where_transitions))
  return(1/n_transitions)
}
