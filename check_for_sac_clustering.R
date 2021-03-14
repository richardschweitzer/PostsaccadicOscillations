# Takes in a EK-algorithm output table (from microsacc.R) and looks for saccade clusters
# by Richard Schweitzer
# NOTE: only for the (most common) case of clusters of two saccades, for now.

check_for_sac_clustering <- function(sac_table, n_cluster_samples=50, do_what=1) {
  # check whether there are saccade offsets and saccade onsets are closer than 50 samples..
  # but if so, do what:
  # 0: nothing
  # 1: remove second saccade of cluster
  # 2: merge the two saccades
  if (is.data.frame(sac_table) || is.matrix(sac_table)) {
    nrows <- nrow(sac_table)
  } else if (is.vector(sac_table)) {
    nrows <- 1
  }
  # do we have more than one saccade present?
  if (!is.null(sac_table) && nrows>1) {
    time_between_sacs <- rep(NaN, nrows)
    time_between_sacs_too_small <- rep(FALSE, nrows)
    for (r in 2:nrows) {
      time_between_sacs[r] <- sac_table[r,1] - sac_table[r-1, 2]
      time_between_sacs_too_small[r] <- time_between_sacs[r] < n_cluster_samples
    }
    # Option 1: simply remove the second saccade of the cluster
    if (do_what==0) {
      sac_table_new <- sac_table
    } else if (do_what==1) { 
      sac_table_new <- sac_table[!time_between_sacs_too_small, , drop = FALSE]
    # Option 2: merge somehow
    } else if (do_what==2) { 
      sac_table_new <- matrix(NaN, 
                              nrow = nrows, # number of okay saccades
                              ncol = 7) # number of columns of the table output of the Engbert-Kliegl algorithm
      for (r in 1:nrows) {
        if (r < nrows && !time_between_sacs_too_small[r] && time_between_sacs_too_small[r+1]) {
          # this is the first case: we have cluster of saccades
          sac_table_new[r, 1] <- as.numeric(sac_table[r, 1])
          sac_table_new[r, 2] <- as.numeric(sac_table[r+1, 2])
          sac_table_new[r, 3] <- max(c(as.numeric(sac_table[r, 3]), 
                                       as.numeric(sac_table[r+1, 3])))
          sac_table_new[r, 4] <- as.numeric(sac_table[r, 4] + sac_table[r+1, 4])
          sac_table_new[r, 5] <- as.numeric(sac_table[r, 5] + sac_table[r+1, 5])
          sac_table_new[r, 6] <- as.numeric(sac_table[r, 6] + sac_table[r+1, 6])
          sac_table_new[r, 7] <- as.numeric(sac_table[r, 7] + sac_table[r+1, 7])
        } else if ((r < nrows  && !time_between_sacs_too_small[r] && !time_between_sacs_too_small[r+1]) ||
                   (r == nrows && !time_between_sacs_too_small[r])) {
          # this is the second case: a saccade has been detected just fine
          sac_table_new[r, 1] <- as.numeric(sac_table[r, 1])
          sac_table_new[r, 2] <- as.numeric(sac_table[r, 2])
          sac_table_new[r, 3] <- as.numeric(sac_table[r, 3])
          sac_table_new[r, 4] <- as.numeric(sac_table[r, 4])
          sac_table_new[r, 5] <- as.numeric(sac_table[r, 5])
          sac_table_new[r, 6] <- as.numeric(sac_table[r, 6])
          sac_table_new[r, 7] <- as.numeric(sac_table[r, 7])
        }
      }
      # remove the NaN rows
      sac_table_new <- sac_table_new[!is.na(sac_table_new[,1]), , drop = FALSE]
    }
  } else {
    sac_table_new <- sac_table
  }
  return(sac_table_new)
}