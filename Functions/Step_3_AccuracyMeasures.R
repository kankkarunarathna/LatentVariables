# Function to compute the DRPS (Discrete Ranked Probability Score) and interval coverage
drps_score = function (truth, fc, interval_width = 0.9, log = FALSE)
{
  # Apply log transformation if specified
  if (log) {
    truth <- log(truth + 0.001)  # Avoid log(0) by adding a small constant
    fc <- log(fc + 0.001)
    # Set upper bound for evaluation range based on max of truth and forecast
    nsum <- max(c(truth, fc), na.rm = TRUE) + 5
  } else {
    # Set upper bound using 99th percentile of forecast and a buffer
    nsum <- max(c(truth, quantile(fc, probs = 0.99, na.rm = TRUE)), na.rm = TRUE) + 1000
  }
  
  # Create empirical cumulative distribution function from forecast samples
  Fy = ecdf(fc)
  
  # Generate a sequence of integer values from 0 to nsum
  ysum <- 0:nsum
  
  # Indicator function: 1 if ysum >= truth, else 0
  indicator <- ifelse(ysum - truth >= 0, 1, 0)
  
  # DRPS score: sum of squared differences between indicator and forecast CDF
  score <- sum((indicator - Fy(ysum))^2)
  
  # Compute prediction interval bounds
  interval <- quantile(fc, probs = c((1 - interval_width)/2,
                                     (interval_width + (1 - interval_width)/2)), na.rm = TRUE)
  
  # Check if truth falls within the prediction interval
  in_interval <- ifelse(truth <= interval[2] & truth >= interval[1], 1, 0)
  
  # Return both DRPS score and interval coverage indicator
  return(c(score, in_interval))
}


# Function to apply DRPS scoring across multiple forecasted time points
drps_mcmc_object = function (truth, fc, interval_width = 0.9, log = FALSE)
{
  # Identify indices where truth is not NA
  indices_keep <- which(!is.na(truth))
  
  # If all truth values are NA, return NA-filled data frame
  if (length(indices_keep) == 0) {
    scores = data.frame(drps = rep(NA, length(truth)), interval = rep(NA, length(truth)))
  } else {
    # Initialize matrix to store scores
    scores <- matrix(NA, nrow = length(truth), ncol = 2)
    
    # Loop through valid indices and compute DRPS score and interval coverage
    for (i in indices_keep) {
      scores[i, ] <- drps_score(truth = as.vector(truth)[i],
                                fc = fc[, i], interval_width, log = log)
    }
  }
  
  # Return matrix of scores
  scores
}
