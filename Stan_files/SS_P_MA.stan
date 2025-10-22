
data {
  int<lower=0> t; // total number of timepoints
  int<lower=0> n_sensors; // number of NDVI sensors
  array[t] int<lower=0> obs_counts; // observed species counts
  array[t, n_sensors] real obs_ndvi; // ndvi observations
}

transformed data {
}

parameters {
  // ndvi process model parameters
  vector[t] true_ndvi;
  real<lower=-1, upper=1> ar1_ndvi;
  real<lower=-1, upper=1> ar2_ndvi;
  real<lower=-1, upper=1> ar12_ndvi;
  real<lower=-1, upper=1> ar24_ndvi;

  // ndvi observation parameters
  // (if sigmas get too small or too large, the ndvi estimates are nonidentified)
  vector<lower=0, upper=1>[n_sensors] sigma_sensor;

  // count process model parameters
  vector[t] mu;
  real<lower=-1, upper=1> ar1_count;
  real<lower=-1, upper=1> ar12_count;
  real<lower=0> sigma_count;
  real alpha_count;
  real beta_ndvi;
  real beta_ndvisq;
}

transformed parameters {
  vector[t] MA_ndvi;
  MA_ndvi[1 : 12] = true_ndvi[1 : 12];
  for (i in 13 : t) 
    { MA_ndvi[i] = mean(true_ndvi[(i-11):i]);
    }
}

model {
  // ndvi process model priors
  ar1_ndvi ~ std_normal();
  ar2_ndvi ~ std_normal();
  ar12_ndvi ~ std_normal();
  ar24_ndvi ~ std_normal();

  // ndvi observation priors
  sigma_sensor ~ beta(5, 15);
  beta_ndvi ~ std_normal();
  beta_ndvisq ~ std_normal();

  // count process model priors
  ar1_count ~ std_normal();
  ar12_count ~ std_normal();
  alpha_count ~ std_normal();
  sigma_count ~ inv_gamma(1.418, 0.452);

  // latent ndvi process (zero-centred for identifiability)
  true_ndvi[1] ~ normal(0, 0.01); // 0.01
  true_ndvi[2] ~ normal(ar1_ndvi * true_ndvi[1], 1.0); // 1
  for (i in 3 : 12) {
    true_ndvi[i] ~ normal(ar1_ndvi * true_ndvi[i - 1]
                          + ar2_ndvi * true_ndvi[i - 2], 1.0);
  }
  
  
  for (i in 13 : 24) {
      true_ndvi[i] ~ normal(ar1_ndvi * true_ndvi[i - 1]
                            + ar2_ndvi * true_ndvi[i - 2]
                            + ar12_ndvi * true_ndvi[i - 12], 1.0);
  }
  
  for (i in 25 : t) {
      true_ndvi[i] ~ normal(ar1_ndvi * true_ndvi[i - 1]
                            + ar2_ndvi * true_ndvi[i - 2]
                            + ar12_ndvi * true_ndvi[i - 12]
                            + ar24_ndvi * true_ndvi[i - 24], 1.0);
  }

  // latent rodent abundance process
  mu[1] ~ normal(alpha_count + beta_ndvi * MA_ndvi[1]
                 + beta_ndvisq * square(MA_ndvi[1]), sigma_count);

  for (i in 2 : 12) {
      mu[i] ~ normal(alpha_count + beta_ndvi * MA_ndvi[i]
                     + beta_ndvisq * square(MA_ndvi[i])
                     + ar1_count
                       * (mu[i - 1]
                          - (alpha_count + beta_ndvi * MA_ndvi[i - 1])
                          + beta_ndvisq * square(MA_ndvi[i - 1])),
                     sigma_count);
  }
  
  for (i in 13 : t) {
      mu[i] ~ normal(alpha_count + beta_ndvi * MA_ndvi[i]
                     + beta_ndvisq * square(MA_ndvi[i])
                     + ar1_count
                       * (mu[i - 1]
                          - (alpha_count + beta_ndvi * MA_ndvi[i - 1]
                             + beta_ndvisq * square(MA_ndvi[i - 1])))
                     + ar12_count
                       * (mu[i - 12]
                          - (alpha_count + beta_ndvi * MA_ndvi[i - 12]
                             + beta_ndvisq * square(MA_ndvi[i - 12]))),
                     sigma_count);
  }
  
  // likelihood functions
  {
    // observed ndvi likelihood contributions
    for (s in 1 : n_sensors) {
      for (i in 1 : t) {
        if (obs_ndvi[i, s] <= 998) {
          obs_ndvi[i, s] ~ normal(true_ndvi[i], sigma_sensor[s]);
        }
      }
    }

    // observed species count likelihood contributions
    for (i in 1 : t) {
      if (obs_counts[i] <= 998) {
        obs_counts[i] ~ poisson_log(mu[i]);
      }
    }
  }
}

generated quantities {
  // posterior predictions
  array[t] int ypred;
  ypred[1 : t] = poisson_log_rng(mu[1 : t]);
}




