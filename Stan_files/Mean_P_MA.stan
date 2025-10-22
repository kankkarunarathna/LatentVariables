
data {
  int<lower=0> t; // total number of timepoints
  array[t]  int<lower=0> obs_counts; // observed species counts
  vector[t]  mean_ndvi; // ndvi observations
}

transformed data {

}
parameters {
  real alpha_count;
  real<lower=-1, upper=1> ar1_count;
  real<lower=-1, upper=1> ar12_count;
  real beta_ndvi;
  real beta_ndvisq;   
  vector[t] mu;
  real<lower=0> sigma_count;
}
transformed parameters {
  vector[t] MA_ndvi;
  MA_ndvi[1 : 12] = mean_ndvi[1 : 12];
  for (i in 13 : t) 
    { MA_ndvi[i] = mean(mean_ndvi[(i-11):i]);
    }
}
model {
  // count process model priors
  alpha_count ~ std_normal();
  beta_ndvi ~ std_normal();
  beta_ndvisq ~ std_normal();
  ar1_count ~ std_normal(); 
  ar12_count ~ std_normal();  
  

  
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
    
  {
    // observed species count likelihood contributions
    for (i in 1: t){
         if (obs_counts[i] <= 998) {
      obs_counts[i] ~ poisson_log(mu[i]); // /[i]);//poisson_log(mu[i]);
    }
    }
  }
}
generated quantities {
  // posterior predictions
  // array[t] int ypred;
  // ypred[1 : t] = poisson_log_rng(mu[1 : t]);
  array[t] int ypred;
  ypred[1 : t] = poisson_log_rng(mu[1 : t]);
}


