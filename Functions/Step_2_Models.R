
# Compile Stan model for latent NDVI using multiple sensors
stan_SS_P <- cmdstan_model(
  "Stan_files/SS_P_MA.stan",
  stanc_options = list(
    "canonicalize=deprecations,braces,parentheses"
  )
)

# Compile Stan model using average NDVI across sensors
stan_Mean_P <- cmdstan_model(
  "Stan_files/Mean_P_MA.stan",
  stanc_options = list(
    "canonicalize=deprecations,braces,parentheses"
  )
)

# Function to fit two Stan models: one with multiple NDVI sensors, one with average NDVI
modfit1 <- function(model_data) {
  
  # Prepare data for the multi-sensor latent NDVI model
  st_data1 <- list(
    t = nrow(model_data),  # Number of time points
    n_sensors = 6,         # Number of NDVI sensors
    obs_counts = model_data$y,  # Rodent abundance observations
    obs_ndvi = as.data.frame(cbind(
      model_data$GIMMSv0, 
      model_data$Landsat5, 
      model_data$Landsat7,
      model_data$MODIS,
      model_data$Landsat8, 
      model_data$Landsat9
    ))
  )
  
  # Prepare data for the average NDVI model
  st_data2 <- list(
    t = nrow(model_data),
    n_sensors = 1,
    obs_counts = model_data$y,
    mean_ndvi = as.vector(model_data$Avg_NDVI)
  )
  
  # Fit the multi-sensor latent NDVI model
  message('This is mod_SS_P')
  mod_SS_P <- stan_SS_P$sample(
    data = st_data1,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.85
  )
  
  # Fit the average NDVI model
  message('This is mod_Mean_P')
  mod_Mean_P <- stan_Mean_P$sample(
    data = st_data2,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.85
  )
  
  # Return both fitted models
  models <- list(
    SS_M = mod_SS_P,
    Mean_M = mod_Mean_P
  )
  return(models)
}

# Function to fit models using randomly sampled NDVI sensors
modfit2 <- function(model_data) {
  
  # Identify NDVI sensors with at least 10 valid observations
  cols_sample <- 5:10  # NDVI sensor columns
  cols_with_samples <- unlist(
    lapply(cols_sample, function(x){
      obs <- model_data[, x]
      length(which(obs != 999)) >= 10
    })
  )
  
  # Sample 1 sensor
  st_data1 <- list(
    t = nrow(model_data),
    n_sensors = 1,
    obs_counts = model_data$y,
    obs_ndvi = as.data.frame(
      model_data[, sample(cols_sample, 1, replace = FALSE, prob = cols_with_samples)]
    )
  )
  
  # Sample 2 sensors
  st_data2 <- list(
    t = nrow(model_data),
    n_sensors = 2,
    obs_counts = model_data$y,
    obs_ndvi = as.data.frame(
      model_data[, sample(cols_sample, 2, replace = FALSE, prob = cols_with_samples)]
    )
  )
  
  # Sample 3 sensors
  st_data3 <- list(
    t = nrow(model_data),
    n_sensors = 3,
    obs_counts = model_data$y,
    obs_ndvi = as.data.frame(
      model_data[, sample(cols_sample, 3, replace = FALSE, prob = cols_with_samples)]
    )
  )
  
  # Sample 4 sensors
  st_data4 <- list(
    t = nrow(model_data),
    n_sensors = 4,
    obs_counts = model_data$y,
    obs_ndvi = as.data.frame(
      model_data[, sample(cols_sample, 4, replace = FALSE, prob = cols_with_samples)]
    )
  )
  
  # Fit models for each sensor combination
  mod_s_1 <- stan_SS_P$sample(
    data = st_data1,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.85,
    show_messages = FALSE
  )
  
  mod_s_2 <- stan_SS_P$sample(
    data = st_data2,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.85,
    show_messages = FALSE
  )
  
  mod_s_3 <- stan_SS_P$sample(
    data = st_data3,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.85,
    show_messages = FALSE
  )
  
  mod_s_4 <- stan_SS_P$sample(
    data = st_data4,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.85,
    show_messages = FALSE
  )
  
  # Return all fitted models
  models <- list(
    s_1 = mod_s_1,
    s_2 = mod_s_2,
    s_3 = mod_s_3,
    s_4 = mod_s_4
  )
  
  return(models)
}

# Function to fit a single model using all six NDVI sensors
modfit3 <- function(model_data) {
  
  # Prepare full NDVI sensor data
  st_data <- list(
    t = nrow(model_data),
    n_sensors = 6,
    obs_counts = model_data$y,
    obs_ndvi = as.data.frame(cbind(
      model_data$GIMMSv0, 
      model_data$Landsat5, 
      model_data$Landsat7,
      model_data$MODIS,
      model_data$Landsat8, 
      model_data$Landsat9
    ))
  )
  
  # Fit the model using all sensors
  message('This is mod_SS_P')
  mod_SS_P <- stan_SS_P$sample(
    data = st_data,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 500,
    adapt_delta = 0.85
  )
  
  # Return the fitted model
  models <- list(
    mod_SS_P = mod_SS_P
  )
  return(models)
}
