
# Loading packages
library(cmdstanr)     # Interface to CmdStan for Bayesian modeling
library(bayesplot)    # Visualization tools for Bayesian models
library(mvgam)        # Multivariate GAMs for time series modeling
library(ggplot2)      # Data visualization
library(forecast)     # Time series forecasting methods
library(fpp2)         # Additional time series datasets and functions
library(patchwork)    # Combine ggplot2 plots
library(dplyr)        # Data manipulation
library(tidyr)        # Data tidying



# cmdstan_path()
# set_cmdstan_path("C:/Users/s4740365/.cmdstan/cmdstan-2.36.0")

# Source functions
source("Functions/Step_1_DataPreparation.R")  # Data preparation steps
source("Functions/Step_2_Models.R")           # Model fitting functions
source("Functions/Step_3_AccuracyMeasures.R")  # Data preparation steps


# Naming models
mod_names <- c(
  'Latent NDVI based model',
  'Average NDVI based model'
)

mod_fc_results1 <- do.call(
  rbind, 
  lapply(1:60, function(j){
    message(paste0('Fitting CV fold ', j, ' of 60'))
    
    # Define training window
    L <- 1
    U <- 447 + j
    k <- U - (L - 1)
    times_keep <- L:U
    
    # Prepare model data
    model_data <- merged_data %>%
      dplyr::arrange(time) %>%
      slice(times_keep)
    
    model_data$time <- 1:length(times_keep)
    times_keep <- times_keep - (L - 1)
    times_test <- tail(times_keep, 12)
    
    # Extract true values for evaluation
    truth <- model_data[times_test, 4]
    truth[truth == 999] <- NA
    
    # Mask test observations
    model_data[times_test, 4] <- 999
    model_data[times_test, 5:10] <- 999  # NDVI sensors
    
    # Forecast average NDVI using ARIMA
    armod <- auto.arima(ts(model_data[1:(k-12), 11], frequency = 12))
    model_data[(k-11):k, 11] <- as.vector(forecast(armod, h = 12)$mean)
    
    # Fit Stan models
    stan_mods <- modfit1(model_data)
    
    # Evaluate model performance
    mod_performances <- do.call(
      rbind, 
      lapply(seq_along(stan_mods), function(i){
        preds <- stan_mods[[i]]$draws('ypred', format = 'matrix')
        mod_drps <- drps_mcmc_object(truth, preds[, times_test])
        
        data.frame(
          last_train = U - 12,
          model = mod_names[i],
          DRPS = sum(mod_drps[,1], na.rm = TRUE),
          Coverage = sum(mod_drps[,2], na.rm = TRUE) / sum(!is.na(mod_drps[,2]))
        )
      })
    )
    
    # Rank models
    mod_performances$DRPS_rank <- rank(mod_performances$DRPS)
    mod_performances$Coverage_rank <- rank(-mod_performances$Coverage)
    
    mod_performances
  })
)

# Save and reload results
# write.csv(mod_fc_results1, "Results/mod_fc_results1.csv")

mod_fc_results1 <- read.csv("Results/mod_fc_results1.csv")

#####################################################################################################

# Naming models
mod_names <- c("1_source", "2_sources", "3_sources", "4_sources")

# Sample evaluation times when 4 sensors are available
eval_times <- sort(sample(212:358, size = 50, replace = FALSE))

mod_fc_results2 <- do.call(
  rbind,
  lapply(1:50, function(j) {
    message('Modelling CV fold ', j, ' of 50')
    
    times_keep <- 1:eval_times[j]
    model_data <- merged_data %>%
      dplyr::arrange(time) %>%
      slice(times_keep)
    
    model_data$time <- 1:length(times_keep)
    last_train <- eval_times[j] - 12
    times_test <- tail(times_keep, 12)
    
    # Extract true values
    truth <- model_data[times_test, 4]
    truth[truth == 999] <- NA
    
    # Mask test observations
    model_data[times_test, 4] <- 999
    model_data[times_test, 5:10] <- 999
    
    # Fit Stan models
    stan_mods <- modfit2(model_data)
    
    # Evaluate model performance
    mod_performances <- do.call(
      rbind,
      lapply(seq_along(stan_mods), function(i) {
        preds <- stan_mods[[i]]$draws("ypred", format = "matrix")
        mod_drps <- mvgam:::drps_mcmc_object(truth, preds[, times_test])
        
        data.frame(
          model = mod_names[i],
          DRPS = sum(mod_drps[, 1], na.rm = TRUE),
          Coverage = sum(mod_drps[, 2], na.rm = TRUE) / sum(!is.na(mod_drps[, 2])),
          last_train = last_train
        )
      })
    )
    
    # Rank models
    mod_performances$DRPS_rank <- rank(mod_performances$DRPS)
    mod_performances$Coverage_rank <- rank(-mod_performances$Coverage)
    
    mod_performances
  })
)

# Save and reload results
# write.csv(mod_fc_results2, file="Results/mod_fc_results2.csv")
mod_fc_results2 <- read.csv("Results/mod_fc_results2.csv")



##########################################################################################################
 

## Rerun model once to get the coefficients estimates
L = 1   
U = 447 + j  
k <- U - (L - 1)        
times_keep <- L : U
model_data <- data.frame(merged_data %>%
                           dplyr::arrange(time))[times_keep,]

model_data <- model_data %>% dplyr::arrange(time)
model_data$time <- 1:length(times_keep)

times_keep <- times_keep - (L-1)
times_test <- tail(times_keep, 12)

model_data[(times_test), 4] -> truth
truth[truth == 999] <- NA

# Now set last 12 obs of y_agg to 999 in model_data
model_data[times_test, 4] <- 999

# Also set last 12 obs of NDVI sensors to 6 as they will be missing
# in the forecast scenario
model_data[times_test, c(5 : 10)] <- 999

# Fit the three Stan models
stan_mods <- modfit3(model_data)

mod_SS <- stan_mods$mod_SS_P

## To obtain summary of selected  parameters
mod_SS$summary(variables = c(
  "ar1_ndvi",
  "ar2_ndvi",
  "ar12_ndvi",
  "ar24_ndvi",
  "alpha_count",
  "beta_ndvi",
  "beta_ndvisq",
  "ar1_count",
  "ar12_count",
  "sigma_sensor"
))



########################################################################################################

##### Figures for manuscript/chapter

## Figure S1

rep_na = function(x){
  x[x > 998] <- NA  # Replace values >998 with NA (used as missing indicators)
  x
}

# Apply rep_na to all columns in merged_data
merged_da <- merged_data %>%
  dplyr::mutate_all(rep_na)


plot_mvgam_series(data = merged_da, y = 'y')  # Visualize cleaned response series

# Reshape NDVI Data for Faceted Plot
df_long <- ndvid %>%
  pivot_longer(
    cols = c("GIMMSv0", "Landsat5", "Landsat7", "MODIS", "Landsat8", "Landsat9"),
    names_to = "series",          
    values_to = "NDVI"            
  )

df_long$series <- as.factor(df_long$series)
df_long$time <- rep(1:(nrow(df_long)/6), each = 6)
df_long$NDVI[df_long$NDVI > 995] <- NA  # Clean NDVI values


# Faceted NDVI Time Series Plot
ggplot2::ggplot(df_long, aes(x=time, y=NDVI)) + 
  ggplot2::facet_wrap(~series, scales = "free_x") +
  ggplot2::labs(x = "Time", y = 'NDVI')  +
  ggplot2::theme_bw() +
  ggplot2::geom_line(colour = "#8F2727", linewidth = 0.75)


## Figure S2: Plot Average NDVI Series
plot_mvgam_series(data = merged_data, y = 'Avg_NDVI')


## Figure S3: Correlation Among NDVI Sensors

library(corrplot)

# Clean scaled NDVI data
ndvid2 <- ndvi_scaled %>%
  dplyr::mutate_all(rep_na)

# Compute correlation matrices for different sensor groups
M1 <- cor(ndvid2[3:6], use = 'complete.obs')  # MODIS and Landsat sensors
M2 <- cor(ndvid2[7:8], use = 'complete.obs')  # Landsat8 and Landsat9

# Plot correlation matrices side by side
layout(matrix(c(1:2), nrow = 1, ncol = 2))
corrplot.mixed(M1, order = 'AOE')
corrplot.mixed(M2, order = 'AOE')


##  Summarize Posterior Estimates of Selected Parameters
mod_SS <- stan_mods$mod_SS_P
mod_SS$summary(variables = c(
  "ar1_ndvi",
  "ar2_ndvi",
  "ar12_ndvi",
  "ar24_ndvi",
  "alpha_count",
  "beta_ndvi",
  "beta_ndvisq",
  "ar1_count",
  "ar12_count",
  "sigma_sensor"
))


## Figure S4: R-hat Histogram (Convergence Diagnostic)
library(bayesplot)
mcmc_rhat_hist(mod_SS$summary()$rhat)


# Figure S5: Effective Sample Size Histogram
mcmc_neff_hist(mod_SS$summary()$ess_bulk)


# Figure S6: Trace Plots for Key Parameters
mcmc_trace(mod_SS$draws(c(
  "ar1_ndvi",
  "ar2_ndvi",
  "ar12_ndvi",
  "ar24_ndvi",
  "alpha_count",
  "beta_ndvi",
  "beta_ndvisq",
  "ar1_count",
  "ar12_count",
  "sigma_sensor"
)))

# Figure S7: Pairwise Posterior Correlations
mcmc_pairs(mod_SS$draws(c(
  "ar1_ndvi",
  "ar2_ndvi",
  "ar12_ndvi",
  "ar24_ndvi",
  "alpha_count",
  "beta_ndvi",
  "beta_ndvisq",
  "ar1_count",
  "ar12_count"
)))


# Figure S8: Compare True vs. Observed NDVI
true_mean <- as.vector(colMeans(ndvi_true))
Avg_NDVI <- ndvid$Avg_NDVI
ndvi_true_obs <- cbind(time=seq(1:492), true_mean, Avg_NDVI)
ndvi_true_obs[ndvi_true_obs >5] <- NA 
ggplot(ndvi_true_obs, aes(x = Avg_NDVI, y = true_mean)) +
  geom_point() +
  theme_classic() +
  labs(x="Average of observed NDVI", y="Average of estimated true NDVI") 


# Figure 2: Predicted vs. Observed Counts

# Extract posterior draws of predicted counts from the fitted Stan model
ypred_draws <- mod_SS$draws('ypred', format = 'matrix')

# Compute 10th, 50th (median), and 90th percentiles for each time point
lims <- apply(ypred_draws, 2, function(x) quantile(x, probs = c(0.1, 0.5, 0.9)))

# Extract observed counts from merged_data
y1 <- merged_data$y[1:492]

# Clean observed counts: replace values >400 with NA (likely outliers or placeholders)
for (i in 1:492) {
  if (y1[i] > 400) {
    y1[i] <- NA
  } else {
    y1[i] <- y1[i]  # Redundant but keeps structure
  }
}


#####  

# Create a data frame for plotting with time, prediction intervals, and median predictions
ggplot(data.frame(
  time = 1:NCOL(lims),         # Time index (e.g., months)
  ylower = lims[1,],           # 10th percentile of predicted counts
  ymean  = lims[2,],           # Median predicted counts
  yupper = lims[3,]            # 90th percentile of predicted counts
),
aes(x = time, y = ymean)) +    # Set x-axis as time and y-axis as median prediction
  # Add shaded ribbon for 10–90% prediction interval
  geom_ribbon(aes(ymin = ylower, ymax = yupper),
              alpha = 0.4) +   # Transparency for the ribbon
  
  # Add line for median predicted counts
  geom_line() +
  
  # Use a clean, classic theme
  theme_classic() + 
  
  # Add axis labels and title
  labs(title = '') +
  labs(x = "Time (in months)", 
       y = "95% prediction band and observed counts") +
  
  # Overlay observed counts as brown open circles
  geom_point(aes(y = y1), 
             col = 'brown', 
             shape = 1, 
             size = 2) +
  
  # Add a vertical dashed red line at time point 480
  geom_vline(xintercept = 480, 
             linetype = "dashed", 
             color = "red", 
             linewidth = 1.5)



# Figure 3: Predicted True NDVI with Observed NDVI from Multiple Sources

# Define a function to replace placeholder values (>998) with NA
scale_2 = function(x) {
  x[x > 998] <- NA
  x
}

# Combine selected NDVI sensor columns and apply cleaning function
as.data.frame(cbind(
  ndvi_scaled$GIMMSv0,
  ndvi_scaled$Landsat5, 
  ndvi_scaled$Landsat7,
  ndvi_scaled$MODIS, 
  ndvi_scaled$Landsat8, 
  ndvi_scaled$Landsat9,
  ndvi_scaled$Avg_NDVI
)) %>%
  dplyr::mutate_all(scale_2) -> ndvi_scaled2

# Restrict to first 492 time points
ndvi_scaled2 <- ndvi_scaled2[1:492,]

# Extract posterior draws of true NDVI from the fitted Stan model
ndvi_true <- mod_SS$draws('true_ndvi', format = 'matrix')

# Compute 95% credible interval (2.5%, 50%, 97.5%) for each time point
lims <- apply(ndvi_true, 2, function(x) 
  quantile(x, probs = c(0.025, 0.5, 0.975)))

# Compute mean of true NDVI across posterior samples
ndvi_true_mean <- cbind(time = seq(1:492), true_mean = as.vector(colMeans(ndvi_true)))

# Create ggplot with prediction intervals and observed NDVI points
ggplot(data.frame(
  time = 1:NCOL(lims),
  ylower = lims[1,],
  ymean  = lims[2,],
  yupper = lims[3,]
), aes(x = time, y = ymean)) +
  
  # Add shaded ribbon for 95% prediction interval
  geom_ribbon(aes(ymin = ylower, ymax = yupper), 
              fill = "gray", alpha = 0.7) +
  
  # Overlay observed NDVI from each sensor with distinct colors and shapes
  geom_point(aes(y = ndvi_scaled2$V1), colour = 'darkmagenta', shape = 20, size = 1) +  # GIMMSv0
  geom_point(aes(y = ndvi_scaled2$V2), colour = 'green', shape = 3, size = 1) +         # Landsat5
  geom_point(aes(y = ndvi_scaled2$V3), colour = 'orange', shape = 15, size = 1) +       # Landsat7
  geom_point(aes(y = ndvi_scaled2$V4), colour = '#8F2727', shape = 16, size = 1) +      # MODIS
  geom_point(aes(y = ndvi_scaled2$V5), colour = 'blue', shape = 17, size = 1) +         # Landsat8
  geom_point(aes(y = ndvi_scaled2$V6), colour = 'chartreuse1', shape = 18, size = 1) +  # Landsat9
  
  # Add vertical dashed line at time point 480 (492 - 12), marking forecast start
  geom_vline(xintercept = (492 - 12), linetype = "dashed", color = "red", size = 1.5) +
  
  # Apply classic theme for clean presentation
  theme_classic() +
  
  # Add axis labels and title
  labs(title = "",
       x = "Time (in months)",
       y = "95% Prediction Band and Observed NDVI")


# Figure 4: Predicted Average Count vs. NDVI

# Extract posterior draws for NDVI coefficients from the fitted Stan model
beta_ndvi <- mod_SS$draws('beta_ndvi', format = 'matrix')      # Linear term
beta_ndvisq <- mod_SS$draws('beta_ndvisq', format = 'matrix')  # Quadratic term

# Simulate NDVI values across a range for prediction
ndvi_sim <- runif(495, -3, 3)  # 495 simulated NDVI values between -3 and 3

# Generate predicted values using posterior samples for each time point
ndvipreds <- do.call(rbind, lapply(1:492, function(x){
  data.frame(
    ndvi = ndvi_sim,
    pred = beta_ndvi[x] * ndvi_sim + beta_ndvisq[x] * (ndvi_sim)^2  # Quadratic prediction
  )
}))

# Summarize predictions: compute mean and 95% credible intervals for each NDVI value
ndvipreds %>%
  group_by(ndvi) %>%
  dplyr::mutate(
    ylower = quantile(pred, probs = 0.025),  # Lower bound (2.5%)
    ymean  = mean(pred),                     # Mean prediction
    yupper = quantile(pred, probs = 0.975)   # Upper bound (97.5%)
  ) %>%
  dplyr::arrange(ndvi) -> plotdat  # Arrange for smooth plotting


layout(1)  # Set layout for single plot

ggplot(plotdat, aes(x = ndvi, y = ymean)) +
  
  # Add shaded ribbon for 95% prediction interval
  geom_ribbon(aes(ymin = ylower, ymax = yupper),
              fill = "grey70", alpha = 2) +
  
  # Add line for mean predicted values
  geom_line(col = "brown") +
  
  # Apply classic theme for clean presentation
  theme_classic() +
  
  # Remove legend (not needed here)
  theme(legend.position = 'none') +
  
  # Add axis labels and title
  labs(title = "",
       x = "NDVI",
       y = "log(average_predicted true counts)")


# Figure 5: Posterior Distributions of NDVI Autoregressive Coefficients

# Extract posterior draws for autoregressive NDVI terms from the fitted Stan model
ar1_ndvi  <- mod_SS$draws('ar1_ndvi', format = 'matrix')   # Lag-1
ar2_ndvi  <- mod_SS$draws('ar2_ndvi', format = 'matrix')   # Lag-2
ar12_ndvi <- mod_SS$draws('ar12_ndvi', format = 'matrix')  # Lag-12
ar24_ndvi <- mod_SS$draws('ar24_ndvi', format = 'matrix')  # Lag-24

# Set layout for 2x2 grid of histograms
layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE))

# Set plot margins: bottom, left, top, right
par(mar = c(5, 2, 2, 2))

# Plot histogram for AR1 coefficient
hist(ar1_ndvi, 
     main = expression(paste("(a) ", AR1)), 
     border = "brown", 
     col = "grey70", 
     xlab = '', 
     ylab = '')

# Plot histogram for AR2 coefficient
hist(ar2_ndvi, 
     main = expression(paste("(b) ", AR2)), 
     border = "brown", 
     col = "grey70", 
     xlab = '', 
     ylab = '')

# Plot histogram for AR12 coefficient
hist(ar12_ndvi, 
     main = expression(paste("(c) ", AR12)), 
     border = "brown", 
     col = "grey70", 
     xlab = '', 
     ylab = '')

# Plot histogram for AR24 coefficient
hist(ar24_ndvi, 
     main = expression(paste("(d) ", AR24)), 
     border = "brown", 
     col = "grey70", 
     xlab = '', 
     ylab = '')

# Reset layout to default (single plot)
layout(1)


# Figure 6: Posterior Distributions of Count Model Parameters

# Extract posterior draws for count model parameters from the fitted Stan model
alpha_count   <- mod_SS$draws('alpha_count', format = 'matrix')     # Intercept term
beta_ndvi     <- mod_SS$draws('beta_ndvi', format = 'matrix')       # NDVI linear effect
beta_ndvisq   <- mod_SS$draws('beta_ndvisq', format = 'matrix')     # NDVI squared effect
ar1_count     <- mod_SS$draws('ar1_count', format = 'matrix')       # Lag-1 autoregressive term
ar12_count    <- mod_SS$draws('ar12_count', format = 'matrix')      # Lag-12 autoregressive term

# Set layout for 2x3 grid of histograms
layout(matrix(c(1:6), nrow = 2, ncol = 3, byrow = TRUE))

# Set plot margins: bottom, left, top, right
par(mar = c(5, 2, 2, 2))

# Plot histogram for intercept (alpha_count)
hist(alpha_count, 
     main = expression(paste("(a) ", beta[0])), 
     border = "brown", 
     col = "grey70", 
     xlab = '', 
     ylab = '')

# Plot histogram for NDVI linear effect
hist(beta_ndvi, 
     main = expression(paste("(b) ", beta[1])), 
     border = "brown", 
     col = "grey70", 
     xlab = '', 
     ylab = '')

# Plot histogram for NDVI squared effect
hist(beta_ndvisq, 
     main = expression(paste("(c) ", beta[2])), 
     border = "brown", 
     col = "grey70", 
     xlab = '', 
     ylab = '')

# Plot histogram for AR1 term in count model
hist(ar1_count, 
     main = expression(paste("(d) ", AR1)), 
     border = "brown", 
     col = "grey70", 
     xlab = '', 
     ylab = '')

# Plot histogram for AR12 term in count model
hist(ar12_count, 
     main = expression(paste("(e) ", AR12)), 
     border = "brown", 
     col = "grey70", 
     xlab = '', 
     ylab = '')

# Reset layout to default (single plot)
layout(1)



# Figure 7(a): DRPS and Coverage Over Time for NDVI-Based

# Load cross-validation results from CSV file
results <- read.csv(file = "mod_fc_results1.csv")

# Create first plot: DRPS vs. time for each model
p1 <- ggplot(results,
             aes(x = last_train, y = DRPS, color = model, shape = model)) +
  geom_point(size = 1.5, show.legend = FALSE) +  # Plot points without legend
  theme_classic() +                              # Apply clean theme
  labs(title = '(a)',                            # Panel label
       x = "Time (in months)", 
       y = "DRPS")                               # DRPS = Dawid–Sebastiani Score

# Create second plot: Coverage vs. time for each model
p2 <- ggplot(results,
             aes(x = last_train, y = Coverage, color = model, shape = model)) +
  geom_point(size = 1.5) +                       # Plot points with legend
  theme_classic() +                              # Apply clean theme
  labs(title = '(b)',                            # Panel label
       x = "Time (in months)", 
       y = "Coverage")                           # Coverage = proportion of true values within prediction interval

# Combine the two plots side by side using patchwork
p1 + p2


# Figure 8: Distribution of DRPS Ranks and Coverage Across Models

# Plot the performance ranking distributions for the various models
mod_fc_results1 <- read.csv("mod_fc_results1.csv")
p1 <- ggplot(
  mod_fc_results1,
  aes(
    x = model,
    y = DRPS_rank,
    fill = model
  )
) +
  geom_violin() +
  geom_jitter(
    height = 0,
    width = 0.15
  ) +
  scale_fill_viridis_d(alpha = 0.6) +
  theme_classic() +
  labs(
    x = "",
    y = "DRPS rank"
  ) +
  theme(legend.position = "None")


p2 <- ggplot(
  mod_fc_results1,
  aes(
    x = model,
    y = Coverage,
    fill = model
  )
) +
  geom_hline(
    yintercept = 0.9,
    linetype = "dashed"
  ) +
  geom_violin() +
  geom_jitter(
    height = 0,
    width = 0.15
  ) +
  scale_fill_viridis_d(alpha = 0.6) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(
    x = "",
    y = "90% coverage"
  ) +
  theme(legend.position = "None")

patchwork::wrap_plots(p1, p2, ncol = 1)



# Figure 9
# Plot the performance ranking distributions for the various models

# Load model performance results from CSV
mod_fc_results1 <- read.csv("mod_fc_results1.csv")

# --- Panel (a): DRPS Rank Distribution ---

# Create violin plot for DRPS ranks across models
p1 <- ggplot(mod_fc_results1,
             aes(x = model, y = DRPS_rank, fill = model)) +
  geom_violin() +                            # Show distribution shape
  geom_jitter(height = 0, width = 0.15) +    # Add individual data points
  scale_fill_viridis_d(alpha = 0.6) +        # Use Viridis color scale with transparency
  theme_classic() +                          # Apply clean theme
  labs(x = "", y = "DRPS rank") +            # Axis labels
  theme(legend.position = "None")            # Hide legend

# --- Panel (b): Coverage Distribution ---

# Create violin plot for coverage across models
p2 <- ggplot(mod_fc_results1,
             aes(x = model, y = Coverage, fill = model)) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +  # Reference line at 90% coverage
  geom_violin() +                                      # Show distribution shape
  geom_jitter(height = 0, width = 0.15) +              # Add individual data points
  scale_fill_viridis_d(alpha = 0.6) +                  # Use Viridis color scale
  scale_y_continuous(limits = c(0, 1)) +               # Set y-axis limits
  theme_classic() +                                    # Apply clean theme
  labs(x = "", y = "90% coverage") +                   # Axis labels
  theme(legend.position = "None")                      # Hide legend



# Load model performance results from CSV
mod_fc_results1 <- read.csv("mod_fc_results1.csv")

# Filter results to focus on evaluation period 456–475 (rows 41 to 80)
mod_fc_results1 <- mod_fc_results1[41:80,]

# --- Panel (b): DRPS Rank Distribution for Period 456–475 ---
p3 <- ggplot(mod_fc_results1,
             aes(x = model, y = DRPS_rank, fill = model)) +
  geom_violin() +                            # Show distribution shape
  geom_jitter(height = 0, width = 0.15) +    # Add individual data points
  scale_fill_viridis_d(alpha = 0.6) +        # Use Viridis color scale
  theme_classic() +                          # Apply clean theme
  labs(x = "", y = "DRPS rank") +            # Axis labels
  theme(legend.position = "None") +          # Hide legend
  ggtitle("(b) Period: 456:475")             # Panel title

# --- Panel (d): Coverage Distribution for Period 456–475 ---
p4 <- ggplot(mod_fc_results1,
             aes(x = model, y = Coverage, fill = model)) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +  # Reference line at 90% coverage
  geom_violin() +                                      # Show distribution shape
  geom_jitter(height = 0, width = 0.15) +              # Add individual data points
  scale_fill_viridis_d(alpha = 0.6) +                  # Use Viridis color scale
  scale_y_continuous(limits = c(0, 1)) +               # Set y-axis limits
  theme_classic() +                                    # Apply clean theme
  labs(x = "", y = "90% coverage") +                   # Axis labels
  theme(legend.position = "None")                      # Hide legend

# --- Combine all panels using patchwork ---
patchwork::wrap_plots(p1, p3, p2, p4, ncol = 2)


# Figure 9: DRPS Rank and Coverage Across Models with Varying Source Counts


# Load model performance results from CSV
mod_fc_results2 <- read.csv("mod_fc_results2.csv")

# --- Panel (a): DRPS Rank Distribution Across Models ---

p1 <- ggplot(mod_fc_results2,
             aes(x = model, y = DRPS_rank, fill = model)) +
  geom_violin() +                            # Show distribution shape of DRPS ranks
  geom_jitter(height = 0, width = 0.15) +    # Add individual data points
  scale_fill_viridis_d(alpha = 0.6) +        # Use Viridis color scale with transparency
  theme_classic() +                          # Apply clean theme
  labs(x = "", y = "DRPS rank") +            # Axis labels
  theme(legend.position = "None")            # Hide legend

# --- Panel (b): Coverage Distribution Across Models ---

p2 <- ggplot(mod_fc_results2,
             aes(x = model, y = Coverage, fill = model)) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +  # Reference line at 90% coverage
  geom_violin() +                                      # Show distribution shape
  geom_jitter(height = 0, width = 0.15) +              # Add individual data points
  scale_fill_viridis_d(alpha = 0.6) +                  # Use Viridis color scale
  scale_y_continuous(limits = c(0, 1)) +               # Set y-axis limits
  theme_classic() +                                    # Apply clean theme
  labs(x = "", y = "90% coverage") +                   # Axis labels
  theme(legend.position = "None")                      # Hide legend

# --- Combine both plots vertically using patchwork ---
patchwork::wrap_plots(p1, p2, ncol = 1)



###  Statistical Testing: Effect of Model Type on DRPS and Coverage

# Fit a linear model to test whether DRPS varies by model type
ft1 <- lm(DRPS ~ model, data = mod_fc_results2)

# Perform ANOVA to assess significance of model effect on DRPS
anova(ft1)

# Fit a linear model to test whether Coverage varies by model type
ft2 <- lm(Coverage ~ model, data = mod_fc_results2)

# Perform ANOVA to assess significance of model effect on Coverage
anova(ft2)




