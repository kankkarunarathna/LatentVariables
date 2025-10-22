# Load required libraries for data manipulation, visualization, and date handling
library(dplyr)
library(ggplot2)
library(lubridate)
library(portalr)
library(tidyr)


### Note 1: The downloaded and arranged data using the following set of codes has been saved in
###  abundance_data.csv


# Download the latest rodent observation data from the Portal Project
# download_observations()

# Summarize rodent abundance data at the treatment level for long-term plots
# rodent_data1 <- summarize_rodent_data(
#   level = "Treatment",
#   plots = "Longterm",
#   shape = "long",
#   time = "all",
#   output = "abundance",
#   na_drop = FALSE,
#   zero_drop = FALSE,
#   min_traps = 1,
#   min_plots = 1,
#   effort = TRUE,
#   download_if_missing = TRUE,
#   quiet = FALSE
# )

# Filter and transform the abundance data for species "PP" under control treatment

# abundance_data <- rodent_data1 |>
#   filter(treatment == "control" | is.na(treatment),  # Include control plots and missing treatment info
#          species == "PP",                            # Focus on species "PP"
#          newmoonnumber > 68,                         # Filter out early observations
#          newmoonnumber < 576) |>
#   mutate(abundance = as.integer(round(abundance * 4 / nplots, 0))) |>  # Normalize abundance by effort
#   mutate(year = lubridate::year(censusdate),        # Extract year from census date
#          month = lubridate::month(censusdate),      # Extract month
#          time = row_number()) |>                     # Create a time index
#   select(newmoonnumber, censusdate, year, month, time, abundance) |>  # Select relevant columns
#   rename(y = abundance)                              # Rename abundance column to 'y'

# Save the cleaned abundance data to a CSV file
  # write.csv(abundance_data, "Data/abundance_data.csv", row.names = FALSE)

# Read the saved abundance data back into R
abundance_data <- read.csv("Data/abundance_data.csv")

# Load NDVI data and filter for relevant years
ndvid <- read.csv(file = "Data/NDVI_arranged_3.csv") %>%
  dplyr::filter(year > 1982) %>%
  dplyr::filter(year < 2024)


# Define a scaling function that handles missing values coded as >998
scale_no_na = function(x){
  x[x > 998] <- NA                      # Replace placeholder values with NA
  x_sc <- as.vector(scale(x))          # Standardize the data
  x_sc[is.na(x_sc)] <- 999             # Recode NA back to placeholder
  x_sc
}

# Apply scaling function to all NDVI columns
ndvi_scaled <- as.data.frame(cbind(ndvid$GIMMSv0,
                                   ndvid$Landsat5,
                                   ndvid$Landsat7,
                                   ndvid$MODIS,
                                   ndvid$Landsat8,
                                   ndvid$Landsat9,
                                   ndvid$Avg_NDVI)) %>%
  dplyr::mutate_all(scale_no_na)

# Assign appropriate column names to the scaled NDVI data
colnames(ndvi_scaled) <- colnames(ndvid[, 3:9])

# Combine year and month columns with scaled NDVI data
ndvi_scaled <- cbind(ndvid[, 1:2], ndvi_scaled)

# Replace missing abundance values with placeholder (999)
abundance_data$y[is.na(abundance_data$y)] <- 999
# Merge abundance and NDVI datasets by year and month
merged_data <- merge(abundance_data,
                     ndvi_scaled,
                     by = c("year", "month"))

# Arrange and select final columns for modeling or analysis
merged_data %>%
  dplyr::arrange(year, month) %>%
  select(year, month, time, y, GIMMSv0, Landsat5, Landsat7, MODIS, Landsat8, Landsat9, Avg_NDVI) -> merged_data


