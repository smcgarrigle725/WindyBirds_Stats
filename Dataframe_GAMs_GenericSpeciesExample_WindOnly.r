# Generic Species at Generic Spatial Grain Example

##### Libraries Needed
# Set generic working directory
setwd("path/to/your/project/data")

# Load packages
# Data manipulation & tidy tools
library(tidyverse)     # includes dplyr, ggplot2, and others
library(lubridate)
# Spatial & raster data
library(sf)
library(raster)
library(dggridR)
library(rnaturalearth)
# Plotting & graphics
library(ggpubr)
library(gridExtra)
library(viridis)
library(fields)
library(corrplot)
# Statistical modeling
library(mgcv)
library(pscl)
library(MASS)
library(FSA)
library(fitdistrplus)
library(hglm)
library(zigam)
library(grplasso)
# Model comparison & diagnostics
library(AICcmodavg)
library(MuMIn)
library(DescTools)
library(pdp)
# eBird specific
library(auk)
library(ebirdst)

# Resolve namespace conflicts
select <- dplyr::select

##### Part 1: eBird Data Extraction and Preparation
# Setup auk eBird data (update paths as appropriate)
ebd <- auk_ebd("raw_data/ebd_species_raw.txt", file_sampling = "raw_data/ebd_sampling_raw.txt")

# Initial filters (replace species name with your generic species)
ebd_filters <- ebd %>% 
  auk_species("species_name_here") %>% 
  auk_year(year = 2012:2021) %>% 
  auk_date(date = c("*-06-14","*-08-03")) %>% 
  auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
  auk_complete()

# Data directory and output paths (generic)
data_dir <- "data"
if (!dir.exists(data_dir)) dir.create(data_dir)

f_ebd <- file.path(data_dir, "ebd_species_breed.txt")
f_sampling <- file.path(data_dir, "ebd_sampling_species_breed.txt")

# Filter data if not already done
if (!file.exists(f_ebd)) {
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling)
}

# Zero-fill to include nondetections
ebd_zf <- auk_zerofill(f_ebd, f_sampling, collapse = TRUE)

# Convert observation start time to decimal hours
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x)/60 + second(x)/3600
}

# Clean and transform variables
ebd_zf <- ebd_zf %>% 
  mutate(
    observation_count = if_else(observation_count == "X", NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    effort_distance_km = if_else(protocol_type != "Traveling", 0, effort_distance_km),
    time_observations_started = time_to_decimal(time_observations_started),
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# Additional filtering for effort, observers
ebd_filtered <- ebd_zf %>% 
  filter(
    duration_minutes <= 300,
    effort_distance_km <= 5,
    number_observers <= 10
  )

# Select needed columns
ebird_data <- ebd_filtered %>% 
  select(
    checklist_id, observer_id, sampling_event_identifier,
    scientific_name, observation_count, species_observed,
    state_code, locality_id, latitude, longitude,
    protocol_type, all_species_reported,
    observation_date, year, day_of_year,
    time_observations_started,
    duration_minutes, effort_distance_km,
    number_observers
  )
# Spatial & Temporal Subsampling
# Load US border shapefile (generic path)
USborder <- st_read("shapefiles/US_border.shp") %>% st_transform(crs = 4326)

# Load species range map (generic path)
rangemap_allseasons <- st_read("range_maps/species_range.gpkg")

# Filter range for breeding season
rangemap_season <- filter(rangemap_allseasons, season == "breeding")

# Clip range to US border
rangemap <- st_intersection(rangemap_season, USborder)

# Convert checklist data to sf points for spatial intersection
ebird_sf <- st_as_sf(ebird_data, coords = c("longitude", "latitude"), crs = 4326)

# Keep only points inside species breeding range
ebird_filtered <- st_intersection(ebird_sf, rangemap) %>% as_tibble()

# Extract lat/lon from geometry
ebird_filtered <- ebird_filtered %>% 
  mutate(
    longitude = map_dbl(geometry, 1),
    latitude = map_dbl(geometry, 2)
  )

# Remove checklists with missing counts
ebird_filtered <- ebird_filtered %>% filter(!is.na(observation_count))

# Create spatial grid with dggridR (spacing in degrees, adjust as needed)
dggs <- dgconstruct(spacing = 9)  # ~4.5km grid approx

# Assign hex grid cell and week number to each checklist
ebird_filtered <- ebird_filtered %>% 
  mutate(
    cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum,
    week = week(observation_date)
  )

# Calculate mean abundance and effort per year per cell
ebird_summary <- ebird_filtered %>%
  mutate(cell_year = paste(cell, year, sep = "_")) %>%
  group_by(cell_year) %>%
  summarise(
    meanCount = mean(observation_count),
    year = max(year),
    meanStartTime = mean(time_observations_started),
    meanSampDur = mean(duration_minutes),
    meanDIST = mean(effort_distance_km),
    meanNumObs = mean(number_observers)
  ) %>%
  ungroup() %>%
  mutate(species_observed = if_else(meanCount > 0, "TRUE", "FALSE"))

# Subsample effort by year to equalize sampling
min_year_count <- min(table(ebird_summary$year))

ebird_summary_ss <- ebird_summary %>%
  group_by(year) %>%
  sample_n(size = min_year_count) %>%
  ungroup()

# Save cleaned and subsampled data (generic path)
write_csv(ebird_summary_ss, file.path(data_dir, "species_breed_finaldata.csv"))

##### Part 2: Importing and Joining Datasets
# Load raw bird data and select relevant columns
spec_raw_spatialgrain <- read_csv(file.choose())  # species_season_finaldata.csv
spec_clean_spatialgrain <- spec_raw_spatialgrain[, c(1:11, 36:45)]  # keep ID, effort, landcover, and response info

# Load land cover data and joing by cell-year
landcover <- read.csv("C:/Users/eller/Documents/OSU/Analysis/Data/Landcover/4.5kmgrid/landcover4.5_cell_PLAND.csv", header=TRUE)
landcover <- landcover  %>% select(-seqnum,-year)
spec_joined_spatialgrain <- left_join(spec_clean_spatialgrain, landcover, by = 'cell_year')

# Load turbine data and join by cell-year
turbines <- read.csv(file.choose(), header = TRUE)  # WRD_foranalysis.csv
spec_joined_spatialgrain <- left_join(spec_joined_spatialgrain, turbines, by = "cell_year") %>%
  mutate(WindCount = ifelse(is.na(WindCount), 0, WindCount))

# Prevalence after subsampling
spec_joined_spatialgrain_ss %>%
  count(species_observed) %>%
  mutate(percent = n / sum(n))

# Final cleaned dataset to be used in downstream analyses
spec_final <- spec_joined_spatialgrain_ss


##### Part 3: Categorizing Wind Variables
# Presence/absence of turbines
spec_joined_spatialgrain$WindPA <- as.factor(ifelse(spec_joined_spatialgrain$WindCount > 0, "Y", "N"))

# Wind age bins
spec_joined_spatialgrain$WindAge_cat <- ifelse(spec_joined_spatialgrain$WindPA == "N", "None",
                                         ifelse(spec_joined_spatialgrain$WindAge <= 2, "0-2",
                                         ifelse(spec_joined_spatialgrain$WindAge <= 4, "3-4",
                                         ifelse(spec_joined_spatialgrain$WindAge <= 6, "5-6",
                                         ifelse(spec_joined_spatialgrain$WindAge <= 8, "7-8", "9+")))))
spec_joined_spatialgrain$WindAge_cat <- factor(spec_joined_spatialgrain$WindAge_cat,
                                                levels = c("None", "0-2", "3-4", "5-6", "7-8", "9+"))

# Wind height bins
spec_joined_spatialgrain$WindHeight_cat <- ifelse(spec_joined_spatialgrain$WindPA == "N", "None",
                                            ifelse(spec_joined_spatialgrain$WindHeight <= 20, "1-20",
                                            ifelse(spec_joined_spatialgrain$WindHeight <= 40, "21-40",
                                            ifelse(spec_joined_spatialgrain$WindHeight <= 60, "41-60",
                                            ifelse(spec_joined_spatialgrain$WindHeight <= 80, "61-80", "81-100+")))))
spec_joined_spatialgrain$WindHeight_cat <- factor(spec_joined_spatialgrain$WindHeight_cat,
                                                  levels = c("None", "1-20", "21-40", "41-60", "61-80", "81-100+"))

# Rotor swept area bins
spec_joined_spatialgrain$WindRSA_cat <- ifelse(spec_joined_spatialgrain$WindPA == "N", "None",
                                         ifelse(spec_joined_spatialgrain$WindRSA <= 5000, "0-5000",
                                         ifelse(spec_joined_spatialgrain$WindRSA <= 10000, "5001-10000",
                                         ifelse(spec_joined_spatialgrain$WindRSA <= 15000, "10001-15000", "15000+"))))
spec_joined_spatialgrain$WindRSA_cat <- factor(spec_joined_spatialgrain$WindRSA_cat,
                                                levels = c("None", "0-5000", "5001-10000", "10001-15000", "15000+"))

# Wind capacity bins
spec_joined_spatialgrain$WindCap_cat <- ifelse(spec_joined_spatialgrain$WindPA == "N", "None",
                                         ifelse(spec_joined_spatialgrain$WindCap <= 1000, "0-1000",
                                         ifelse(spec_joined_spatialgrain$WindCap <= 2000, "1001-2000",
                                         ifelse(spec_joined_spatialgrain$WindCap <= 3000, "2001-3000",
                                         ifelse(spec_joined_spatialgrain$WindCap <= 4000, "3001-4000", "4001+")))))
spec_joined_spatialgrain$WindCap_cat <- factor(spec_joined_spatialgrain$WindCap_cat,
                                                levels = c("None", "0-1000", "1001-2000", "2001-3000", "3001-4000", "4001+"))
# Drop missing rows for key wind attributes
spec_filtered_spatialgrain <- spec_joined_spatialgrain %>%
  drop_na(WindHeight_cat, WindCap_cat, WindRSA_cat)

##### Part 4: Correlation Matrices
# Select variables to assess multicollinearity
cortest <- spec_filtered_spatialgrain %>%
  select(WindCount,
         Developed, Cropland, GrassShrub, TreeCover, Water, Wetland, Barren, IceSnow,
         meanStartTime, meanSampDur, meanDIST, meanNumObs)

# Correlation matrix
corr_matrix <- cor(cortest, use = "complete.obs")

# Reorder columns for visual consistency
corr_matrix <- corr_matrix[, c(
  "Developed", "Cropland", "GrassShrub", "TreeCover", "Water", "Wetland", "Barren", "IceSnow",
  "meanStartTime", "meanSampDur", "meanDIST", "meanNumObs", "WindCount"
)]

# Save
write.csv(corr_matrix, "spec_corr_matrix.csv", row.names = TRUE)

##### Part 5: Standardization
# Select covariates to scale
scaled_cols <- spec_filtered_spatialgrain %>%
  select(meanStartTime, meanSampDur, meanDIST, meanNumObs,
         Developed, Cropland, GrassShrub, TreeCover, Water, Wetland, Barren, IceSnow,
         WindCount, Lat, Long)

# Save scaling parameters
scaling_params <- scaled_cols %>%
  summarise(across(where(is.numeric), list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE)))) %>%
  pivot_longer(everything(), names_to = c("variable", "statistic"), names_sep = "_", values_to = "value") %>%
  pivot_wider(names_from = "statistic", values_from = "value")

# Apply scaling
scaled_data <- scaled_cols %>% 
  mutate(across(where(is.numeric), scale))

# Combine with categorical and response columns
spec_scaled_spatialgrain <- cbind(
  meanCount = spec_filtered_spatialgrain$meanCount,
  meanCount_round = spec_filtered_spatialgrain$meanCount_round,
  WindAge = spec_filtered_spatialgrain$WindAge_cat,
  WindHeight = spec_filtered_spatialgrain$WindHeight_cat,
  WindCap = spec_filtered_spatialgrain$WindCap_cat,
  WindRSA = spec_filtered_spatialgrain$WindRSA_cat,
  WindPA = spec_filtered_spatialgrain$WindPA,
  scaled_data
)

##### Part 6: Split Dataset into Training and Testing Sets
# Select model inputs
spec_split_spatialgrain <- spec_scaled_spatialgrain %>% 
  select(meanCount, meanCount_round,
         meanStartTime, meanSampDur, meanDIST, meanNumObs,
         Developed, Cropland, GrassShrub, TreeCover, Water, Wetland, Barren, IceSnow,
         WindCount, WindHeight, WindCap, WindAge, WindPA, WindRSA,
         Lat, Long)

# Random 80/20 split
spec_split_spatialgrain <- spec_split_spatialgrain %>%
  split(if_else(runif(nrow(.)) <= 0.8, "train", "test"))

map_int(spec_split_spatialgrain, nrow)


##### Part 7: Evaluating Distribution Options using Null Model
# GAM parameters
k <- 5        # knots for smooths
k_time <- 7   # knots for cyclic time
time_knots <- list(meanStartTime = seq(0, 24, length.out = k_time))

# Null model formula
Null <- meanCount_round ~ 1

# Fit null models with different families
start.time <- Sys.time()
gam_nb <- gam(Null, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time(); print(round(end.time - start.time, 2))

start.time <- Sys.time()
gam_pois <- gam(Null, data = spec_split_spatialgrain$train, family = "poisson", knots = time_knots)
end.time <- Sys.time(); print(round(end.time - start.time, 2))

start.time <- Sys.time()
gam_zinb <- zinbgam(Null, pi.formula = ~1, data = spec_split_spatialgrain$train, knots = time_knots)
end.time <- Sys.time(); print(round(end.time - start.time, 2))

start.time <- Sys.time()
gam_zip <- gam(Null, data = spec_split_spatialgrain$train, family = "ziP", knots = time_knots)
end.time <- Sys.time(); print(round(end.time - start.time, 2))

# Also run a quasi-Poisson
gam_qp <- gam(Null, data = spec_split_spatialgrain$train, family = "quasipoisson", knots = time_knots)

# Compare model fit
AIC(gam_nb, gam_pois, gam_zip)
BIC(gam_nb, gam_pois, gam_zip)

# Manually calculate BIC for zinb
AIC_zinb <- gam_zinb$aic
BIC_zinb <- AIC_zinb + log(nrow(spec_split_spatialgrain$train)) * (length(coef(gam_zinb)) - length(fitted(gam_zinb)) + 1)
print(AIC_zinb)
print(BIC_zinb)

# Quasi-Poisson summary
sum.gam <- summary(gam_qp)
sum.gam$p.pv  # parametric p-values
sum.gam$s.pv  # smooth term p-values

# Residuals
resid <- residuals(gam_qp, type = "response")
sum(resid)


##### Part 8: Base Model Plus Additional Variable (assessing edf scores)
# Model formulas for effort covariates
base <- meanCount ~ s(Lat, k = 5) + s(Long, k = 5)
effort_SampDur <- meanCount ~ s(meanSampDur, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
effort_NumObs <- meanCount ~ s(meanNumObs, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
effort_Dist <- meanCount ~ s(meanDIST, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
effort_StartTime <- meanCount ~ s(meanStartTime, bs = "cc", k = 7) + s(Lat, k = 5) + s(Long, k = 5)
effort_Global <- meanCount ~ s(meanSampDur, k = 5) + s(meanNumObs, k = 5) + s(meanDIST, k = 5) + 
                 s(meanStartTime, bs = "cc", k = 7) + s(Lat, k = 5) + s(Long, k = 5)

# Fit models for effort covariates
start.time <- Sys.time()
m_effort_base <- gam(base, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_effort_base)

start.time <- Sys.time()
m_effort_SampDur <- gam(effort_SampDur, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_effort_SampDur)

start.time <- Sys.time()
m_effort_NumObs <- gam(effort_NumObs, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_effort_NumObs)

start.time <- Sys.time()
m_effort_Dist <- gam(effort_Dist, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_effort_Dist)

start.time <- Sys.time()
m_effort_StartTime <- gam(effort_StartTime, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_effort_StartTime)

start.time <- Sys.time()
m_effort_global <- gam(effort_Global, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_effort_global)

MuMIn::AICc(m_effort_base, m_effort_SampDur, m_effort_NumObs, m_effort_Dist, m_effort_StartTime, m_effort_global)

# Land cover model formulas
lc_developed <- meanCount ~ s(Developed, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
lc_cropland <- meanCount ~ s(Cropland, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
lc_grasshrub <- meanCount ~ s(GrassShrub, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
lc_treecover <- meanCount ~ s(TreeCover, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
lc_water <- meanCount ~ s(Water, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
lc_wetland <- meanCount ~ s(Wetland, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
lc_barren <- meanCount ~ s(Barren, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
lc_icesnow <- meanCount ~ s(IceSnow, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
lc_base <- meanCount ~ s(Lat, k = 5) + s(Long, k = 5)

# Fit land cover models
start.time <- Sys.time()
m_lc_base <- gam(lc_base, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_lc_base)

start.time <- Sys.time()
m_lc_developed <- gam(lc_developed, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_lc_developed)

start.time <- Sys.time()
m_lc_cropland <- gam(lc_cropland, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_lc_cropland)

start.time <- Sys.time()
m_lc_grasshrub <- gam(lc_grasshrub, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_lc_grasshrub)

start.time <- Sys.time()
m_lc_treecover <- gam(lc_treecover, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_lc_treecover)

start.time <- Sys.time()
m_lc_water <- gam(lc_water, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_lc_water)

start.time <- Sys.time()
m_lc_wetland <- gam(lc_wetland, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_lc_wetland)

start.time <- Sys.time()
m_lc_barren <- gam(lc_barren, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_lc_barren)

start.time <- Sys.time()
m_lc_icesnow <- gam(lc_icesnow, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_lc_icesnow)
MuMIn::AICc(m_lc_base, m_lc_developed, m_lc_cropland, m_lc_grasshrub, m_lc_treecover, m_lc_water, m_lc_wetland, m_lc_barren, m_lc_icesnow)

# Wind Models Formula
Wind_windcount_simple <- meanCount ~ s(WindCount, k = 5) + s(Lat, k = 5) + s(Long, k = 5)
Wind_windPA_simple <- meanCount ~ factor(WindPA) + s(Lat, k = 5) + s(Long, k = 5)
Wind_windrsa_simple <- meanCount ~ factor(WindRSA) + s(Lat, k = 5) + s(Long, k = 5)
Wind_windheight_simple <- meanCount ~ factor(WindHeight) + s(Lat, k = 5) + s(Long, k = 5)
Wind_windage_simple <- meanCount ~ factor(WindAge) + s(Lat, k = 5) + s(Long, k = 5)
Wind_windcap_simple <- meanCount ~ factor(WindCap) + s(Lat, k = 5) + s(Long, k = 5)

# Fit Wind Models
start.time <- Sys.time()
m_wind_count_simple <- gam(Wind_windcount_simple, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_wind_count_simple)

start.time <- Sys.time()
m_wind_cap_simple <- gam(Wind_windcap_simple, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_wind_cap_simple)

start.time <- Sys.time()
m_Wind_windPA_simple <- gam(Wind_windPA_simple, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_Wind_windPA_simple)

start.time <- Sys.time()
m_wind_rsa_simple <- gam(Wind_windrsa_simple, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_wind_rsa_simple)

start.time <- Sys.time()
m_wind_height_simple <- gam(Wind_windheight_simple, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_wind_height_simple)

start.time <- Sys.time()
m_wind_age_simple <- gam(Wind_windage_simple, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 2)
time.taken
summary(m_wind_age_simple)

MuMIn::AICc(m_wind_count_simple, m_wind_cap_simple, m_Wind_windPA_simple, m_wind_rsa_simple, m_wind_height_simple, m_wind_age_simple, m_lc_base)


##### Part 9: Wind Only Models with Effort and Land Cover Variables
# Model formulas
Wind_windcount <- meanCount ~ 
  s(meanSampDur, k = 5) + s(meanNumObs, k = 5) + s(meanDIST, k = 5) + 
  s(meanStartTime, bs = "cc", k = 7) + 
  s(Developed, k = 5) + s(TreeCover, k = 5) + s(Water, k = 5) + 
  Wetland + s(Barren, k = 5) + s(IceSnow, k = 5) + 
  s(Lat, k = 5) + s(Long, k = 5) + 
  s(WindCount, k = 5) 
Wind_windcap <- meanCount ~ 
  s(meanSampDur, k = 5) + s(meanNumObs, k = 5) + s(meanDIST, k = 5) + 
  s(meanStartTime, bs = "cc", k = 7) + 
  s(Developed, k = 5) + s(TreeCover, k = 5) + s(Water, k = 5) + 
  Wetland + s(Barren, k = 5) + s(IceSnow, k = 5) + 
  s(Lat, k = 5) + s(Long, k = 5) + 
  factor(WindCap) 
Wind_windPA <- meanCount ~ 
  s(meanSampDur, k = 5) + s(meanNumObs, k = 5) + s(meanDIST, k = 5) + 
  s(meanStartTime, bs = "cc", k = 7) + 
  s(Developed, k = 5) + s(TreeCover, k = 5) + s(Water, k = 5) + 
  Wetland + s(Barren, k = 5) + s(IceSnow, k = 5) + 
  s(Lat, k = 5) + s(Long, k = 5) + 
  factor(WindPA)  
Wind_windrsa <- meanCount ~ 
  s(meanSampDur, k = 5) + s(meanNumObs, k = 5) + s(meanDIST, k = 5) + 
  s(meanStartTime, bs = "cc", k = 7) + 
  s(Developed, k = 5) + s(TreeCover, k = 5) + s(Water, k = 5) + 
  Wetland + s(Barren, k = 5) + s(IceSnow, k = 5) + 
  s(Lat, k = 5) + s(Long, k = 5) + 
  factor(WindRSA)  
Wind_windheight <- meanCount ~ 
  s(meanSampDur, k = 5) + s(meanNumObs, k = 5) + s(meanDIST, k = 5) + 
  s(meanStartTime, bs = "cc", k = 7) + 
  s(Developed, k = 5) + s(TreeCover, k = 5) + s(Water, k = 5) + 
  Wetland + s(Barren, k = 5) + s(IceSnow, k = 5) + 
  s(Lat, k = 5) + s(Long, k = 5) + 
  factor(WindHeight)  
Wind_windage <- meanCount ~ 
  s(meanSampDur, k = 5) + s(meanNumObs, k = 5) + s(meanDIST, k = 5) + 
  s(meanStartTime, bs = "cc", k = 7) + 
  s(Developed, k = 5) + s(TreeCover, k = 5) + s(Water, k = 5) + 
  Wetland + s(Barren, k = 5) + s(IceSnow, k = 5) + 
  s(Lat, k = 5) + s(Long, k = 5) + 
  factor(WindAge)  
Wind_base <- meanCount ~ 
  s(meanSampDur, k = 5) + s(meanNumObs, k = 5) + s(meanDIST, k = 5) + 
  s(meanStartTime, bs = "cc", k = 7) + 
  s(Developed, k = 5) + s(TreeCover, k = 5) + s(Water, k = 5) + 
  Wetland + s(Barren, k = 5) + s(IceSnow, k = 5) + 
  s(Lat, k = 5) + s(Long, k = 5)

# Fit GAMs and time runs
start.time <- Sys.time()
m_wind_count <- gam(Wind_windcount, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
cat("Time for m_wind_count: ", round(end.time - start.time, 2), "seconds\n")
summary(m_wind_count)

start.time <- Sys.time()
m_wind_cap <- gam(Wind_windcap, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
cat("Time for m_wind_cap: ", round(end.time - start.time, 2), "seconds\n")
summary(m_wind_cap)

start.time <- Sys.time()
m_Wind_windPA <- gam(Wind_windPA, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
cat("Time for m_Wind_windPA: ", round(end.time - start.time, 2), "seconds\n")
summary(m_Wind_windPA)

start.time <- Sys.time()
m_wind_rsa <- gam(Wind_windrsa, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
cat("Time for m_wind_rsa: ", round(end.time - start.time, 2), "seconds\n")
summary(m_wind_rsa)

start.time <- Sys.time()
m_wind_height <- gam(Wind_windheight, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
cat("Time for m_wind_height: ", round(end.time - start.time, 2), "seconds\n")
summary(m_wind_height)

start.time <- Sys.time()
m_wind_age <- gam(Wind_windage, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
cat("Time for m_wind_age: ", round(end.time - start.time, 2), "seconds\n")
summary(m_wind_age)

start.time <- Sys.time()
m_wind_base <- gam(Wind_base, data = spec_split_spatialgrain$train, family = "nb", knots = time_knots)
end.time <- Sys.time()
cat("Time for m_wind_base: ", round(end.time - start.time, 2), "seconds\n")
summary(m_wind_base)

MuMIn::AICc(m_wind_count, m_wind_cap, m_Wind_windPA, m_wind_rsa, m_wind_height, m_wind_age, m_wind_base)

##### Part 10: Top Model Diagnostics
# Variance Components (rescaled)
gam.vcomp(m_wind, rescale = TRUE)

# Model ANOVA
anova(m_wind)

# Default Diagnostic Plot
plot(m_wind)


##### Part 11: Figure Generation
# Extract scaling parameters for continuous wind variable (assumes standardized inputs)
scaling_params_sub <- scaling_params %>%
  filter(variable %in% c("wind_variable"))

# Create unstandardized variable in test set for plotting
spec_split_spatialgrain$test$wind_variable_unscaled <- spec_split_spatialgrain$test$wind_variable *
  scaling_params_sub$sd[scaling_params_sub$variable == "wind_variable"] +
  scaling_params_sub$mean[scaling_params_sub$variable == "wind_variable"]

# 1. Continuous Wind Variable Only Figure
pred_wind_nl <- predict.gam(m_wind_nonlinear, newdata = spec_split_spatialgrain$test,
                            type = "response", se.fit = TRUE, allow.new.levels = TRUE)

pred_df_wind_nl <- spec_split_spatialgrain$test
pred_df_wind_nl$predicted <- pred_wind_nl$fit

ggplot(pred_df_wind_nl, aes(x = wind_variable, y = predicted)) +
  geom_smooth() +
  labs(x = "Wind Variable", y = "Predicted Abundance") +
  theme_minimal()

# 2. Categorical Wind Variable Only Figure
# Predict and plot from categorical wind-only model
pred_wind_cat <- predict.gam(m_wind_categorical, newdata = spec_split_spatialgrain$test,
                             type = "response", se.fit = TRUE, allow.new.levels = TRUE)

pred_df_wind_cat <- spec_split_spatialgrain$test
pred_df_wind_cat$predicted <- pred_wind_cat$fit

ggplot(pred_df_wind_cat, aes(x = wind_variable, y = predicted)) +
  geom_boxplot() +
  labs(x = "Wind Variable", y = "Predicted Abundance") +
  theme_minimal()
