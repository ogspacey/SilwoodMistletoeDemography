### Non-discrete IPM without stage distinctions
### Oliver G. Spacey, Owen R. Jones, Sydne Record, Arya Y. Yue, Wenyi Liu, Alice Rosen, Michael Crawley, Chris J. Thorogood, Roberto Salguero-Gómez
### 02.04.2025

# Pre-amble ---------------------------------------------------------------
# Clear environment
rm(list = ls())

# Load packages - install if necessary
library(ggplot2)      # Data visualisation
library(patchwork)    # Plot arrangement
library(tidyverse)    # Data wrangling
library(lme4)         # Regression models
library(lmerTest)     # Regression models
library(RColorBrewer) # Colours for visualisation
library(khroma)       # Colours for visualisation
library(popbio)       # Demographic analysis
library(Rage)         # Demographic analysis
library(reshape2)     # Data restructuring
library(ggnewscale)   # Data visualisation
library(scales)       # Data visualisation
library(cowplot)      # Data visualisation
library(formattable)  # Data visualisation

# Set working directory as required

# Set colours for visualisation
colours <- color("vibrant")

# Load data ---------------------------------------------------------------
# Load mistletoe demographic data
mst_df_raw <- read.csv("mistletoe_data_raw.csv")

# Load host metadata
hst_df_raw <- read.csv("host_metadata_raw.csv")

# Wrangle data ------------------------------------------------------------
# Remove hosts from host data frame in which mistletoes were not measured
hst_df_raw <- filter(hst_df_raw, Host_ID %in% unique(mst_df_raw$Host_ID))

# Create genus variable for hosts
hst_df_raw <- mutate(hst_df_raw, Genus = word(Genus_spp))

# Count number of hosts
nrow(hst_df_raw)
# 24 host trees

# Combine mistletoe demographic data and host data into master dataframe
mst_hst_df <- full_join(mst_df_raw, hst_df_raw, by = "Host_ID")

# Create individual mistletoe ID combining mistletoe and host IDs
mst_hst_df <- mutate(mst_hst_df, Indiv_ID = paste(Mistletoe_ID, Host_ID, sep = "_"))

# Remove columns measuring perimeter as area is a better indicator of size
mst_hst_df <- select(mst_hst_df, -contains('Perim'))

# Remove mistletoes that are never seen
mst_hst_df <- filter(mst_hst_df, Notes != "Never seen")

# Remove record from clump containing multiple individual (ID is longer than 4 characters)
mst_hst_df <- filter(mst_hst_df, nchar(Mistletoe_ID) == 4)

# Calculate number of mistletoe individuals measured
nrow(mst_hst_df)
# 740 mistletoes

# Calculate number of observations - where mistletoes either first censused (FC), survived (Surv) or dead (Dead)
sum(mst_hst_df == "FC", na.rm = TRUE) + sum(mst_hst_df == "Surv", na.rm = TRUE) + sum(mst_hst_df == "Dead", na.rm = TRUE)
# 3394 observations made

# Count number of times mistletoes ("FC" or "Surv") were not measured (NA) but considered present - too blurry or obscured
known_presences <- 0
NAs <- 0
mst_hst_df$Missed <- 0
for(i in 1:nrow(mst_hst_df)) {
  for(j in 1:(ncol(mst_hst_df)-1)){
    if(mst_hst_df[i,j] %in% c("FC", "Surv") == TRUE) {
      known_presences <- known_presences + 1
      if(is.na(mst_hst_df[i,j+1]) == TRUE){
        NAs <- NAs + 1
        mst_hst_df[i,"Missed"] <- mst_hst_df[i,"Missed"] + 1
      }
    }
    else { count <- count }
  }
}
nrow(subset(mst_hst_df, Missed == 0))
nrow(subset(mst_hst_df, Missed == 0))/nrow(mst_hst_df)
# 454 NAs out of 3183 detections affecting 430 out of 740 (58.1%) mistletoes in at least one year

# Estimate detectability
NAs/known_presences
# 14.3% of known presences not measured

# Estimate proportion of mistletoes for which height is calculated directly
nrow(filter(mst_hst_df, Height > 0))/nrow(mst_hst_df)
# 39.1% of mistletoe heights measured directly

# Create figure to show observation history
obs_area_df <- select(mst_hst_df, starts_with("Area"))
obs_status_df <- select(mst_hst_df, starts_with("Status"))

# Create Fig_S1 data frame
fig_s1_df <- obs_area_df
# Assign values where area measured
fig_s1_df[!is.na(fig_s1_df)] <- "black"
# Assign values where not yet seen
fig_s1_df[obs_status_df == "NYS"] <- "white"
# Assign values where dead and after dead
fig_s1_df[obs_status_df == "Dead"] <- "Dead"

# Function to make after dead white row-wise
propagate_white <- function(row) {
  dead_pos <- which(row == "Dead")
  if (length(dead_pos) > 0) {
    min_white <- min(dead_pos)
    row[min_white:length(row)] <- "white"
  }
  return(row)
}

# Apply row-wise
fig_s1_df[] <- t(apply(fig_s1_df, 1, propagate_white))

# Anything remaining was not measured but was not dead or NYS, so assign as grey
fig_s1_df[is.na(fig_s1_df)] <- "lightgrey"

# Make names years
names(fig_s1_df) <- as.character(c(2014:2023))

# Convert to long format for ggplot
fig_s1_df_long <- fig_s1_df %>%
  mutate(row = nrow(.) : 1) %>%  # Flip y-axis so first row is at the top
  pivot_longer(-row, names_to = "col", values_to = "fill_color") %>%
  mutate(col = factor(col, levels = names(fig_s1_df)))  # Keep column order

# Plot with geom_tile
ggplot(fig_s1_df_long, aes(x = col, y = row, fill = fill_color)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  scale_fill_identity() +
  scale_x_discrete(name = NULL, position = "top") +  # Labels on top
  scale_y_continuous(name = NULL, breaks = NULL) +   # Remove row numbers
  theme_void() +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(size = 12, face = "bold")
  )


## Correct areas for distance and skew ----------------------------------------------------
# Estimate distance to mistletoes (hypotenuse) for mistletoes with imputed height via Pythagoras, to rectify size perception
mst_hst_df <- mutate(mst_hst_df, 
                     Estimated_hypotenuse = sqrt((mst_hst_df$Imputed_height)^2 + (mst_hst_df$Host_horizontal_distance)^2))

# Make "#VALUE!" NA for height variable and make numeric
mst_hst_df$Height[mst_hst_df$Height == "#VALUE!"] <- NA
mst_hst_df$Height <- as.numeric(mst_hst_df$Height)

# For each mistletoe, combine field estimates of height, horizontal distance and hypotenuse,
# With those measured digitally, estimating horizontal distance to mistletoe as distance to the host for mistletoes where it was not measured
mst_hst_df <- mutate(mst_hst_df, 
                     Combined_hypotenuse = case_when(is.na(Hypotenuse) == FALSE ~ Hypotenuse,
                                                     is.na(Hypotenuse) == TRUE ~ Estimated_hypotenuse),
                     Combined_height = case_when(is.na(Height) == FALSE ~ Height,
                                                 is.na(Height) == TRUE ~ Imputed_height),
                     Combined_hor_dist = case_when(is.na(Horizontal) == FALSE ~ Horizontal,
                                                   is.na(Horizontal) == TRUE ~ Host_horizontal_distance))

# Check relationship between area and height before correction
# Make long format so there is a single column for areas
area_long_df <- mst_hst_df %>%
  pivot_longer(cols = starts_with("Area"), 
               names_to = "Year", 
               values_to = "Area") %>%
  mutate(Year = as.numeric(sub("Area", "", Year)) + 2000)

# Plot height against area
ggplot(data = area_long_df, aes(x = Area, y = Combined_height)) +
  geom_point() +
  theme_bw()

# Plot height against log(area) to reduce skew
ggplot(data = area_long_df, aes(x = log(Area), y = Combined_height)) +
  geom_point() +
  theme_bw()

# Mixed effect model of area and height with individual ID as random effect
area_height_mem <- lmer(Combined_height ~ log(Area) + (1 | Indiv_ID), data = area_long_df)
summary(area_height_mem)
# Significant positive relationship between height and area

# Rectify size based on horizontal distance to base of tree (standard)
# If an object is twice as far away relative to the standard, it will appear twice as small in length, and four times as small in area
# Calculate the area adjusting for distance (hypotenuse) relative to the standard (horizontal distance) as Adjusted area = (Distance to object/Distance to standard)^2 * Measured area

# Select columns of areas, hypotenuses and horizontal distance to standard (host), and Indiv_ID
adjusted_df <- dplyr::select(mst_hst_df, Indiv_ID, starts_with("Area"), Combined_hypotenuse, Host_horizontal_distance)

# For M001_H0142, height was measured directly but horizontal distance was not
# It is at the same height as the observer, so assume hypotenuse/height = 1
# As a placeholder, set hypotenuse = height = 1
adjusted_df[adjusted_df$Indiv_ID == "M001_H0142",]$Combined_hypotenuse <- 1
adjusted_df[adjusted_df$Indiv_ID == "M001_H0142",]$Host_horizontal_distance <- 1

# Adjust sizes
adjusted_df[,2:11] <- adjusted_df[,2:11] * (adjusted_df$Combined_hypotenuse/adjusted_df$Host_horizontal_distance)^2

# Stack area columns
stacked_areas_df <- pivot_longer(adjusted_df, cols = starts_with("Area"), names_to = "Year", values_to = "Area")

# Visualise areas in histogram 
ggplot(data = stacked_areas_df, aes(x = Area)) +
  geom_histogram() +
  theme_bw()
# Highly skewed distribution

# Log-transform to reduce skew and re-plot areas
ggplot(data = stacked_areas_df, aes(x = log(Area))) +
  geom_histogram() +
  theme_bw()
# More normally distributed, less skewed
# No clear effect of particular years

# Perform Shapiro-Wilk test for normality
shapiro.test(log(stacked_areas_df$Area))
# Reject null that log(Areas) is normally distributed
# Even though not normally distributed, log transformation facilitates plotting and reduces skew

# Log-transform areas
log_area_df <- adjusted_df
log_area_df[,2:11] <- log(log_area_df[,2:11])
colnames(log_area_df)[2:11] <- paste0("log", colnames(log_area_df)[2:11])

## Convert to long format --------------------------------------------------
# Convert table to long format to show status, area and fruit in time t0 and status, area and fruit in time t1
# This will be used for models across a time step, e.g., survival and growth
# Use log-transformed, adjusted growth rates
years <- c(14:22) # vector of years to run for loop through
long_df_ls <- list()   # list of long data frames to store each year
for(year in years){ # run through years sequentially
  Status_t0  <- paste("Status", year, sep = "") 
  logArea_t0 <- paste("logArea", year, sep = "")
  Fruit_t0   <- paste("Fruit", year, sep = "")
  Status_t1  <- paste("Status", year + 1, sep = "")
  logArea_t1 <- paste("logArea", year + 1, sep = "")
  Fruit_t1   <- paste("Fruit", year + 1, sep = "")
  mst_df_year <- data.frame(Host_ID = mst_hst_df$Host_ID,
                            Mistletoe_ID = mst_hst_df$Mistletoe_ID,
                            Indiv_ID = mst_hst_df$Indiv_ID, 
                            Year_t0 = year,
                            Status_t0 = mst_hst_df[[Status_t0]],
                            logArea_t0 = log_area_df[[logArea_t0]],
                            Fruit_t0 = mst_hst_df[[Fruit_t0]],
                            Year_t1 = year + 1,
                            Status_t1 = mst_hst_df[[Status_t1]],
                            logArea_t1 = log_area_df[[logArea_t1]],
                            Fruit_t1 = mst_hst_df[[Fruit_t1]],
                            Height = mst_hst_df$Combined_height
  )
  long_df_ls[[year]] <- mst_df_year
}

# Merge rows to create long data frame
long_df <- bind_rows(long_df_ls)

# Create long data frame for each year separately (just t0)
# This will be used for models within one year, e.g., reproduction
# Use log-transformed, adjusted growth rates
years <- c(14:23) # vector of years to run for loop through
long_df_ls <- list()   # list of long data frames to store each year
for(year in years){ # run through years sequentially
  Status  <- paste("Status", year, sep = "") 
  logArea <- paste("logArea", year, sep = "")
  Fruit   <- paste("Fruit", year, sep = "")
  mst_df_year <- data.frame(Host_ID = mst_hst_df$Host_ID,
                            Indiv_ID = mst_hst_df$Indiv_ID, 
                            Year = year,
                            Status = mst_hst_df[[Status]],
                            logArea = log_area_df[[logArea]],
                            Fruit = mst_hst_df[[Fruit]],
                            Height = mst_hst_df$Combined_height
  )
  long_df_ls[[year]] <- mst_df_year
}


# Merge rows to create long data frame
long_by_yr_df <- bind_rows(long_df_ls)

## Calculate intensities ------------------------------------------------
# Sum intensity (parasite load) for each host each year
# Create list of host IDs
hosts <- hst_df_raw$Host_ID

# Create empty data frame with only Host IDs, columns for intensities
intens_df <- data.frame(Host_ID = hst_df_raw$Host_ID)
row.names(intens_df) <- hst_df_raw$Host_ID
intens_df[c("I_14", "I_15", "I_16", "I_17", "I_18",
            "I_19", "I_20", "I_21", "I_22", "I_23")] <- NA
years <- as.character(14:23)
for(year in years){
  for(host in hosts){
    filtered_df <- filter(mst_hst_df, Host_ID == host) %>%
      dplyr::select(ends_with(year))
    intensity <- length(which(filtered_df[1] == "FC" | filtered_df[1] == "Surv"))
    I_col_name <- paste("I_", year, sep = "")
    intens_df[host, I_col_name] <- intensity
  }
}

# Calculate how many mistletoes seen each year
colSums(intens_df[,-1])

# Mean number of mistletoes seen each year - mean population size
mean_pop_size <- mean(colSums(intens_df[,-1]))

# Calculate maximum parasite load observed for each tree
intens_df$Max_I <- pmax(intens_df$I_14,
                        intens_df$I_15,
                        intens_df$I_16,
                        intens_df$I_17,
                        intens_df$I_18,
                        intens_df$I_19,
                        intens_df$I_20,
                        intens_df$I_21,
                        intens_df$I_22,
                        intens_df$I_23)

# Add host metadata to PL data frame
intens_df <- full_join(intens_df, hst_df_raw)

# Separate stages ---------------------------------------------------------
# Separate stages based on size at minimum reproduction
# This will only be used for regressions where only adults are being considered, and for designating recruits
# Calculate minimum size at reproduction
min_size_rep <- min(filter(long_by_yr_df, Fruit == 1)$logArea, na.rm = TRUE)

# Assign individuals as juveniles with less than min size at reproduction, adults as larger
wrangled_by_yr_df <- mutate(long_by_yr_df, 
                            Stage = case_when(logArea < min_size_rep  ~ "Juvenile",
                                              logArea >= min_size_rep ~ "Adult"))
wrangled_df <- mutate(long_df, 
                      Stage = case_when(logArea_t0 < min_size_rep  ~ "Juvenile",
                                        logArea_t0 >= min_size_rep ~ "Adult"))

# Explore data --------------------------------------------------------
# Visualise distribution of mistletoe heights
Fig_S2b <- ggplot(data = mst_hst_df, aes(x = Combined_height)) +
  geom_histogram() +
  labs(x = "Height on the host (m)", y = "") +
  theme_bw()
Fig_S2b

# Test for normality
shapiro.test(mst_hst_df$Combined_height)
# Not normally distributed

# Log-transform and test
shapiro.test(log(mst_hst_df$Combined_height))
# Normality reduced after log-transformation - keep raw values

# Calculate heights of mistletoe relative to tree
mst_hst_df <- mutate(mst_hst_df, 
                     Relative_height = Combined_height/Host_height)

# Visualise relative heights in histogram
ggplot(data = mst_hst_df, aes(x = Relative_height)) +
  geom_histogram() +
  theme_bw()
# Some estimated at >1 due to parallax with mistletoes, making them appear higher than the top of the tree
# Also, raw heights are a better proxy for resource access
# Use raw heights instead

# Visualise relationship between host height and mistletoe height
Fig_S4a <- ggplot(data = mst_hst_df, aes(x = Host_height, y = Combined_height)) +
  geom_point(aes(col = Genus)) +
  theme_bw() +
  geom_smooth(method = "lm", col = "darkgrey") +
  labs(x = "Host height (m)", y = "Mistletoe height on the host (m)") +
  theme(legend.position = "none", 
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))
Fig_S4a

# Visualise relationship between host height and mistletoe intensity
Fig_S4b <- ggplot(data = intens_df, aes(x = Host_height, y = Max_I)) +
  geom_point(aes(col = Genus)) +
  theme_bw() +
  geom_smooth(method = "lm", col = "darkgrey") +
  labs(x = "Host height (m)", y = "Maximum mistletoe intensity") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15, face = "italic"))
Fig_S4b
# Taller the host, higher up mistletoes generally are

# Plot Figure S4
Fig_S4 <- Fig_S4a + Fig_S4b + 
  plot_annotation(tag_levels = "a") +
  theme(plot.tag = element_text(size = 15))
Fig_S4

# Visualise relationship between host and mistletoe relative heights
ggplot(data = mst_hst_df, aes(x = Host_height, y = Relative_height)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(x = "Host height", y = "Relative mistletoe height")
# But mistletoe height does not increase relative to the host tree height

# Visualise distribution of host heights and species
ggplot(data = hst_df_raw, aes(x = Host_height, fill = Genus)) +
  geom_histogram() +
  theme_bw()
# Highly collinear - Tilia are tallest trees and majority

# Visualise distribution of mistletoe heights and individual hosts
ggplot(data = mst_hst_df, aes(x = Host_ID, y = Combined_height, col = Genus)) +
  geom_point() +
  theme_bw()
# Highly collinear - hard to separate effects

# Visualise distribution of host height, species and maximum intensity
ggplot(data = intens_df, aes(x = Host_height, y = Max_I, col = Genus)) +
  geom_point() +
  theme_bw()
# Any size tree can have low intensity, but only huge trees can have high intensity

# Visualise distribution of host species and maximum intensity
ggplot(data = intens_df, aes(x = Max_I, fill = Genus)) +
  geom_histogram() +
  theme_bw()
# Intensity overlaps slightly between genera

# Calculate mean and SD of log(Area)
mean_logArea <- mean(wrangled_by_yr_df$logArea, na.rm = TRUE)
sd_logArea <- sd(wrangled_by_yr_df$logArea, na.rm = TRUE)

# Visualise distribution of mistletoe areas
ggplot(data = wrangled_by_yr_df, aes(x = logArea)) +
  geom_histogram() +
  geom_vline(xintercept = mean_logArea - 3 * sd_logArea) +
  geom_vline(xintercept = mean_logArea + 3 * sd_logArea) +
  theme_bw()

# Remove outliers +- 3SD
wrangled_by_yr_df <- filter(wrangled_by_yr_df, 
                            logArea > mean_logArea - 3 * sd_logArea & 
                              logArea < mean_logArea + 3 * sd_logArea)

wrangled_df <- filter(wrangled_df, 
                      logArea_t0 > (mean_logArea - 3 * sd_logArea) & 
                        logArea_t0 < (mean_logArea + 3 * sd_logArea))

# Create new x label
xl <- expression(log ~ "area" ~ cm^2)

# Re-visualise distribution of mistletoe areas
Fig_S2a <- ggplot(data = wrangled_by_yr_df, aes(x = logArea)) +
  geom_histogram() +
  labs(x = xl, y = "Frequency") +
  theme_bw()
Fig_S2a

# Plot Figure S2
Fig_S2 <- Fig_S2a + Fig_S2b + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10))
Fig_S2

# Test for normality in area
shapiro.test(wrangled_by_yr_df$logArea)
# Still not normal

# Outliers of area probably due to measurement error, unlike height outliers
# Do not remove outliers for height

# Visualise relationship between areas and height
Fig_S3 <- ggplot(data = wrangled_by_yr_df, aes(x = logArea, y = Height)) +
  geom_point() +
  labs(x = xl, y = "Height on the host (m)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))
Fig_S3

# Re-test relationship between height and area

# Mixed effect model of area and height  with individual ID as random effect
area_height_mem <- lmer(Height ~ logArea + (1 | Indiv_ID), data = wrangled_by_yr_df)
summary(area_height_mem)
# No longer significant relationship between height and area after correction

# Construct and select vital rate regressions --------------------------------------------------------
# Create function to compare models for a particular vital rate, creating empty lists
calls <- list()
sig <- list()
est_p <- list()
n_obs <- list()
n_ind <- list()
aics <- list()
mod.sel.table <- function(mods){
  for(i in 1:length(mods)){
    est_ps <- list()
    sigs <- list()
    calls[i] <- as.character(mods[[i]]@call)[2]
    summary_tb <- summary(mods[[i]])$coefficients
    for(j in 2:nrow(summary_tb)){
      est_ps[j] <- paste(row.names(summary_tb)[j], 
                         "β=", round(summary_tb[j,1], digits = 3), 
                         ", P=", round(summary_tb[j,ncol(summary_tb)], digits = 3),
                         " ")
      if(summary_tb[j, ncol(summary_tb)] <= 0.05){
        sigs[j] <- row.names(summary_tb)[j]
      } else { sigs[j] <- "" }
    }
    est_p[i] <- paste(unlist(est_ps), collapse = " ")
    sig[i] <- paste(unlist(sigs), collapse = " ")
    n_obs[i] <- as.numeric(summary(mods[[i]])$devcomp$dims[1])
    n_ind[i] <- as.numeric(summary(mods[[i]])$ngrps)
    aics[i] <- round(AIC(mods[[i]]), digits = 2)
    mod_sel_table <- data.frame(Model = unlist(calls),
                                Parameter_estimates_p = unlist(est_p),
                                Significant = unlist(sig),
                                N_obs = unlist(n_obs),
                                N_ind = unlist(n_ind),
                                AIC = unlist(aics))
    mod_sel_table$Significant <- trimws(mod_sel_table$Significant)
    mod_sel_table <- data.frame(lapply(mod_sel_table, function(x){
      gsub("P= 0 ", "P<0.001 ", x)
    }))
  }
  return(mod_sel_table) 
}

## Model survival ----------------------------------------------------------
# Create binary variable for survival 
wrangled_df <- mutate(wrangled_df, Survive = case_when(
  Status_t1 == "Dead" ~ 0,
  Status_t1 == "Surv" ~ 1
))

# Plot Survival ~ log(Area)
ggplot(data = wrangled_df, aes(x = logArea_t0, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ log(Area) with individual ID as random effect
sur_area_glmm <- glmer(Survive ~ logArea_t0 + (1 | Indiv_ID), 
                       data = wrangled_df, family = binomial)
summary(sur_area_glmm)
# No significant effect of log(Area)

# Figure 1a - Plot Survival ~ log(Area)
Fig_1a <- ggplot(data = wrangled_df, aes(x = logArea_t0, y = Survive)) +
  geom_point(size = 1, alpha = 0.5, col = colours(6)[1]) +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial), 
              col = colours(6)[1]) +
  labs(x = "", y = "Survival probability (t to t+1)") +
  theme_bw() + 
  theme(axis.title.y = element_text(size = 10))
Fig_1a

# Plot Survival ~ Height
ggplot(data = wrangled_df, aes(x = Height, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "Height above ground (m)", y = "Probability of survival to t1") +
  theme_bw()

# Model Survival ~ Height with individual ID as random effect
sur_ht_glmm <- glmer(Survive ~ Height + (1 | Indiv_ID), 
                     data = wrangled_df, family = binomial)
summary(sur_ht_glmm)
# No significant effect of Height

# Figure 1b - Plot Survival ~ Height
Fig_1b <- ggplot(data = wrangled_df, aes(x = Height, y = Survive)) +
  geom_point(size = 1, alpha = 0.5, col = colours(6)[1]) +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial),
              col = colours(6)[1]) +
  labs(x = "", y = "") + 
  theme_bw()
Fig_1b

# Model Survival ~ log(Area) + Height with individual ID as random effect
sur_area_ht_glmm <- glmer(Survive ~ logArea_t0 + Height + (1 | Indiv_ID), 
                          data = wrangled_df, family = binomial)
summary(sur_area_ht_glmm)
# No significant effect of either Area or Height

# Model Survival ~ log(Area) * Height with individual ID as random effect
sur_area_ht_int_glmm <- glmer(Survive ~ logArea_t0 * Height + (1 | Indiv_ID), 
                          data = wrangled_df, family = binomial)
summary(sur_area_ht_int_glmm)
# Significant interaction between area and height, BUT has singular fit, so cannot interpret accurately - ignore

# List models for survival
sur_mods <- c(sur_area_glmm,
              sur_ht_glmm,
              sur_area_ht_glmm)

# Create model selection table for survival
sur_mod_sel <- bind_cols(data.frame(Vital_rate = "Survival"),
                         mod.sel.table(sur_mods))
# Neither area nor height predict survival across all mistletoes
# For IPM, use Survival ~ logArea as simplest model which incorporates the common state variable

## Model growth ------------------------------------------------------------
# Plot log(Area) in time t1 against log(Area) in time t0
ggplot(data = wrangled_df, aes(x = logArea_t0, y = logArea_t1)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  lims(x = c(0, max(long_df$logArea_t0)), y = c(0, max(long_df$logArea_t1))) +
  geom_smooth(method = "lm") +
  theme_bw()

# Calculate per-time-step growth rate for each individual in each time step
wrangled_df$Growth_t0_t1 <- (wrangled_df$logArea_t1 - wrangled_df$logArea_t0)/wrangled_df$logArea_t0

# Plot relative growth from t0 to t1 against log(Area) in time t0
ggplot(data = wrangled_df, aes(x = logArea_t0, y = Growth_t0_t1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() 

# Model Growth ~ log(Area) with Indiv_ID as random effect
gro_area_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + (1 | Indiv_ID), 
                     data = wrangled_df)
summary(gro_area_mem)

# Plot relative growth from t0 to t1 against log(Area) with quadratic function in time t0
ggplot(data = wrangled_df, aes(x = logArea_t0, y = Growth_t0_t1)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  theme_bw() 

# Model Growth ~ log(Area) with quadratic term and Indiv_ID as random effect
gro_area_quad_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + I(logArea_t0^2) + (1 | Indiv_ID), 
                          data = wrangled_df)
summary(gro_area_quad_mem)

# Figure 1c - Plot Growth ~ log(Area)^2 + log(Area)
Fig_1c <- ggplot(data = wrangled_df, aes(x = logArea_t0, y = Growth_t0_t1)) +
  geom_point(size = 1, alpha = 0.2, col = colours(6)[6]) +
  geom_smooth(data = wrangled_df, method = "lm", formula = y ~ poly(x, 2),
              col = colours(6)[6]) +
  labs(x = "", y = "RGR (t to t+1)") +
  theme_bw()
Fig_1c

# Plot relative growth from t0 to t1 against Height
ggplot(data = wrangled_df, aes(x = Height, y = Growth_t0_t1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() 

# Model Growth ~ Height with Indiv_ID as random effect
gro_ht_mem <- lmer(Growth_t0_t1 ~ Height + (1 | Indiv_ID), 
                   data = wrangled_df)
summary(gro_ht_mem)
# Singular fit - ignore

# Figure 1d - Plot Growth ~ Height
Fig_1d <- ggplot(data = wrangled_df, aes(x = Height, y = Growth_t0_t1)) +
  geom_point(size = 1, alpha = 0.2, col = colours(6)[6]) +
  geom_smooth(data = wrangled_df, method = "lm", formula = y ~ x,
              col = colours(6)[6]) +
  labs(x = "", y = "") +
  theme_bw()
Fig_1d

# Model Growth ~ log(Area) + Height with Indiv_ID as random effect
gro_area_ht_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + Height + (1 | Indiv_ID), 
                        data = wrangled_df)
summary(gro_area_ht_mem)

# Model Growth ~ Height + Stage with Indiv_ID as random effect
gro_ht_stage_mem <- lmer(Growth_t0_t1 ~ Height + Stage + (1 | Indiv_ID), 
                         data = wrangled_df)
summary(gro_ht_stage_mem)

# Model Growth ~ log(Area) * Height with Indiv_ID as random effect
gro_area_ht_int_mem <- lmer(Growth_t0_t1 ~ logArea_t0 * Height + (1 | Indiv_ID), 
                              data = wrangled_df)
summary(gro_area_ht_int_mem)

# List models for growth
gro_mods <- c(gro_area_mem,
              gro_area_quad_mem,
              gro_area_ht_mem,
              gro_area_ht_int_mem)

# Create model selection table for growth
gro_mod_sel <- bind_cols(data.frame(Vital_rate = "Growth"),
                         mod.sel.table(gro_mods))
# Quadratic function of area best to predict growth across all mistletoes
# Find difference in AIC for quadratic vs linear model
AIC(gro_area_quad_mem) - AIC(gro_area_mem)

## Model fruiting ------------------------------------------------------
# Figure 1e
# Plot fruiting ~ logArea
Fig_1e <- ggplot(data = wrangled_by_yr_df, aes(x = logArea, y = Fruit)) +
  geom_point(size = 1, alpha = 0.5, col = colours(6)[3]) +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial),
              col = colours(6)[3]) +
  labs(x = xl, y = "Fruiting probability (t)") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 10))
Fig_1e

# Model fruiting ~ logArea with Indiv_ID as random effect
fru_area_glmm <- glmer(Fruit ~ logArea + (1 | Indiv_ID), 
                       data = wrangled_by_yr_df, family = binomial)
summary(fru_area_glmm)

# Plot fruiting ~ Height
ggplot(data = wrangled_by_yr_df, aes(x = Height, y = Fruit)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "Height above ground (m)", y = "Probability of fruiting in t0") +
  theme_bw()

# Figure 1f
# Plot fruiting ~ Height
Fig_1f <- ggplot(data = wrangled_by_yr_df, aes(x = Height, y = Fruit)) +
  geom_point(size = 1, alpha = 0.5, col = colours(6)[3]) +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial),
              col = colours(6)[3]) +
  labs(x = "Height above ground (m)", y = "") +
  theme_bw()
Fig_1f

# Model fruiting ~ Height with Indiv_ID as random effect
fru_ht_glmm <- glmer(Fruit ~ Height + (1 | Indiv_ID), 
                     data = wrangled_by_yr_df, family = binomial)
summary(fru_ht_glmm)

# Plot fruiting ~ log(Area) + Height
ggplot(data = wrangled_by_yr_df, aes(x = logArea, y = Fruit, col = Height)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "log(Area) in t0", y = "Probability of fruiting in t0") +
  theme_bw()

# Model fruiting ~ logArea + Height with Indiv_ID as random effect
fru_area_ht_glmm <- glmer(Fruit ~ logArea + Height + (1 | Indiv_ID), 
                          data = wrangled_by_yr_df, family = binomial)
summary(fru_area_ht_glmm)

# Model fruiting ~ logArea * Height with Indiv_ID as random effect
fru_area_ht_int_glmm <- glmer(Fruit ~ logArea * Height + (1 | Indiv_ID), 
                              data = wrangled_by_yr_df, family = binomial)
summary(fru_area_ht_int_glmm)
# Does not converge - ignore model

# List models for fruiting
fru_mods <- c(fru_area_glmm,
              fru_ht_glmm,
              fru_area_ht_glmm)

# Create model selection table for fruiting
fru_mod_sel <- bind_cols(data.frame(Vital_rate = "Fruiting"),
                         mod.sel.table(fru_mods))
# Fruiting best modelled as additive model of area and height

## Overall model selection -------------------------------------------------
# Table S1
# Combine all model selection tables
mod_sel <- bind_rows(sur_mod_sel, gro_mod_sel, fru_mod_sel)

# Find minimum AIC values
min_AIC <- mod_sel %>%
  group_by(Vital_rate) %>%
  filter(AIC == case_when(AIC > 0 ~ min(AIC),
                          AIC < 0 ~ max(AIC))) %>%
  select(AIC, Vital_rate)

# Show which models selected
mod_sel <- mod_sel %>%
  mutate(Lowest_AIC = case_when(AIC %in% min_AIC$AIC ~ TRUE, 
                              .default = FALSE))

# Export as .csv
write.csv(mod_sel, "Table_S1.csv")


# Figure 1 ----------------------------------------------------------------
# Combine all vital rate regressions into Figure 1
Fig_1 <- (Fig_1a + Fig_1b) / (Fig_1c + Fig_1d) / (Fig_1e + Fig_1f) + 
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10))
Fig_1


# Trade-offs --------------------------------------------------------------
# Create data frame to store models of vital rate trade-offs
vr_to_df <- data.frame(VR_1 = NA,
                       VR_2 = NA,
                       Model = NA, 
                       Est_p = NA,
                       N_obs = NA,
                       N_ind = NA)

# To measure trade-offs with future survival, create extra variable for survival from t+1 to t+2
wrangled_to_df <- wrangled_df %>%
  arrange(Indiv_ID, Year_t0) %>%
  group_by(Indiv_ID) %>%
  mutate(Survive_t1_t2 = lead(Survive, 1)) %>% # Shift the Survive column forward by 1 row within each individual
  ungroup()

# Rename Survive to Survive_t0_t1 for clarity
wrangled_to_df <- wrangled_to_df %>%
  rename(Survive_t0_t1 = Survive)

# Replace Survive_t1_t2 with NA for cases where Year_t1 is the last year of measurement
wrangled_to_df <- wrangled_to_df %>%
  mutate(Survive_t1_t2 = ifelse(Year_t1 == max(Year_t0), NA, Survive_t1_t2))

# Plot survival from t+1 to t+2 as a function of growth from t to t+1 for all individuals
ggplot(data = wrangled_to_df, aes(x = Growth_t0_t1, y = Survive_t1_t2)) +
  geom_point() +
  stat_smooth(method = "glm", method.args = list(family = binomial)) +
  labs(x = "", y = "Survival probability (t+1 to t+2)") +
  theme_bw()

# Model survival from t+1 to t+2 as a function of growth from t to t+1 for all individuals with Indiv_ID as random effect
sur_gro_to_glmm <- glmer(Survive_t1_t2 ~ Growth_t0_t1 + (1 | Indiv_ID), 
                         data = wrangled_to_df, family = binomial)
summary(sur_gro_to_glmm)
# When individual ID controlled for, no significant relationship - is singular, cannot measure

# Examine trade-offs involving reproduction for adults
# Filter wrangled data for adults known to be female
wrangled_fem_to_df <- filter(wrangled_to_df, Stage == "Adult") %>%
  group_by(Indiv_ID) %>%
  filter(any(Fruit_t0 == 1)) %>%
  ungroup()

# Plot survival from t to t+1 as a function of fruiting in t for females
ggplot(data = wrangled_fem_to_df, aes(x = Fruit_t0, y = Survive_t0_t1)) +
  geom_jitter(width = 0.2, height = 0.1, alpha = 0.3) +
  stat_smooth(method = "glm", method.args = list(family = binomial)) +
  labs(x = "", y = "Survival probability (t to t+1)") +
  theme_bw()

# Model survival from t to t+1 as a function of fruiting in t for females with Indiv_ID as random effect
sur_fru_to_glmm <- glmer(Survive_t0_t1 ~ Fruit_t0 + (1 | Indiv_ID), 
                         data = wrangled_fem_to_df, family = binomial)
summary(sur_fru_to_glmm)
# Significant negative relationship between fruiting and survival
# Singular fit means cannot interpret

# Plot fruiting in t+1 as a function of growth from t to t+1 for adults
Fig_2 <- ggplot(data = wrangled_fem_to_df, aes(x = Growth_t0_t1, y = Fruit_t1)) +
  geom_point(col = colours(6)[2]) +
  stat_smooth(method = "glm", method.args = list(family = binomial), col = colours(6)[2]) +
  labs(x = "RGR (t to t+1)", y = "Fruiting probability (t+1)") +
  lims(x = c(-0.15, max(wrangled_fem_to_df$Growth_t0_t1, na.rm = TRUE))) +
  theme_bw()
Fig_2

# Model fruiting in t+1 as a function of growth from t to t+1 for adults with Indiv_ID as random effect
fru_gro_to_glmm <- glmer(Fruit_t1 ~ Growth_t0_t1 + (1 | Indiv_ID), 
                         data = wrangled_fem_to_df, family = binomial)
summary(fru_gro_to_glmm)
# significant negative relationship; more growth, less likely to fruit
# reflects investment in reproduction as mistletoes grow larger

# Plot growth from t to t+1 as a function of fruiting in t between reproduction in t and reproduction in t+1 for adults
ggplot(data = wrangled_fem_to_df, aes(x = Fruit_t0, y = Growth_t0_t1)) +
  geom_jitter(width = 0.2, height = 0.1, alpha = 0.3) +
  stat_smooth(method = "lm", se = TRUE) +
  labs(x = "Fruiting probability (t)", y = "RGR (t to t+1)") +
  theme_bw()

# Model trade-off between reproduction in t and growth from t to t+1 for adults
gro_fru_to_mem <- lmer(Growth_t0_t1 ~ Fruit_t0 + (1 | Indiv_ID), 
                       data = wrangled_fem_to_df)
summary(gro_fru_to_mem)
# Positive relationship between reproduction and fruiting - singular fit, cannot interpret

# Input model data into data frame
to_mods <- c(sur_gro_to_glmm,
             sur_fru_to_glmm,
             fru_gro_to_glmm,
             gro_fru_to_mem)

for(i in 1:length(to_mods)){
  summ <- summary(to_mods[[i]])
  vr_to_df[i,"VR_1"] <- as.character(summ$call$formula[[3]][[2]])
  vr_to_df[i,"VR_2"] <- as.character(summ$call$formula[[2]])
  vr_to_df[i,"Model"] <- as.character(to_mods[[i]]@call)[2]
  vr_to_df[i,"Est_p"] <- paste("β=", round(summ$coefficients[2,1], digits = 3), 
                               " (P=", signif(summ$coefficients[2,ncol(summ$coefficients)], digits = 3),")", sep = "")
  vr_to_df[i,"N_obs"] <- as.numeric(summ$devcomp$dims[1])
  vr_to_df[i,"N_ind"] <- as.numeric(summ$ngrps)
}

# # Export model summary as .csv - still do this?
# write.csv(vr_to_df, "Table S2.csv")

# Recruitment -------------------------------------------------------------
# Mistletoes are recruited through spread of seed from fruiting individuals to a suitable branch
# Seeds germinate and establish for 3 years before they can be observed (in their 3rd year)
# I.e., it takes three time-steps for fruiting individuals to contribute to the juvenile pool
# As we lack data on germination and establishment success, we estimate an initial value of s0 linking fruiting in t to juvenile recruitment in t+3
# Originally, individuals were considered to produce berries directly proportional to the number of terminal shoots, which is directly proportional to their size - b = k * e^(log(a))
# According to Mellado and Zamora, 2014, 1m^2 (10000cm^2) of mistletoe produces ~2000 berries, suggesting that 2000 = k * exp(ln(10000)) => k = 0.2
k <- 0.2

# Predict number of berries produced
wrangled_by_yr_df$Berries <- k * exp(wrangled_by_yr_df$logArea)

# Plot berry production
Fig_S5 <- ggplot(data = wrangled_by_yr_df, aes(x = logArea, y = Berries)) +
  geom_line(col = colours(6)[3]) +
  theme_bw() +
  labs(x = xl, y = "Number of berries") +
  theme(axis.title = element_text(size = 15))
Fig_S5

# As per Lucas et al., 2008, Fertility from t to t + 3 = (total no. offspring in t+3)/(total no. berries in t) * (no. berries for individual in t) = Σrec/Σb * b(z)
# Via two intermediary stages: 1-year-old seedlings (S1) and 2-year-old seedlings (S2)
# Σrec(t+3)/Σb(t) = S1(t+1)/b(t) * S2(t+2)/S1(t+1) * rec(t+3)/S2(t+2) = s0 * s1 * s2
# Assuming constant survival as we cannot measure establishment s0, survival of 1-year-old seedlings s1 or survival of 2-year-old seedlings s2, s0 = s1 = s2
# therefore s0 = (Σrec/Σb)^(1/3)
# Σb is estimated using the berry production constant k, where b = k*exp(z)
# s0 will later be estimated more realistically through iteration to obtain λ=1.1

# When the IPM is built using a particular berry production constant k, s0 will be estimated as follows
get.s0 <- function(k, s1, s2){
  # Filter long df for only fruiting individuals and estimate berry production
  fru_df <- filter(wrangled_by_yr_df, Fruit == 1) %>%
    mutate(Berry_no = k * exp(logArea))
  # mutate(Berry_no = k/(1+exp((-1)*(logArea-0.5*(max_size + min_size_rep)))))
  # Calculate total number of berries each year
  berries_df <- fru_df %>%
    group_by(Year) %>%
    summarise(Total_berries = sum(Berry_no, na.rm = TRUE))
  
  # Estimate number of recruits in t+3
  rec_df <- filter(wrangled_by_yr_df, 
                   Status == "FC" & Stage == "Juvenile") %>%
    count(Year)
  
  # Run through each year and estimate fertility from t to t+3: Σrec(t+3)/Σb(t)
  f_t_t3 <- list()
  for(year in 14:20){
    b_t    <- berries_df[berries_df$Year == year, ]$Total_berries
    rec_t3 <- rec_df[rec_df$Year == year+3, ]$n
    f_t_t3[year - 13] <- rec_t3/b_t
  }
  # Estimate mean fertility across all 3-year periods
  mean_f_t_t3 <- mean(unlist(f_t_t3))
  
  # Estimate s0
  s0 <- mean_f_t_t3 / (s1 * s2)
  
  # Return s0
  return(s0)
}

# Find distribution of new recruit (3-year-old) sizes
# Filter only new recruits
rec_only_df <- filter(wrangled_by_yr_df, 
                      Status == "FC" & Stage == "Juvenile")

# Plot distribution of recruit sizes
ggplot(data = rec_only_df, aes(x = logArea)) +
  geom_histogram(aes(y = ..density..), alpha = 0.3, fill = colours(6)[3], bins = 10) +
  stat_function(fun = function(x) {
    dnorm(x, mean = mean(rec_only_df$logArea), sd = sd(rec_only_df$logArea))
  }, color = colours(6)[3]) +
  labs(x = xl, y = "Density") + 
  theme_bw() +
  theme(text = element_text(size = 14),           
        axis.title = element_text(size = 16),    
        axis.text = element_text(size = 14))

# Plot distribution of heights for all individuals to get a better height distribution
Fig_S6b <- ggplot(data = wrangled_by_yr_df, aes(x = Height)) +
  geom_histogram(aes(y = ..density..), alpha = 0.3, fill = colours(6)[3], bins = 10) +
  stat_function(fun = function(x) {
    dnorm(x, mean = mean(wrangled_by_yr_df$Height), sd = sd(wrangled_by_yr_df$Height))
  }, color = colours(6)[3]) +
  labs(x = "Height (m)", y = "") + 
  theme_bw() 



# Build IPM ---------------------------------------------------------------
## Extract parameters
# Extract parameters from selected models
# Cannot used random effect because there is no way of tracking individual ID in an IPM
# Survival ~ log(Area)
sur_area_glm <- glm(Survive ~ logArea_t0, 
                    data = wrangled_df, family = binomial)
summary(sur_area_glm)

# Plot Survival ~ log(Area)
ggplot(data = wrangled_df, aes(x = logArea_t0, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm",
              method.args = list(family = binomial)) +
  theme_bw()

# Growth ~ log(Area) - quantify as absolute growth rather than relative rate to project next size
gro_abs_area_quad_lm <- lm(logArea_t1 ~ logArea_t0 + I(logArea_t0^2), 
                           data = wrangled_df)
summary(gro_abs_area_quad_lm)

# Plot log(Area) in t+1 ~ log(Area) in t
ggplot(data = wrangled_df, aes(x = logArea_t0, y = logArea_t1)) +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  theme_bw()

# Find mean 3-yo juvenile size in time t
S_3yo <- mean(rec_only_df$logArea)

# Seedlings below 3 years cannot be measured, so their vital rates are extrapolated to complete the life cycle
# Extrapolate to estimate 2-yo juvenile size at t-1, by taking inverse of quadratic function
# Define the quadratic function parameters
a <- coef(gro_abs_area_quad_lm)[3]
b <- coef(gro_abs_area_quad_lm)[2]
c <- coef(gro_abs_area_quad_lm)[1]

# Function to compute inverse of quadratic
quad_inverse <- function(S_next, a, b, c) {
  discriminant <- b^2 - 4 * a * (c - S_next)
  if (discriminant < 0) {
    return(NA)  # No real solution
  }
  S_prev1 <- (-b + sqrt(discriminant)) / (2 * a)
  S_prev2 <- (-b - sqrt(discriminant)) / (2 * a)
  return(c(S_prev1, S_prev2))  # Returns both possible solutions
}

# Estimate 2-yo juvenile size at t-1
S_2yo <- quad_inverse(S_3yo, a, b, c)
S_2yo <- max(S_2yo)

# Estimate 1-yo juvenile size at t-2
S_1yo <- quad_inverse(S_2yo, a, b, c)
S_1yo <- max(S_1yo)

# Plot distribution of 1yo seedling sizes
S_1yo_sd <- sd(rec_only_df$logArea)
rec_distrib <- data.frame(x = seq(S_1yo - 4 * S_1yo_sd, S_1yo + 4 * S_1yo_sd, length.out = 100))
Fig_S6a <- ggplot(data = rec_distrib, aes(x = x)) +
  stat_function(fun = dnorm, args = list(mean = S_1yo, sd = S_1yo_sd), color = colours(6)[3], size = 1) +
  labs(x = xl, y = "Density") +
  theme_bw()

# Plot Figure S6
Fig_S6 <- (Fig_S6a + Fig_S6b) + plot_annotation(tag_levels = "a")
Fig_S6

# Estimate survival for 2yo and 1yo juveniles
data.frame(logArea_t0 = c(S_2yo, S_1yo))
s2 <- exp(as.numeric(coef(sur_area_glm)[1]) + as.numeric(coef(sur_area_glm)[2]) * S_2yo) / (1 + exp((as.numeric(coef(sur_area_glm)[1]) + as.numeric(coef(sur_area_glm)[2]) * S_2yo)))
s1 <- exp(as.numeric(coef(sur_area_glm)[1]) + as.numeric(coef(sur_area_glm)[2]) * S_1yo) / (1 + exp((as.numeric(coef(sur_area_glm)[1]) + as.numeric(coef(sur_area_glm)[2]) * S_1yo)))

# Chosen k suggests s0 ~ 0.05, but survival difficult to estimate, so s0 will be fixed such that λ=1.1
get.s0(k, s1, s2)

# Adult fruiting ~ log(Area) + Height
fru_area_ht_glm <- glm(Fruit ~ logArea + Height, 
                       data = wrangled_by_yr_df, family = binomial)
summary(fru_area_ht_glm)

# Set up model parameters from fitted models
params <- data.frame(
  sur.int    = summary(sur_area_glm)$coefficients[1,1],
  sur.slope  = summary(sur_area_glm)$coefficients[2,1],
  gro.int    = summary(gro_abs_area_quad_lm)$coefficients[1,1],
  gro.slope  = summary(gro_abs_area_quad_lm)$coefficients[2,1],
  gro.slope2 = summary(gro_abs_area_quad_lm)$coefficients[3,1],
  gro.sd     = sd(resid(gro_abs_area_quad_lm)),
  fru.int        = summary(fru_area_ht_glm)$coefficients[1,1],
  fru.area.slope = summary(fru_area_ht_glm)$coefficients[2,1],
  fru.ht.slope   = summary(fru_area_ht_glm)$coefficients[3,1],
  k.berries      = k, 
  s0             = get.s0(k, s1, s2),  
  rec.area.mean  = S_1yo,
  rec.area.sd    = S_1yo_sd,
  rec.ht.mean    = mean(wrangled_by_yr_df$Height),
  rec.ht.sd      = sd(wrangled_by_yr_df$Height)
)

#Create function to build IPM and plot it at a reference height
build.ipm <- function(params, mesh, ref_ht){
  
  # Create outputs list
  output_ls <- list()
  
  # Use model coefficients to build vital rate functions
  # Keep mistletoe height constant unless new height generated for a new recruit
  ht_fun <- function (h1, h) {
    if (h1 == h) {return(1)}
    if (h1 != h) {return(0)}
  }
  
  ## Define vital rate functions
  # Survival function
  sur_fun <- function(z){
    u = exp(params$sur.int + params$sur.slope * z)
    return(u/(1+u))
  }
  
  # Growth function
  gro_fun <- function(z1, z) {
    dnorm(z1, mean = params$gro.int + params$gro.slope * z + params$gro.slope2 * z * z, 
          sd = params$gro.sd)
  }
  
  # Fruiting function
  fru_fun <- function(z, h) {
    u = exp(params$fru.int + params$fru.area.slope * z + params$fru.ht.slope * h)
    return(u/(1+u))
  }
  
  # Fecundity (berries produced per fruiting individual)
  fec_fun <- function(z) {
    params$k.berries * exp(z)
    # params$k.berries / ((1+exp((-1)*(z-0.5*(max_size + min_size_rep)))))
  }
  
  # 1-year-old seedling heights
  S1_ht_fun <- function(h1){
    dnorm(h1, mean = params$rec.ht.mean, sd = params$rec.ht.sd)
  }
  
  # Recruit sizes
  R_z_fun <-function(z1){
    dnorm(z1, mean = params$rec.area.mean, sd = params$rec.area.sd)
  }
  
  ## Construct sub-kernels
  ### P sub-kernel
  # Describes probability of an individual of size z and height h becoming a individual of size z1 and height h1
  
  # P sub-kernel - sur * gro 
  P_zh <- function(z1, h1, z, h) {
    return(sur_fun(z) * gro_fun(z1, z) * ht_fun(h1, h))
  }
  
  # Number of meshpoints
  meshz <- mesh
  meshh <- mesh
  
  # Set limits for integration - may impact eviction
  Lz <- 0.95 * (params$rec.area.mean - 4 * params$rec.area.sd)
  Uz <- 1.05 * max(wrangled_by_yr_df$logArea, na.rm = TRUE)
  Lh <- 0.95 * min(wrangled_by_yr_df$Height, na.rm = TRUE)
  Uh <- 1.05 * max(wrangled_by_yr_df$Height, na.rm = TRUE)
  
  # Set bin size and midpoints for each trait variable
  
  # Log area bin size
  hz <- (Uz - Lz) / meshz
  
  # Log area midpoints
  yz <- Lz + hz * ((1:meshz) - 0.5)
  
  # Height bin size
  hh <- (Uh - Lh) / meshh
  
  # Height midpoints
  yh <- Lh + hh * ((1:meshh) - 0.5)
  
  # Create function to fun through mesh points
  eta_ij <- function(i, j, meshz) {(j - 1) * meshz + i}
  
  # Create matrix of indices with which to evaluate mega matrix
  Eta <- outer(1:meshz, 1:meshh, eta_ij, meshz = meshz); 
  
  # Create mega matrix of 0s in which to put transition probabilities at each area-height combination
  P <- matrix(0, meshz*meshh, meshz*meshh)
  
  # Create array of 0s in which to transition probabilities at each area-height combination
  Pvals <- array(0, c(meshz, meshh, meshz, meshh))
  
  # Run through each set of meshpoints and evaluate transition probabilities
  for(i in 1:meshz){
    for(j in 1:meshh){
      for(k in 1:meshh){
        pvals = P_zh(yz, yh[k], yz[i], yh[j])
        P[Eta[,k], Eta[i,j]] = pvals
        Pvals[,k,i,j] = pvals
      }
    }
  }
  
  # Multiply matrix by bin sizes to standardise transition probabilities
  P <- hz * hh * P
  
  # Plot P kernel at mean height and area
  image(t(log(P[(meshz*meshh/2+1):(meshz*meshh/2+meshz), (meshz*meshh/2+1):(meshz*meshh/2+meshz)] + 0.01)))
  
  ### F sub-kernel
  # Describes number of S1 from an individual of size z and height h establishing at height h1
  # CS1 sub-kernel - fru prob * fec * estab prob * new ht
  F_zh <- function(z1, h1, z, h) {
    return(fru_fun(z, h) * fec_fun(z) * params$s0 * S1_ht_fun(h1) * R_z_fun(z1))
  }
  
  # Create function to fun through mesh points
  eta_ij <- function(i, j, meshz) {(j - 1) * meshz + i}
  
  # Create matrix of indices with which to evaluate mega matrix
  Eta <- outer(1:meshz, 1:meshh, eta_ij, meshz = meshz); 
  
  # Modify kernel sizes as z1 is discrete
  # Create mega matrix of 0s in which to put transition probabilities at each area-height combination
  F <- matrix(0, meshz * meshh, meshz * meshh)
  
  # Create array of 0s in which to transition probabilities at each area-height combination
  Fvals <- array(0, c(meshz, meshh, meshz, meshh))
  
  # Run through each set of meshpoints and evaluate transition probabilities
  for(i in 1:meshz){
    for(j in 1:meshh){
      for(k in 1:meshh){
        fvals = F_zh(yz, yh[k], yz[i], yh[j])
        F[Eta[,k], Eta[i,j]] = fvals
        Fvals[,k,i,j] = fvals
      }
    }
  }
  
  # Multiply matrix by bin sizes to standardise transition probabilities
  F <- hz * hh * F
  
  # Plot F kernel at mean height and area
  image(t(log(F[(meshz*meshh/2+1):(meshz*meshh/2+meshz), (meshz*meshh/2+1):(meshz*meshh/2+meshz)] + 0.01)))
  
  # Find bins of min size at reproduction and set non-reproductive bins to 0
  min_size_rep_bin <- ceiling((min_size_rep - Lz) / (Uz - Lz) * mesh)
  
  # Non-reproductive bins
  non_rep_bins <- c(sapply(seq(1, meshz*meshh, by = mesh), function(x) seq(x, x + (min_size_rep_bin - 1))))
  
  # Designate non-reproductive bins as 0
  for(bin in non_rep_bins){
    F[, bin] <- 0
  }
  
  # Combine P and F into mega matrix kernel K
  K <- P + F
  
  # Output P kernel
  output_ls[[1]] <- P
  # Output F kernel
  output_ls[[2]] <- F
  # Output full kernel
  output_ls[[3]] <- K
  # Output meshes
  output_ls[[4]] <- mesh
  
  # Rename outputs
  names(output_ls) <- c("P", "F", "K", "mesh")
  
  # Return outputs
  return(output_ls)
  
  # Find bin of reference height
  ref_ht_bin <- as.integer((ref_ht - Lh)/(Uh - Lh) * meshh)
  
  # Plot K kernel at reference height
  image(t(log(K[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)] + 0.01)))
  
  # Create P kernel at reference height
  ref_P <- P[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)]
  
  # Create F kernel at reference height
  ref_F <- F[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)]
  
  # Plot full kernel
  ref_K <- ref_P + ref_F
  image(t(ref_K))

}

# Build initial IPM with 50 meshpoints per state variable and reference height as mean height
mean_ht <- params$rec.ht.mean
outputs <- build.ipm(params = params, mesh = 50, ref_ht = mean_ht)
lambda(outputs[["K"]])

# Extract lambda and life history traits from IPM
get.traits <- function(outputs){

  # Calculate lambda
  outputs[[5]] <- eigen.analysis(outputs[["K"]])$lambda1
  names(outputs)[5] <- "Lambda"
  
  # Calculate R0 from P and F sub-kernels
  outputs[[6]] <- net_repro_rate(outputs[["P"]], outputs[["F"]])
  names(outputs)[6] <- "R0"
  
  # Calculate generation time from R0 and lambda
  outputs[[7]] <- abs(log(outputs[["R0"]])/log(outputs[["Lambda"]]))
  names(outputs)[7] <- "GenTfun"
  
  # Calculate mean life expectancy from mature distribution
  # Identify which stages are reproductive
  mat_dist <- mature_distrib(outputs[["P"]], repro_stages = repro_stages(outputs[["F"]]))
  outputs[[8]] <- life_expect_mean(matU = outputs[["P"]], mixdist = mat_dist, start = NULL)
  names(outputs)[8] <- "Lmean"
  
  # Calculate mean age at maturity, from first non-propagule stage
  outputs[[9]] <- mature_age(matU = outputs[["P"]], matF = outputs[["F"]])
  names(outputs)[9] <- "La"
  
  # Calculate mean reproductive window
  outputs[[10]] <- outputs[["Lmean"]] - outputs[["La"]]
  names(outputs)[10] <- "Lamean"
  
  # Return outputs
  return(outputs)
}

# Extract initial life history traits
traits <- get.traits(outputs)

# Fix establishment constant s0 -------------------------------------------
# Create lists to store lambdas and s0 values
lambdas <- list()
s0s <- list()

# Begin with initial lambda
lambdas[[1]] <- lambda(outputs[["K"]])

# Start with s0 derived from k = 0.2
s0s[[1]] <- get.s0(params$k.berries, s1, s2)

# Initialise loop, set number of iterations and set progress bar
i <- 1
iter <- 100
pb <- txtProgressBar(min = 0, max = iter, style = 3)

# If lambda > 1.11, decrease s0; if lambda < 1.09, increase s0 
while(abs(1.1 - as.numeric(lambdas[[i]])) > 0.01 & i < iter) {
  if(lambdas[[i]] > 1.1){
    params$s0 <- 0.9 * s0s[[i]]
  } else {
    params$s0 <- 1.1 * s0s[[i]]
  }
  s0s[[i+1]] <- params$s0
  outputs <- build.ipm(params = params, mesh = 50, ref_ht = params$rec.ht.mean)
  names(outputs) <- c("P", "F", "K", "mesh")
  lambdas[[i+1]] <- lambda(outputs[["K"]])
  paste("s0 = ", s0s[[i+1]], "; Lambda = ", lambdas[[i+1]])
  setTxtProgressBar(pb, i)
  i <- i + 1
}
close(pb)

# Set s0 as iteration output
params$s0 <- s0s[[length(s0s)]]

# Build final IPM ---------------------------------------------------------
# Final params
final_params <- params

# Build and plot final IPM with 50 meshpoints
outputs <- build.ipm(params = final_params, mesh = 50, ref_ht = final_params$rec.ht.mean)

# Set maxima for plotting
max_P <- max(outputs[["P"]])
max_F <- max(outputs[["F"]])

# Extract final life history traits
final_traits <- get.traits(outputs)

# Export final parameters for mistletoe LHTs
write.csv(final_params, "IPM_parameters.csv")

# Create new x and y labels
xl_t <- expression(log ~ "area" ~ cm^2 ~ "in t")
yl_t1 <- expression(log ~ "area" ~ cm^2 ~ "in t+1")

# Plot final IPM at mean height 
plot.ipm <- function(P_kernel, F_kernel, ref_ht, mesh){
  # Set limits to specify bins and for plotting
  Lz <- 0.95 * (params$rec.area.mean - 4 * params$rec.area.sd)
  Uz <- 1.05 * max(wrangled_by_yr_df$logArea, na.rm = TRUE)
  Lh <- 0.95 * min(wrangled_by_yr_df$Height, na.rm = TRUE)
  Uh <- 1.05 * max(wrangled_by_yr_df$Height, na.rm = TRUE)
  
  # Find bin of reference height
  ref_ht_bin <- as.integer((ref_ht - Lh)/(Uh - Lh) * mesh)
  
  # Create P kernel at reference height
  ref_P <- P_kernel[((ref_ht_bin-1)*mesh + 1):(ref_ht_bin*mesh), ((ref_ht_bin-1)*mesh + 1):(ref_ht_bin*mesh)]
  
  # Create F kernel at reference height
  ref_F <- F_kernel[((ref_ht_bin-1)*mesh + 1):(ref_ht_bin*mesh), ((ref_ht_bin-1)*mesh + 1):(ref_ht_bin*mesh)]
  
  # Divide by maximum value for plotting
  plot_P <- t(ref_P/max_P)
  plot_F <- t(ref_F/max_F)
  
  # Convert matrices into long format
  dfP <- melt(plot_P)
  dfF <- melt(plot_F)
  
  # Custom axis labels (replace original scale with these)
  custom_labels <- signif(seq(Lz, Uz, length.out = mesh), 2)
  
  # Create the plot
  main_plot <- ggplot() +
    # First matrix with gradient
    geom_tile(data = dfP, aes(x = Var1, y = Var2, fill = value), alpha = max(ref_P)/max_P) +
    scale_fill_gradient(name = "Transition probability", low = "white", high = colours(6)[6], guide = "none") +
    
    # Second matrix with a different gradient
    geom_tile(data = dfF, aes(x = Var1, y = Var2, alpha = value), fill = colours(6)[3]) +
    scale_alpha_continuous(name = "Fecundity (scaled)", range = c(0, max(ref_F)/max_F), 
                           guide = "none") +
    
    scale_x_continuous(name = xl_t, expand = c(0, 0), 
                       breaks = seq(1, mesh, by = 10), 
                       labels = custom_labels[seq(1, 50, by = 10)]) +
    scale_y_continuous(name = yl_t1, expand = c(0, 0), 
                       breaks = seq(1, mesh, by = 10), 
                       labels = custom_labels[seq(1, 50, by = 10)]) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          # Ensure axis labels are closer to the graph
          axis.title.x = element_text(margin = margin(t = 5), size = 20),  # Moves x-axis label closer
          axis.title.y = element_text(margin = margin(r = 5), size = 20)  # Moves y-axis label closer
          )

  #  Fake Legend Plot for P
  P_legend_plot <- ggplot(dfP, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(name = "Transition probability", low = "white", high = colours(6)[6], limits = c(0, max_P)) + 
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(size = 15))  # Keep only legend
  
  #  Fake Legend Plot for F (Uses Fill to Simulate Alpha) 
  F_legend_plot <- ggplot(dfF, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient(name = "Fecundity", low = "white", high = colours(6)[3], limits = c(0, max_F)) +  # Mimic alpha
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(size = 15))  # Keep only legend
  
  # Extract the Legends 
  P_legend <- get_legend(P_legend_plot)
  F_legend <- get_legend(F_legend_plot)
  
  # Combine Main Plot and Separate Legend 
  final_plot <- plot_grid(main_plot, P_legend, F_legend, ncol = 3, rel_widths = c(4, 1, 0.5))
  return(final_plot)
  }

Fig_3b <- plot.ipm(P_kernel = outputs[["P"]], F_kernel = outputs[["F"]], ref_ht = mean_ht, mesh = outputs[["mesh"]])

Fig_3b

# Extract output LHTs -----------------------------------------------------
# Create data frame for mistletoe LHTs
Mst_LHTs <- data.frame(SpeciesAccepted = "Viscum album",
                       Lambda = traits[["Lambda"]],
                       R0 = traits[["R0"]],
                       GenTfun = traits[["GenTfun"]],
                       Lmean = traits[["Lmean"]],
                       La = traits[["La"]],
                       Lamean = traits[["Lamean"]])

# Export .csv for mistletoe LHTs
write.csv(Mst_LHTs, "Mst_LHTs.csv")

# Sensitivity of traits to parameters (including k) ---------------------------------
traits_sens_params <- list()
traits_sens_params[[1]] <- final_traits
traits_sens_params[[1]]$P <- NA
traits_sens_params[[1]]$F <- NA
traits_sens_params[[1]]$K <- NA

for(i in 2:(ncol(params)+1)){
  paste("Parameter = ", names(params)[i-1])
  params[,i-1] <- final_params[,i-1] + 0.001
  outputs <- build.ipm(params = params, mesh = 50, ref_ht = params$rec.mean.ht)
  traits_sens_params[[i]] <- get.traits(outputs)
  traits_sens_params[[i]]$P <- NA
  traits_sens_params[[i]]$F <- NA
  traits_sens_params[[i]]$K <- NA
  params[,i-1] <- params[,i-1] - 0.001
}

# Create data frame in which to store raw values
traits_sens_params_df <- data.frame(Parameter = c("none", names(params)), 
                                    Lambda = NA,
                                    R0 = NA,
                                    GenTfun = NA,
                                    Lmean = NA,
                                    La = NA, 
                                    Lamean = NA)
# Input raw outputs
for(i in 1:length(traits_sens_params)){
  traits_sens_params_df[i,]$Lambda <- traits_sens_params[[i]]$Lambda
  traits_sens_params_df[i,]$R0 <- traits_sens_params[[i]]$R0
  traits_sens_params_df[i,]$GenTfun <- traits_sens_params[[i]]$GenTfun
  traits_sens_params_df[i,]$Lmean <- traits_sens_params[[i]]$Lmean
  traits_sens_params_df[i,]$La <- traits_sens_params[[i]]$La
  traits_sens_params_df[i,]$Lamean <- traits_sens_params[[i]]$Lamean
}

# Export raw trait values
write.csv(traits_sens_params_df, "Raw parameter sensitivity values")

# Create data frame in which to store sensitivities 
traits_sens_df <- data.frame(Parameter = c("none", names(params)), 
                                    Lambda = NA,
                                    R0 = NA,
                                    GenTfun = NA,
                                    Lmean = NA,
                                    La = NA, 
                                    Lamean = NA)

# Calculate sensitivities
for(i in 1:nrow(traits_sens_params_df)){
  for(j in 2:ncol(traits_sens_params_df)){
    traits_sens_df[i,j] <- (traits_sens_params_df[1,j] - traits_sens_params_df[i,j])/0.001
  }
}

# Export trait sensitivity
write.csv(traits_sens_df, "Parameter sensitivities")

# Plot sensitivities
# Convert to long format
traits_sens_df_long <- traits_sens_df %>%
  pivot_longer(cols = -Parameter,
               names_to = "Output",
               values_to = "Sensitivity")
# Define custom colors for parameters
param_colors <- c("sur.int"   = colours(6)[1],
                  "sur.slope"  = colours(6)[1],
                  "gro.int"    = colours(6)[6],
                  "gro.slope"  = colours(6)[6],
                  "gro.slope2" = colours(6)[6],
                  "gro.sd"     = colours(6)[6],
                  "fru.int"        = colours(6)[3],
                  "fru.area.slope" = colours(6)[3],
                  "fru.ht.slope"   = colours(6)[3],
                  "k.berries"      = colours(6)[3], 
                  "s0"             = colours(6)[3],  
                  "rec.area.mean"  = colours(6)[3],
                  "rec.area.sd"    = colours(6)[3],
                  "rec.ht.mean"    = colours(6)[3],
                  "rec.ht.sd"      = colours(6)[3])

# Give parameters labels and order
param_labels <- c("sur.int"   = "βs0",
                  "sur.slope"  = "βs1",
                  "gro.int"    = "βg0",
                  "gro.slope"  = "βg1",
                  "gro.slope2" = "βg2",
                  "gro.sd"     = "σg",
                  "fru.int"        = "βf0",
                  "fru.area.slope" = "βf1",
                  "fru.ht.slope"   = "βf2",
                  "k.berries"      = "k", 
                  "s0"             = "s0",  
                  "rec.area.mean"  = "μz",
                  "rec.area.sd"    = "σz",
                  "rec.ht.mean"    = "μh",
                  "rec.ht.sd"      = "σh")
param_order <- c("βs0",
                 "βs1",
                  "βg0",
                  "βg1",
                  "βg2",
                  "σg",
                  "βf0",
                  "βf1",
                  "βf2",
                  "k", 
                  "s0",  
                  "μz",
                  "σz",
                  "μh",
                  "σh")


# Give output labels and order
output_labels <- c("Lambda" = "λ", 
                   "R0" = "R0", 
                   "GenTfun" = "Generation time", 
                   "Lmean" = "Mean life expectancy",
                   "La" = "Mean age at maturity",
                   "Lamean" = "Reproductive window")
output_order <- c("Lambda", 
                  "R0", 
                  "GenTfun", 
                  "Lmean",
                  "La",
                  "Lamean")

traits_sens_df_long <- traits_sens_df_long %>%
  mutate(
    Parameter_label = recode(Parameter, !!!param_labels),
    FillGroup = Parameter,
    Parameter_label = factor(Parameter_label, levels = param_order),  # keep original names for color mapping
    Output = factor(Output, levels = output_order, labels = output_labels)
  )

# Plot
Fig_3c <- ggplot(filter(traits_sens_df_long, Parameter != "none"), aes(x = Parameter_label, y = Sensitivity, fill = FillGroup)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Output, scales = "free_y") +
  theme_minimal() +
  labs(x = "Parameter", y = "Sensitivity of trait to parameter") +
  theme(
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # Optional: hide legend if bars are labeled
  ) +
  scale_fill_manual(values = param_colors)
Fig_3c

# Sensitivity of traits to mesh points ----------------------------------------------
# Set mesh sizes
meshes <- seq(30, 50, by = 5)

# Create list to store lambda values
traits_sens_mesh <- list()

# Extract traits for each mesh value
for(mesh in meshes){
  outputs <- build.ipm(params = final_params, mesh = mesh, ref_ht = final_params$rec.ht.mean)
  traits_sens_mesh[[mesh/5 - 5]] <- get.traits(outputs)
}

# Wrangle into dataframe
traits_sens_mesh_df <- data.frame(mesh = meshes, 
                             Lambda = NA,
                             R0 = NA,
                             GenTfun = NA,
                             Lmean = NA,
                             La = NA, 
                             Lamean = NA)
for(mesh in meshes){
    traits_sens_mesh_df[which(traits_sens_mesh_df$mesh == mesh),]$Lambda <- traits_sens_mesh[[mesh/5 - 5]]$Lambda
    traits_sens_mesh_df[which(traits_sens_mesh_df$mesh == mesh),]$R0 <- traits_sens_mesh[[mesh/5 - 5]]$R0
    traits_sens_mesh_df[which(traits_sens_mesh_df$mesh == mesh),]$GenTfun <- traits_sens_mesh[[mesh/5 - 5]]$GenTfun
    traits_sens_mesh_df[which(traits_sens_mesh_df$mesh == mesh),]$Lmean <- traits_sens_mesh[[mesh/5 - 5]]$Lmean
    traits_sens_mesh_df[which(traits_sens_mesh_df$mesh == mesh),]$La <- traits_sens_mesh[[mesh/5 - 5]]$La
    traits_sens_mesh_df[which(traits_sens_mesh_df$mesh == mesh),]$Lamean <- traits_sens_mesh[[mesh/5 - 5]]$Lamean
}

# Export sensitivities to meshpoints
write.csv(traits_sens_mesh_df, "Mesh sensitivities.csv")

# Plot traits against number of mesh points
Fig_S8a <- ggplot(data = traits_sens_mesh_df, aes(x = (mesh*mesh), y = Lambda)) +
  labs(x = "Total meshpoints", y = "λ") +
  geom_line() +
  theme_bw()
Fig_S8a
Fig_S8b <- ggplot(data = traits_sens_mesh_df, aes(x = (mesh*mesh), y = R0)) +
  labs(x = "Total meshpoints", y = expression("R"["0"])) +
  geom_line() +
  theme_bw()
Fig_S8b
Fig_S8c <- ggplot(data = traits_sens_mesh_df, aes(x = (mesh*mesh), y = GenTfun)) +
  labs(x = "Total meshpoints", y = "T") +
  geom_line() +
  theme_bw()
Fig_S8c
Fig_S8d <- ggplot(data = traits_sens_mesh_df, aes(x = (mesh*mesh), y = Lmean)) +
  labs(x = "Total meshpoints", y = expression(eta["e"])) +
  geom_line() +
  theme_bw()
Fig_S8d
Fig_S8e <- ggplot(data = traits_sens_mesh_df, aes(x = (mesh*mesh), y = La)) +
  labs(x = "Total meshpoints", y = expression("L"[alpha])) +
  geom_line() +
  theme_bw()
Fig_S8e
Fig_S8f <- ggplot(data = traits_sens_mesh_df, aes(x = (mesh*mesh), y = Lamean)) +
  labs(x = "Total meshpoints", y = expression(paste("L"[paste(alpha,"-",omega)]))) +
  geom_line() +
  theme_bw()
Fig_S8f

Fig_S8 <- (Fig_S8a + Fig_S8b) / (Fig_S8c + Fig_S8d) / (Fig_S8e + Fig_S8f)
Fig_S8

# Sensitivity of choice of lambda -------------------------------------------
# Make data frame for sensitivity of traits to lambda choice
sens_lambda_df <- data.frame(lambda = c(1.1, 1.05, 1.15), 
                             s0 = signif(c(params$s0, NA, NA), 3),
                             GenTfun = NA,
                             Lmean = NA,
                             La = NA,
                             Lamean = NA)

# Input values from final model
outputs <- build.ipm(params = final_params, mesh = 50, ref_ht = final_params$rec.ht.mean)
traits <- get.traits(outputs)
sens_lambda_df[which(sens_lambda_df$lambda == 1.1), ]$GenTfun <- traits$GenTfun
sens_lambda_df[which(sens_lambda_df$lambda == 1.1), ]$Lmean <- traits$Lmean
sens_lambda_df[which(sens_lambda_df$lambda == 1.1), ]$La <- traits$La
sens_lambda_df[which(sens_lambda_df$lambda == 1.1), ]$Lamean <- traits$Lamean

# Create lists to store lambdas and s0 values
lambdas <- list()
s0s <- list()

# Begin with initial lambda
lambdas[[1]] <- traits$Lambda

# Start with s0 derived from k = 0.2
s0s[[1]] <- final_params$s0

# Initialise loop, set number of iterations and set progress bar
i <- 1
iter <- 10
pb <- txtProgressBar(min = 0, max = iter, style = 3)

# If lambda > 1.06, decrease s0; if lambda < 1.04, increase s0 
while(abs(1.05 - as.numeric(lambdas[[i]])) > 0.01 & i < iter) {
  if(lambdas[[i]] > 1.05){
    params$s0 <- 0.9 * s0s[[i]]
  } else {
    params$s0 <- 1.1 * s0s[[i]]
  }
  s0s[[i+1]] <- params$s0
  outputs <- build.ipm(params = params, mesh = 50, ref_ht = params$rec.ht.mean)
  names(outputs) <- c("P", "F", "K", "mesh")
  lambdas[[i+1]] <- lambda(outputs[["K"]])
  paste("s0 = ", s0s[[i+1]], "; Lambda = ", lambdas[[i+1]])
  setTxtProgressBar(pb, i)
  i <- i + 1
}
close(pb)

# Set s0 as iteration output
params$s0 <- s0s[[length(s0s)]]

# Get traits for IPM
traits <- get.traits(build.ipm(params = params, mesh = 50, ref_ht = params$rec.ht.mean))
sens_lambda_df[which(sens_lambda_df$lambda == 1.05),]$s0 <- params$s0
sens_lambda_df[which(sens_lambda_df$lambda == 1.05),]$GenTfun <- traits$GenTfun
sens_lambda_df[which(sens_lambda_df$lambda == 1.05), ]$Lmean <- traits$Lmean
sens_lambda_df[which(sens_lambda_df$lambda == 1.05), ]$La <- traits$La
sens_lambda_df[which(sens_lambda_df$lambda == 1.05), ]$Lamean <- traits$Lamean

# Repeat for λ=1.15
# Input values from final model
outputs <- build.ipm(params = final_params, mesh = 50, ref_ht = final_params$rec.ht.mean)
traits <- get.traits(outputs)
# Create lists to store lambdas and s0 values
lambdas <- list()
s0s <- list()

# Begin with initial lambda
lambdas[[1]] <- traits$Lambda

# Start with s0 derived from k = 0.2
s0s[[1]] <- final_params$s0

# Initialise loop, set number of iterations and set progress bar
i <- 1
iter <- 10
pb <- txtProgressBar(min = 0, max = iter, style = 3)

# If lambda > 1.16, decrease s0; if lambda < 1.14, increase s0 
while(abs(1.15 - as.numeric(lambdas[[i]])) > 0.01 & i < iter) {
  if(lambdas[[i]] > 1.15){
    params$s0 <- 0.9 * s0s[[i]]
  } else {
    params$s0 <- 1.1 * s0s[[i]]
  }
  s0s[[i+1]] <- params$s0
  outputs <- build.ipm(params = params, mesh = 50, ref_ht = params$rec.ht.mean)
  names(outputs) <- c("P", "F", "K", "mesh")
  lambdas[[i+1]] <- lambda(outputs[["K"]])
  paste("s0 = ", s0s[[i+1]], "; Lambda = ", lambdas[[i+1]])
  setTxtProgressBar(pb, i)
  i <- i + 1
}
close(pb)

# Set s0 as iteration output
params$s0 <- s0s[[length(s0s)]]

# Get traits for IPM
traits <- get.traits(build.ipm(params = params, mesh = 50, ref_ht = params$rec.ht.mean))
sens_lambda_df[which(sens_lambda_df$lambda == 1.15),]$s0 <- params$s0
sens_lambda_df[which(sens_lambda_df$lambda == 1.15),]$GenTfun <- traits$GenTfun
sens_lambda_df[which(sens_lambda_df$lambda == 1.15), ]$Lmean <- traits$Lmean
sens_lambda_df[which(sens_lambda_df$lambda == 1.15), ]$La <- traits$La
sens_lambda_df[which(sens_lambda_df$lambda == 1.15), ]$Lamean <- traits$Lamean

# Save sensitivities
# Export final parameters for mistletoe LHTs
write.csv(sens_lambda_df, "Parameters_lambdas.csv")

# Plot sensitivities
# Reshape the data into long format
sens_lambda_df_long <- sens_lambda_df %>%
  pivot_longer(cols = c(s0, GenTfun, Lmean, La, Lamean),
               names_to = "Parameter",
               values_to = "Value")

# Reorder the parameters manually
sens_lambda_df_long$Parameter <- factor(sens_lambda_df_long$Parameter, levels = c("s0", "GenTfun", "Lmean", "La", "Lamean"))

# Convert lambda to a factor (categorical variable)
sens_lambda_df_long$lambda <- factor(sens_lambda_df_long$lambda, levels = c(1.05, 1.10, 1.15))

# Add parameter labels
parameter_labels <- c(
  "s0" = "s0",
  "GenTfun" = "T",
  "Lmean" = "ηe",
  "La" = "Lα",
  "Lamean" = "Lα-ω"
)

# Create the bar plot
Fig_S9 <- ggplot(sens_lambda_df_long, aes(x = lambda, y = Value, fill = as.factor(lambda))) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Parameter, scales = "free_y", labeller = labeller(Parameter = parameter_labels)) + 
  theme_minimal() +
  labs(x = expression(lambda), y = "Value", fill = expression(lambda)) +
  theme(legend.position = "none",
        strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14)) +
  scale_fill_manual(values = c("1.10" = "grey", "1.05" = colours(6)[4], "1.15" = colours(6)[1]))

Fig_S9

# Sensitivity to berry production form ------------------------------------

# To test the sensitivity to berry production form, a logistic relationship was assumed instead, with midpoint of half the largest size
max_size <- max(wrangled_by_yr_df$logArea, na.rm = TRUE)
min_size <- min(wrangled_by_yr_df$logArea, na.rm = TRUE)
# According to Mellado and Zamora, 2014, 1m^2 (10000cm^2) of mistletoe produces ~2000 berries, suggesting that 2000 = k/(1+exp((-1)*(ln(10000)-0.5*(max_size+min_size_rep))))
k <- 2000*(1+exp((-1)*(log(10000)-0.5*(max_size + min_size_rep))))

# Plot berry number against size for logistic relationship
sizes <- seq(min_size, max_size, length.out = 50)
berry_nos <- k / (1 + exp((-1) * (sizes - 0.5 * (max_size + min_size_rep))))
Fig_S11 <- ggplot(data = data.frame(sizes, berry_nos), aes(x = sizes, y = berry_nos)) +
  geom_line(col = colours(6)[3]) +
  theme_bw() +
  labs(x = xl, y = "Number of berries") +
  theme(axis.title = element_text(size = 15))
Fig_S11

# Amend IPM function to include logistic berry production
build.ipm.logistic <- function(params, mesh, ref_ht){
  
  # Create outputs list
  output_ls <- list()
  
  # Use model coefficients to build vital rate functions
  # Keep mistletoe height constant unless new height generated for a new recruit
  ht_fun <- function (h1, h) {
    if (h1 == h) {return(1)}
    if (h1 != h) {return(0)}
  }
  
  ## Define vital rate functions
  # Survival function
  sur_fun <- function(z){
    u = exp(params$sur.int + params$sur.slope * z)
    return(u/(1+u))
  }
  
  # Growth function
  gro_fun <- function(z1, z) {
    dnorm(z1, mean = params$gro.int + params$gro.slope * z + params$gro.slope2 * z * z, 
          sd = params$gro.sd)
  }
  
  # Fruiting function
  fru_fun <- function(z, h) {
    u = exp(params$fru.int + params$fru.area.slope * z + params$fru.ht.slope * h)
    return(u/(1+u))
  }
  
  # Fecundity (berries produced per fruiting individual), amended for logistic relationship
  fec_fun <- function(z) {
    params$k.berries / ((1+exp((-1)*(z-0.5*(max_size + min_size_rep)))))
  }
  
  # 1-year-old seedling heights
  S1_ht_fun <- function(h1){
    dnorm(h1, mean = params$rec.ht.mean, sd = params$rec.ht.sd)
  }
  
  # Recruit sizes
  R_z_fun <-function(z1){
    dnorm(z1, mean = params$rec.area.mean, sd = params$rec.area.sd)
  }
  
  ## Construct sub-kernels
  ### P sub-kernel
  # Describes probability of an individual of size z and height h becoming a individual of size z1 and height h1
  
  # P sub-kernel - sur * gro 
  P_zh <- function(z1, h1, z, h) {
    return(sur_fun(z) * gro_fun(z1, z) * ht_fun(h1, h))
  }
  
  # Number of meshpoints
  meshz <- mesh
  meshh <- mesh
  
  # Set limits for integration - may impact eviction
  Lz <- 0.95 * (params$rec.area.mean - 4 * params$rec.area.sd)
  Uz <- 1.05 * max(wrangled_by_yr_df$logArea, na.rm = TRUE)
  Lh <- 0.95 * min(wrangled_by_yr_df$Height, na.rm = TRUE)
  Uh <- 1.05 * max(wrangled_by_yr_df$Height, na.rm = TRUE)
  
  # Set bin size and midpoints for each trait variable
  
  # Log area bin size
  hz <- (Uz - Lz) / meshz
  
  # Log area midpoints
  yz <- Lz + hz * ((1:meshz) - 0.5)
  
  # Height bin size
  hh <- (Uh - Lh) / meshh
  
  # Height midpoints
  yh <- Lh + hh * ((1:meshh) - 0.5)
  
  # Create function to fun through mesh points
  eta_ij <- function(i, j, meshz) {(j - 1) * meshz + i}
  
  # Create matrix of indices with which to evaluate mega matrix
  Eta <- outer(1:meshz, 1:meshh, eta_ij, meshz = meshz); 
  
  # Create mega matrix of 0s in which to put transition probabilities at each area-height combination
  P <- matrix(0, meshz*meshh, meshz*meshh)
  
  # Create array of 0s in which to transition probabilities at each area-height combination
  Pvals <- array(0, c(meshz, meshh, meshz, meshh))
  
  # Run through each set of meshpoints and evaluate transition probabilities
  for(i in 1:meshz){
    for(j in 1:meshh){
      for(k in 1:meshh){
        pvals = P_zh(yz, yh[k], yz[i], yh[j])
        P[Eta[,k], Eta[i,j]] = pvals
        Pvals[,k,i,j] = pvals
      }
    }
  }
  
  # Multiply matrix by bin sizes to standardise transition probabilities
  P <- hz * hh * P
  
  # Plot P kernel at mean height and area
  image(t(log(P[(meshz*meshh/2+1):(meshz*meshh/2+meshz), (meshz*meshh/2+1):(meshz*meshh/2+meshz)] + 0.01)))
  
  ### F sub-kernel
  # Describes number of S1 from an individual of size z and height h establishing at height h1
  # CS1 sub-kernel - fru prob * fec * estab prob * new ht
  F_zh <- function(z1, h1, z, h) {
    return(fru_fun(z, h) * fec_fun(z) * params$s0 * S1_ht_fun(h1) * R_z_fun(z1))
  }
  
  # Create function to fun through mesh points
  eta_ij <- function(i, j, meshz) {(j - 1) * meshz + i}
  
  # Create matrix of indices with which to evaluate mega matrix
  Eta <- outer(1:meshz, 1:meshh, eta_ij, meshz = meshz); 
  
  # Modify kernel sizes as z1 is discrete
  # Create mega matrix of 0s in which to put transition probabilities at each area-height combination
  F <- matrix(0, meshz * meshh, meshz * meshh)
  
  # Create array of 0s in which to transition probabilities at each area-height combination
  Fvals <- array(0, c(meshz, meshh, meshz, meshh))
  
  # Run through each set of meshpoints and evaluate transition probabilities
  for(i in 1:meshz){
    for(j in 1:meshh){
      for(k in 1:meshh){
        fvals = F_zh(yz, yh[k], yz[i], yh[j])
        F[Eta[,k], Eta[i,j]] = fvals
        Fvals[,k,i,j] = fvals
      }
    }
  }
  
  # Multiply matrix by bin sizes to standardise transition probabilities
  F <- hz * hh * F
  
  # Plot F kernel at mean height and area
  image(t(log(F[(meshz*meshh/2+1):(meshz*meshh/2+meshz), (meshz*meshh/2+1):(meshz*meshh/2+meshz)] + 0.01)))
  
  # Find bins of min size at reproduction and set non-reproductive bins to 0
  min_size_rep_bin <- ceiling((min_size_rep - Lz) / (Uz - Lz) * mesh)
  
  # Non-reproductive bins
  non_rep_bins <- c(sapply(seq(1, meshz*meshh, by = mesh), function(x) seq(x, x + (min_size_rep_bin - 1))))
  
  # Designate non-reproductive bins as 0
  for(bin in non_rep_bins){
    F[, bin] <- 0
  }
  
  # Combine P and F into mega matrix kernel K
  K <- P + F
  
  # Output P kernel
  output_ls[[1]] <- P
  # Output F kernel
  output_ls[[2]] <- F
  # Output full kernel
  output_ls[[3]] <- K
  # Output meshes
  output_ls[[4]] <- mesh
  
  # Rename outputs
  names(output_ls) <- c("P", "F", "K", "mesh")
  
  # Return outputs
  return(output_ls)
  
  # Find bin of reference height
  ref_ht_bin <- as.integer((ref_ht - Lh)/(Uh - Lh) * meshh)
  
  # Plot K kernel at reference height
  image(t(log(K[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)] + 0.01)))
  
  # Create P kernel at reference height
  ref_P <- P[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)]
  
  # Create F kernel at reference height
  ref_F <- F[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)]
  
  # Plot full kernel
  ref_K <- ref_P + ref_F
  image(t(ref_K))
  
}

# Reset parameters
params <- final_params
params$k <- k

# Reconstruct IPM
outputs <- build.ipm.logistic(params = params, mesh = 50, ref_ht = mean_ht)

# Extract traits
logistic_traits <- get.traits(outputs)

# Remove kernels
traits[["P"]] <- NA
traits[["F"]] <- NA
traits[["K"]] <- NA
logistic_traits[["P"]] <- NA
logistic_traits[["F"]] <- NA
logistic_traits[["K"]] <- NA

# Create data frame of traits for exponential and logistic functional forms
exp_vs_log_df <- data.frame(Berry_model = c(rep("Exponential", 6), rep("Logistic", 6)),
           Trait = rep(c("Lambda", "R0", "GenTfun", "Lmean", "La", "Lamean")),
           Value = c(unlist(traits)[5:10], unlist(logistic_traits)[5:10]))

# Plot traits for exponential and logistic functional forms
Fig_S12 <- ggplot(data = exp_vs_log_df, aes(x = Trait, y = Value, fill = Berry_model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  labs(y = "Value", x = "IPM output", fill = "Berry production function") +
  scale_x_discrete(
    limits = c("Lambda", "R0", "GenTfun", "Lmean", "La", "Lamean"),  # order
    labels = c(
      "Lambda" = expression(lambda),
      "R0" = expression(R[0]),
      "GenTfun" = expression("T"),
      "Lmean" = expression(eta["e"]),
      "La" = expression("L"[alpha]),
      "Lamean" = expression("L"[paste(alpha, "-", omega)])
    )) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  scale_fill_manual(values = c("Exponential" = colours(6)[3], "Logistic" = colours(6)[2]))
Fig_S12

# Plot IPM at different heights -------------------------------------------
# Reset parameters
params <- final_params

# Specify heights at which to plot IPM
hts <- c(params$rec.ht.mean - 2 * params$rec.ht.sd,
         params$rec.ht.mean - params$rec.ht.sd,
         params$rec.ht.mean,
         params$rec.ht.mean + params$rec.ht.sd,
         params$rec.ht.mean + 2 * params$rec.ht.sd)

# Plot IPMs
outputs <- build.ipm(params = params, mesh = 50, ref_ht = hts[1])
Fig_S7a <- plot.ipm(P_kernel = outputs[["P"]], F_kernel = outputs[["F"]], ref_ht = hts[1], mesh = outputs[["mesh"]])
Fig_S7a
outputs <- build.ipm(params = params, mesh = 50, ref_ht = hts[2])
Fig_S7b <- plot.ipm(P_kernel = outputs[["P"]], F_kernel = outputs[["F"]], ref_ht = hts[2], mesh = outputs[["mesh"]])
Fig_S7b
outputs <- build.ipm(params = params, mesh = 50, ref_ht = hts[3])
Fig_S7c <- plot.ipm(P_kernel = outputs[["P"]], F_kernel = outputs[["F"]], ref_ht = hts[3], mesh = outputs[["mesh"]])
Fig_S7c
outputs <- build.ipm(params = params, mesh = 50, ref_ht = hts[4])
Fig_S7d <- plot.ipm(P_kernel = outputs[["P"]], F_kernel = outputs[["F"]], ref_ht = hts[4], mesh = outputs[["mesh"]])
Fig_S7d
outputs <- build.ipm(params = params, mesh = 50, ref_ht = hts[5])
Fig_S7e <- plot.ipm(P_kernel = outputs[["P"]], F_kernel = outputs[["F"]], ref_ht = hts[5], mesh = outputs[["mesh"]])
Fig_S7e
