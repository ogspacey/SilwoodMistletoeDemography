### Silwood Mistletoe Vital Rate Regressions, Trade-offs and IPM
### Oliver G. Spacey, Owen R. Jones, Sydne Record, Arya Y. Yue, Wenyi Liu, Alice Rosen, Michael Crawley, Chris J. Thorogood, Roberto Salguero-Gómez
### 21.03.2024

# Pre-amble ---------------------------------------------------------------
# Clear environment
rm(list = ls())

# Load packages
library(ggplot2)      # Data visualisation
library(patchwork)    # Plot arrangement
library(tidyverse)    # Data wrangling
library(lme4)         # Regression models
library(lmerTest)     # Regression models
library(RColorBrewer) # Colours for visualisation
library(khroma)       # Colours for visualisation
library(popbio)       # Demographic analysis
library(Rage)         # Demographic analysis

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
for(i in 1:nrow(mst_hst_df)) {
  for(j in 1:(ncol(mst_hst_df)-1)){
    if(mst_hst_df[i,j] %in% c("FC", "Surv") == TRUE) {
      known_presences <- known_presences + 1
      if(is.na(mst_hst_df[i,j+1]) == TRUE){
      NAs <- NAs + 1
      }
    }
    else { count <- count }
}
}
# 454 NAs out of 3183 detections

# Estimate detectability
NAs/known_presences
# 14.3% of known presences not measured

### Sensitivity when these records removed - remove individuals or just these years?

# Estimate proportion of mistletoes for which height is calculated directly
nrow(filter(mst_hst_df, Height > 0))/nrow(mst_hst_df)
# 39.1% of mistletoe heights measured directly

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
ggplot(data = mst_hst_df, aes(x = Combined_height)) +
  geom_histogram() +
  labs(x = "Mistletoe height above ground (m)") +
  theme_bw()

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

# Visualise relationship between host and mistletoe heights
ggplot(data = mst_hst_df, aes(x = Host_height, y = Combined_height)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm") +
  labs(x = "Host height", y = "Mistletoe height")
# Taller the host, higher up mistletoes generally are

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

# Re-visualise distribution of mistletoe areas
ggplot(data = wrangled_by_yr_df, aes(x = logArea)) +
  geom_histogram() +
  theme_bw()

# Outliers of area probably due to measurement error, unlike height outliers
# Do not remove outliers for height

# Visualise relationship between areas and height
ggplot(data = wrangled_by_yr_df, aes(x = logArea, y = Height)) +
  geom_point() +
  theme_bw()

# Re-test relationship between height and area

# Mixed effect model of area and height  with individual ID as random effect
area_height_mem <- lmer(Height ~ logArea + (1 | Indiv_ID), data = wrangled_by_yr_df)
summary(area_height_mem)
# No more significant relationship between height and area after correction

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

# Model Survival ~ log(Area) + Height with individual ID as random effect
sur_area_ht_glmm <- glmer(Survive ~ logArea_t0 + Height + (1 | Indiv_ID), 
                     data = wrangled_df, family = binomial)
summary(sur_area_ht_glmm)
# No significant effect of either Area or Height

# Plot Survival ~ Stage
ggplot(data = wrangled_df, aes(x = Stage, y = Survive)) +
  geom_jitter(width = 0.2, height = 0.1) +
  labs(x = "Stage", y = "Probability of survival to t1") +
  theme_bw()

# Model Survival ~ Stage with individual ID as random effect
sur_stage_glmm <- glmer(Survive ~ Stage + (1 | Indiv_ID), 
                     data = wrangled_df, family = binomial)
summary(sur_stage_glmm)
# No significant effect of stage

# Plot Survival ~ log(Area) + Stage
ggplot(data = wrangled_df, aes(x = logArea_t0, y = Survive, col = Stage)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Figure 2a - Plot Survival ~ log(Area) + Stage
Fig_2a <- ggplot(data = wrangled_df, aes(x = logArea_t0, y = Survive, col = Stage)) +
          geom_point(size = 1, alpha = 0.5) +
          stat_smooth(method = "glm", 
                      method.args = list(family = binomial)) +
          labs(x = "", y = "Survival probability (t to t+1)") +
          scale_color_manual(values = c("Juvenile" = colours(6)[6], "Adult" = colours(6)[5])) + 
          geom_vline(xintercept = min_size_rep, color = "grey", linetype = "dashed") +
          theme_bw() + 
          theme(legend.position = c(0.6,0.4),
                legend.title = element_text(size = 6),
                legend.text = element_text(size = 6),
                legend.key.size = unit(0.4, units = "cm"),
                axis.title.y = element_text(size = 7))
Fig_2a

# Model Survival ~ log(Area) with individual ID as random effect
sur_area_stage_glmm <- glmer(Survive ~ logArea_t0 + Stage + (1 | Indiv_ID), 
                       data = wrangled_df, family = binomial)
summary(sur_area_stage_glmm)

# Plot Survival ~ Height + Stage
ggplot(data = wrangled_df, aes(x = Height, y = Survive, col = Stage)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Figure 2b - Plot Survival ~ Height + Stage
Fig_2b <- ggplot(data = wrangled_df, aes(x = Height, y = Survive, col = Stage)) +
          geom_point(size = 1, alpha = 0.5) +
          stat_smooth(method = "glm", 
                      method.args = list(family = binomial)) +
          labs(x = "", y = "") +
          scale_color_manual(values = c("Juvenile" = colours(6)[6], "Adult" = colours(6)[5])) + 
          theme_bw() +
          theme(legend.position = "none")
Fig_2b

# Model Survival ~ Height with individual ID as random effect
sur_ht_stage_glmm <- glmer(Survive ~ Height + Stage + (1 | Indiv_ID), 
                             data = wrangled_df, family = binomial)
summary(sur_ht_stage_glmm)


# Model Survival ~ log(Area) + Height with individual ID as random effect
sur_area_ht_stage_glmm <- glmer(Survive ~ logArea_t0 + Height + Stage + (1 | Indiv_ID), 
                           data = wrangled_df, family = binomial)
summary(sur_area_ht_stage_glmm)

# List models for survival
sur_mods <- c(sur_area_glmm,
              sur_ht_glmm,
              sur_area_ht_glmm,
              sur_stage_glmm, 
              sur_area_stage_glmm,
              sur_ht_stage_glmm,
              sur_area_ht_stage_glmm)

# Create model selection table for survival
sur_mod_sel <- bind_cols(data.frame(Vital_rate = "Survival"),
                             mod.sel.table(sur_mods))
# Neither area, height nor stage predict survival across all mistletoes

# Model juvenile and adult survival separately to better reflect life cycle
# Filter wrangled data for juveniles only
wrangled_juv_df <- filter(wrangled_df, 
                          Stage == "Juvenile")

# Filter wrangled data for adults only
wrangled_adu_df <- filter(wrangled_df, 
                          Stage == "Adult")

### Juvenile survival -------------------------------------------------------
# Plot Survival ~ log(Area) for juveniles
ggplot(data = wrangled_juv_df, aes(x = logArea_t0, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ log(Area) for juveniles with individual ID as random effect
sur_juv_area_glmm <- glmer(Survive ~ logArea_t0 + (1 | Indiv_ID), 
                           data = wrangled_juv_df, family = binomial)
summary(sur_juv_area_glmm)

# Plot Survival ~ Height for juveniles
ggplot(data = wrangled_juv_df, aes(x = Height, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ Height for juveniles with individual ID as random effect
sur_juv_ht_glmm <- glmer(Survive ~ Height + (1 | Indiv_ID), 
                         data = wrangled_juv_df, family = binomial)
summary(sur_juv_ht_glmm)

# Plot Survival ~ log(Area) + Height for juveniles
ggplot(data = wrangled_juv_df, aes(x = logArea_t0, y = Survive, col = Height)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ log(Area) + Height for juveniles with Indiv_ID as random effect
sur_juv_area_ht_glmm <- glmer(Survive ~ logArea_t0 + Height + (1 | Indiv_ID), 
                              data = wrangled_juv_df, family = binomial)
summary(sur_juv_area_ht_glmm)

# Model Survival ~ log(Area) * Height for juveniles with Indiv_ID as random effect
sur_juv_area_ht_int_glmm <- glmer(Survive ~ logArea_t0 * Height + (1 | Indiv_ID), 
                                  data = wrangled_juv_df, family = binomial)
summary(sur_juv_area_ht_int_glmm)

# List models for juvenile survival
sur_juv_mods <- c(sur_juv_area_glmm,
                  sur_juv_ht_glmm,
                  sur_juv_area_ht_glmm,
                  sur_juv_area_ht_int_glmm)

# Create model selection table for juvenile survival
sur_juv_mod_sel <- bind_cols(data.frame(Vital_rate = "Juvenile survival"),
                             mod.sel.table(sur_juv_mods))

### Adult survival -------------------------------------------------------
# Plot Survival ~ log(Area) for adults
ggplot(data = wrangled_adu_df, aes(x = logArea_t0, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "log(Area)", y = "Survival from t0 to t1") +
  theme_bw()

# Model Survival ~ log(Area) for adults with individual ID as random effect
sur_adu_area_glmm <- glmer(Survive ~ logArea_t0 + (1 | Indiv_ID), 
                           data = wrangled_adu_df, family = binomial)
summary(sur_adu_area_glmm)

# Plot Survival ~ Height for adults
ggplot(data = wrangled_adu_df, aes(x = Height, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "Height", y = "Survival from t0 to t1") +
  theme_bw()

# Model Survival ~ Height for adults with Indiv_ID as random effect
sur_adu_ht_glmm <- glmer(Survive ~ Height + (1 | Indiv_ID), 
                         data = wrangled_adu_df, family = binomial)
summary(sur_adu_ht_glmm)

# Plot Survival ~ log(Area) + Height for adults
ggplot(data = wrangled_adu_df, aes(x = logArea_t0, y = Survive, col = Height)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ log(Area) + Height for adults with Indiv_ID as random effect
sur_adu_area_ht_glmm <- glmer(Survive ~ logArea_t0 + Height + (1 | Indiv_ID), 
                              data = wrangled_adu_df, family = binomial)
summary(sur_adu_area_ht_glmm)

# Model Survival ~ log(Area) * Height for adults with Indiv_ID as random effect
sur_adu_area_ht_int_glmm <- glmer(Survive ~ logArea_t0 * Height + (1 | Indiv_ID), 
                                  data = wrangled_adu_df, family = binomial)
summary(sur_adu_area_ht_int_glmm)

# List models for adult survival
sur_adu_mods <- c(sur_adu_area_glmm,
                  sur_adu_ht_glmm,
                  sur_adu_area_ht_glmm,
                  sur_adu_area_ht_int_glmm)

# Create model selection table for adult survival
sur_adu_mod_sel <- bind_cols(data.frame(Vital_rate = "Adult survival"),
                             mod.sel.table(sur_adu_mods))

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

# Plot relative growth from t0 to t1 against Height
ggplot(data = wrangled_df, aes(x = Height, y = Growth_t0_t1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() 

# Model Growth ~ Height with Indiv_ID as random effect
gro_ht_mem <- lmer(Growth_t0_t1 ~ Height + (1 | Indiv_ID), 
                   data = wrangled_df)
summary(gro_ht_mem)

# Model Growth ~ Height with Indiv_ID as random effect
gro_area_ht_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + Height + (1 | Indiv_ID), 
                   data = wrangled_df)
summary(gro_area_ht_mem)

# Plot relative growth from t0 to t1 against stage
ggplot(data = wrangled_df, aes(x = Stage, y = Growth_t0_t1)) +
  geom_violin() +
  geom_point() +
  theme_bw() 

# Model Growth ~ Stage with Indiv_ID as random effect
gro_stage_mem <- lmer(Growth_t0_t1 ~ Stage + (1 | Indiv_ID), 
                      data = wrangled_df)
summary(gro_stage_mem)

# Figure 2c - Plot Growth ~ log(Area) + Stage
Fig_2c <- ggplot(data = wrangled_df, aes(x = logArea_t0, y = Growth_t0_t1, col = Stage)) +
          geom_point(size = 1, alpha = 0.5) +
          geom_smooth(data = subset(wrangled_df, Stage == "Juvenile"), method = "lm", formula = y ~ x) +
          geom_smooth(data = subset(wrangled_df, Stage == "Adult"), method = "lm", formula = y ~ poly(x, 2)) +
          labs(x = "", y = "RGR (t to t+1)") +
          scale_color_manual(values = c("Juvenile" = colours(6)[6], "Adult" = colours(6)[5])) + 
          geom_vline(xintercept = min_size_rep, color = "grey", linetype = "dashed") +
          theme_bw() +
          theme(legend.position = "none")
Fig_2c

# Model Growth ~ log(Area) + Stage with Indiv_ID as random effect
gro_area_stage_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + Stage + (1 | Indiv_ID), 
                         data = wrangled_df)
summary(gro_area_stage_mem)

# Figure 2d - Plot Growth ~ Height + Stage
Fig_2d <- ggplot(data = wrangled_df, aes(x = Height, y = Growth_t0_t1, col = Stage)) +
          geom_point(size = 1, alpha = 0.5) +
          geom_smooth(data = subset(wrangled_df, Stage == "Juvenile"), method = "lm", formula = y ~ x) +
          geom_smooth(data = subset(wrangled_df, Stage == "Adult"), method = "lm", formula = y ~ x) +
          labs(x = "", y = "") +
          scale_color_manual(values = c("Juvenile" = colours(6)[6], "Adult" = colours(6)[5])) + 
          theme_bw() +
          theme(legend.position = "none")
Fig_2d

# Model Growth ~ Height + Stage with Indiv_ID as random effect
gro_ht_stage_mem <- lmer(Growth_t0_t1 ~ Height + Stage + (1 | Indiv_ID), 
                     data = wrangled_df)
summary(gro_ht_stage_mem)

# Model Growth ~ log(Area) + Height + Stage with Indiv_ID as random effect
gro_area_ht_stage_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + Height + Stage + (1 | Indiv_ID), 
                   data = wrangled_df)
summary(gro_area_ht_stage_mem)

# List models for growth
gro_mods <- c(gro_area_mem,
              gro_area_quad_mem,
              gro_ht_mem,
              gro_area_ht_mem,
              gro_stage_mem, 
              gro_area_stage_mem,
              gro_ht_stage_mem,
              gro_area_ht_stage_mem)

# Create model selection table for growth
gro_mod_sel <- bind_cols(data.frame(Vital_rate = "Growth"),
                         mod.sel.table(gro_mods))
# Quadratic function of area best to predict growth across all mistletoes

# Model juvenile and adult growth separately to better reflect life cycle
# Filter wrangled data for juveniles only
wrangled_juv_df <- filter(wrangled_df, 
                          Stage == "Juvenile")

# Filter wrangled data for adults only
wrangled_adu_df <- filter(wrangled_df, 
                          Stage == "Adult")

### Juvenile growth -------------------------------------------------------
# Plot Growth ~ log(Area) for juveniles
ggplot(data = wrangled_juv_df, aes(x = logArea_t0, y = Growth_t0_t1)) +
  geom_point() +
  labs(x = "log(Area)", y = "Relative change in size from t0 to t1") +
  stat_smooth(method = "lm") +
  theme_bw()

# Model Growth ~ log(Area) for juveniles with individual ID as random effect
gro_juv_area_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + (1 | Indiv_ID), 
                           data = wrangled_juv_df)
summary(gro_juv_area_mem)

# Plot quadratic fit for area for juveniles
ggplot(data = wrangled_juv_df, aes(x = logArea_t0, y = Growth_t0_t1)) +
  geom_point() +
  labs(x = "log(Area)", y = "Relative change in size from t0 to t1") +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  theme_bw()

# Model Growth ~ log(Area) + log(Area)^2 for juveniles with individual ID as random effect
gro_juv_area_quad_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + I(logArea_t0^2) + (1 | Indiv_ID), 
                              data = wrangled_juv_df)
summary(gro_juv_area_quad_mem)

# Plot Growth ~ Height for juveniles
ggplot(data = wrangled_juv_df, aes(x = Height, y = Growth_t0_t1)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_bw()

# Model Growth ~ Height for juveniles with individual ID as random effect
gro_juv_ht_mem <- lmer(Growth_t0_t1 ~ Height + (1 | Indiv_ID), 
                         data = wrangled_juv_df)
summary(gro_juv_ht_mem)

# Plot Growth ~ log(Area) + Height for juveniles
ggplot(data = wrangled_juv_df, aes(x = logArea_t0, y = Growth_t0_t1, col = Height)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_bw()

# Model Growth ~ log(Area) + Height for juveniles with individual ID as random effect
gro_juv_area_ht_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + Height + (1 | Indiv_ID), 
                         data = wrangled_juv_df)
summary(gro_juv_area_ht_mem)

# Model Growth ~ log(Area) + Height for juveniles with individual ID as random effect
gro_juv_area_ht_int_mem <- lmer(Growth_t0_t1 ~ logArea_t0 * Height + (1 | Indiv_ID), 
                            data = wrangled_juv_df)
summary(gro_juv_area_ht_int_mem)

# List models for juvenile growth
gro_juv_mods <- c(gro_juv_area_mem,
                  gro_juv_area_quad_mem,
                  gro_juv_ht_mem,
                  gro_juv_area_ht_mem,
                  gro_juv_area_ht_int_mem)

# Create model selection table for juvenile growth
gro_juv_mod_sel <- bind_cols(data.frame(Vital_rate = "Juvenile growth"),
                             mod.sel.table(gro_juv_mods))

### Adult growth -------------------------------------------------------
# Plot Growth ~ log(Area) for adults
ggplot(data = wrangled_adu_df, aes(x = logArea_t0, y = Growth_t0_t1)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_bw()

# Model Growth ~ log(Area) for adults with individual ID as random effect
gro_adu_area_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + (1 | Indiv_ID), 
                         data = wrangled_adu_df)
summary(gro_adu_area_mem)

# Plot quadratic fit for area for adults
ggplot(data = wrangled_adu_df, aes(x = logArea_t0, y = Growth_t0_t1)) +
  geom_point() +
  labs(x = "log(Area)", y = "Relative change in size from t0 to t1") +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2)) +
  theme_bw()

# Model Growth ~ log(Area) + log(Area)^2 for adults with individual ID as random effect
gro_adu_area_quad_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + I(logArea_t0^2) + (1 | Indiv_ID), 
                         data = wrangled_adu_df)
summary(gro_adu_area_quad_mem)

# Plot Growth ~ Height for adults
ggplot(data = wrangled_adu_df, aes(x = Height, y = Growth_t0_t1)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_bw()

# Model Growth ~ Height for adults with individual ID as random effect
gro_adu_ht_mem <- lmer(Growth_t0_t1 ~ Height + (1 | Indiv_ID), 
                       data = wrangled_adu_df)
summary(gro_adu_ht_mem)

# Plot Growth ~ log(Area) + Height for adults
ggplot(data = wrangled_adu_df, aes(x = logArea_t0, y = Growth_t0_t1, col = Height)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_bw()

# Model Growth ~ log(Area) + Height for adults with individual ID as random effect
gro_adu_area_ht_mem <- lmer(Growth_t0_t1 ~ logArea_t0 + Height + (1 | Indiv_ID), 
                            data = wrangled_adu_df)
summary(gro_adu_area_ht_mem)

# Model Growth ~ log(Area) + Height for adults with individual ID as random effect
gro_adu_area_ht_int_mem <- lmer(Growth_t0_t1 ~ logArea_t0 * Height + (1 | Indiv_ID), 
                                data = wrangled_adu_df)
summary(gro_adu_area_ht_int_mem)

# List models for adult growth
gro_adu_mods <- c(gro_adu_area_mem,
                  gro_adu_area_quad_mem,
                  gro_adu_ht_mem,
                  gro_adu_area_ht_mem,
                  gro_adu_area_ht_int_mem)

# Create model selection table for adult growth
gro_adu_mod_sel <- bind_cols(data.frame(Vital_rate = "Adult growth"),
                             mod.sel.table(gro_adu_mods))

# Improvement of model fit by addition of quadratic term
AIC(gro_adu_area_mem) - AIC(gro_adu_area_quad_mem)

## Model maturation --------------------------------------------------------
# Split surviving juveniles into those which matured (i.e., of reproductive size in t1) and those which did not
wrangled_juv_sur_df <- mutate(wrangled_juv_df[wrangled_juv_df$Survive == 1,], 
                              Mature = case_when(logArea_t1 > min_size_rep ~ 1,
                                                 .default = 0))

# Figure S4
# Plot maturation probability as a function of log(Area) 
Fig_S4 <- ggplot(data = wrangled_juv_sur_df, aes(x = logArea_t0, y = Mature)) +
          geom_point(alpha = 0.3, col = colours(6)[6]) +
          stat_smooth(method = "glm", 
                      method.args = list(family = binomial), 
                      col = colours(6)[6]) +
          labs(x = "log(Area) in t", y = "probability of maturing in t+1") +
          theme_bw() +
          theme(text = element_text(size = 14),           # Increase overall text size
                axis.title = element_text(size = 16),     # Increase axis titles size
                axis.text = element_text(size = 14))
Fig_S4

# Model Mature ~ log(Area) for juveniles with Indiv_ID as random effect
mat_area_glmm <- glmer(Mature ~ logArea_t0 + (1 | Indiv_ID), 
                       data = wrangled_juv_sur_df, family = binomial)
summary(mat_area_glmm)

# Plot Maturation ~ Height for juveniles
ggplot(data = wrangled_juv_sur_df, aes(x = Height, y = Mature)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Mature ~ Height for juveniles with Indiv_ID as random effect
mat_ht_glmm <- glmer(Mature ~ Height + (1 | Indiv_ID), 
                     data = wrangled_juv_sur_df, family = binomial)
summary(mat_ht_glmm)

# Plot Maturation ~ log(Area) + Height for juveniles
ggplot(data = wrangled_juv_sur_df, aes(x = logArea_t0, y = Mature, col = Height)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Mature ~ log(Area) + Height for juveniles with Indiv_ID as random effect
mat_area_ht_glmm <- glmer(Mature ~ logArea_t0 + Height + (1 | Indiv_ID), 
                          data = wrangled_juv_sur_df, family = binomial)
summary(mat_area_ht_glmm)

# Model Mature ~ log(Area) * Height for juveniles with Indiv_ID as random effect
mat_area_ht_int_glmm <- glmer(Mature ~ logArea_t0 * Height + (1 | Indiv_ID), 
                              data = wrangled_juv_sur_df, family = binomial)
summary(mat_area_ht_int_glmm)

# List models for juvenile maturation
mat_mods <- c(mat_area_glmm,
              mat_ht_glmm,
              mat_area_ht_glmm,
              mat_area_ht_int_glmm)

# Create model selection table for maturation
mat_mod_sel <- bind_cols(data.frame(Vital_rate = "Maturation"),
                         mod.sel.table(mat_mods))

## Model fruiting ------------------------------------------------------
# Filter wrangled data by year for adults only - fruiting and non-fruiting
wrangled_by_yr_adu_df <- filter(wrangled_by_yr_df, 
                          Stage == "Adult")

# Figure 2e
# Plot fruiting ~ logArea for adults only 
Fig_2e <- ggplot(data = wrangled_by_yr_adu_df, aes(x = logArea, y = Fruit)) +
          geom_point(size = 1, alpha = 0.5, col = colours(6)[3]) +
          stat_smooth(method = "glm", 
                      method.args = list(family = binomial),
                      col = colours(6)[3]) +
          labs(x = "log(Area) in t", y = "Fruiting probability (t)") +
          theme_bw() +
          theme(axis.title.y = element_text(size = 7))
Fig_2e

# Model fruiting ~ logArea with Indiv_ID as random effect
fru_area_glmm <- glmer(Fruit ~ logArea + (1 | Indiv_ID), 
                       data = wrangled_by_yr_adu_df, family = binomial)
summary(fru_area_glmm)

# Plot fruiting ~ Height for adults only 
ggplot(data = wrangled_by_yr_adu_df, aes(x = Height, y = Fruit)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "Height above ground (m)", y = "Probability of fruiting in t0") +
  theme_bw()

# Figure 2f
# Plot fruiting ~ Height for adults only 
Fig_2f <- ggplot(data = wrangled_by_yr_adu_df, aes(x = Height, y = Fruit)) +
          geom_point(size = 1, alpha = 0.5, col = colours(6)[3]) +
          stat_smooth(method = "glm", 
                      method.args = list(family = binomial),
                      col = colours(6)[3]) +
          labs(x = "Height above ground (m)", y = "") +
          theme_bw()
Fig_2f

# Model fruiting ~ Height with Indiv_ID as random effect
fru_ht_glmm <- glmer(Fruit ~ Height + (1 | Indiv_ID), 
                     data = wrangled_by_yr_adu_df, family = binomial)
summary(fru_ht_glmm)

# Plot fruiting ~ log(Area) + Height for adults only 
ggplot(data = wrangled_by_yr_adu_df, aes(x = logArea, y = Fruit, col = Height)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "log(Area) in t0", y = "Probability of fruiting in t0") +
  theme_bw()

# Model fruiting ~ logArea + Height with Indiv_ID as random effect
fru_area_ht_glmm <- glmer(Fruit ~ logArea + Height + (1 | Indiv_ID), 
                          data = wrangled_by_yr_adu_df, family = binomial)
summary(fru_area_ht_glmm)

# Model fruiting ~ logArea * Height with Indiv_ID as random effect
fru_area_ht_int_glmm <- glmer(Fruit ~ logArea * Height + (1 | Indiv_ID), 
                              data = wrangled_by_yr_adu_df, family = binomial)
summary(fru_area_ht_int_glmm)

# List models for fruiting
fru_mods <- c(fru_area_glmm,
              fru_ht_glmm,
              fru_area_ht_glmm,
              fru_area_ht_int_glmm)

# Create model selection table for fruiting
fru_mod_sel <- bind_cols(data.frame(Vital_rate = "Fruiting"),
                         mod.sel.table(fru_mods))


## Overall model selection -------------------------------------------------
# Table S1
# Combine all model selection tables
mod_sel <- bind_rows(sur_mod_sel, sur_juv_mod_sel, sur_adu_mod_sel,
                     gro_mod_sel, gro_juv_mod_sel, gro_adu_mod_sel,
                     mat_mod_sel, fru_mod_sel)

# Find minimum AIC values
min_AIC <- mod_sel %>%
           group_by(Vital_rate) %>%
           filter(AIC == case_when(AIC > 0 ~ min(AIC),
                                   AIC < 0 ~ max(AIC))) %>%
           select(AIC, Vital_rate)

# Show which models selected
mod_sel <- mod_sel %>%
           mutate(Selected = case_when(AIC %in% min_AIC$AIC ~ TRUE, 
                                       .default = FALSE))

# Export as .csv
write.csv(mod_sel, "Table_S1.csv")


# Figure 2 ----------------------------------------------------------------
# Combine all vital rate regressions into Figure 2
Fig_2 <- (Fig_2a + Fig_2b) / (Fig_2c + Fig_2d) / (Fig_2e + Fig_2f) + 
         plot_annotation(tag_levels = "a") &
         theme(plot.tag = element_text(size = 10))
Fig_2


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
Fig_3a <- ggplot(data = wrangled_to_df, aes(x = Growth_t0_t1, y = Survive_t1_t2)) +
          geom_point() +
          stat_smooth(method = "glm", method.args = list(family = binomial)) +
          labs(x = "", y = "Survival probability (t+1 to t+2)") +
          theme_bw()
Fig_3a

# Model survival from t+1 to t+2 as a function of growth from t to t+1 for all individuals with Indiv_ID as random effect
sur_gro_to_glmm <- glmer(Survive_t1_t2 ~ Growth_t0_t1 + (1 | Indiv_ID), 
                              data = wrangled_to_df, family = binomial)
summary(sur_gro_to_glmm)
# When individual ID controlled for, no significant relationship

# Examine trade-offs involving reproduction for adults
# Filter wrangled data for adults known to be female
wrangled_fem_to_df <- filter(wrangled_to_df, Stage == "Adult") %>%
                      group_by(Indiv_ID) %>%
                      filter(any(Fruit_t0 == 1)) %>%
                      ungroup()

# Plot survival from t to t+1 as a function of fruiting in t for females
Fig_3b <- ggplot(data = wrangled_fem_to_df, aes(x = Fruit_t0, y = Survive_t0_t1)) +
          geom_jitter(width = 0.2, height = 0.1, alpha = 0.3) +
          stat_smooth(method = "glm", method.args = list(family = binomial)) +
          labs(x = "", y = "Survival probability (t to t+1)") +
          theme_bw()
Fig_3b

# Model survival from t to t+1 as a function of fruiting in t for females with Indiv_ID as random effect
sur_fru_to_glmm <- glmer(Survive_t0_t1 ~ Fruit_t0 + (1 | Indiv_ID), 
                             data = wrangled_fem_to_df, family = binomial)
summary(sur_fru_to_glmm)
# Significant negative relationship between fruiting and survival

# Plot fruiting in t+1 as a function of growth from t to t+1 for adults
Fig_3c <- ggplot(data = wrangled_fem_to_df, aes(x = Growth_t0_t1, y = Fruit_t1)) +
          geom_point() +
          stat_smooth(method = "glm", method.args = list(family = binomial)) +
          labs(x = "RGR (t to t+1)", y = "Fruiting probability (t+1)") +
          lims(x = c(-0.15, max(wrangled_fem_to_df$Growth_t0_t1, na.rm = TRUE))) +
          theme_bw()
Fig_3c

# Model fruiting in t+1 as a function of growth from t to t+1 for adults with Indiv_ID as random effect
fru_gro_to_glmm <- glmer(Fruit_t1 ~ Growth_t0_t1 + (1 | Indiv_ID), 
                             data = wrangled_fem_to_df, family = binomial)
summary(fru_gro_to_glmm)
# significant negative relationship; more growth, less likely to fruit
# reflects investment in reproduction as mistletoes grow larger

# Plot growth from t to t+1 as a function of fruiting in t between reproduction in t and reproduction in t+1 for adults
Fig_3d <- ggplot(data = wrangled_fem_to_df, aes(x = Fruit_t0, y = Growth_t0_t1)) +
          geom_jitter(width = 0.2, height = 0.1, alpha = 0.3) +
          stat_smooth(method = "lm", se = TRUE) +
          labs(x = "Fruiting probability (t)", y = "RGR (t to t+1)") +
          theme_bw()
Fig_3d

# Model trade-off between reproduction in t and growth from t to t+1 for adults
gro_fru_to_mem <- lmer(Growth_t0_t1 ~ Fruit_t0 + (1 | Indiv_ID), 
                             data = wrangled_fem_to_df)
summary(gro_fru_to_mem)
# Positive relationship between reproduction and fruiting

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

# Export model summary as .csv
write.csv(vr_to_df, "Table S2.csv")

# Combine vital rate trade-off graphs into Figure 3
Fig_3 <- (Fig_3a + Fig_3b) / (Fig_3c + Fig_3d) + plot_annotation(tag_levels = "a")
Fig_3

# Recruitment -------------------------------------------------------------
# Mistletoes are recruited through spread of seed from fruiting individuals to a suitable branch
# Seeds germinate and establish for 3 years before they can be observed (in their 3rd year)
# I.e., it takes three time-steps for fruiting individuals to contribute to the juvenile pool
# As we lack data on germination and establishment success, we estimate an initial value of s0 linking fruiting in t to juvenile recruitment in t+3
# Individuals produce berries directly proportional to the number of terminal shoots, which is directly proportional to their size - b = k * e^(log(a))
# According to Mellado and Zamora, 2014, 1m^2 (10000cm^2) of mistletoe produces ~2000 berries, suggesting that 2000 = k * e(log(10000)) -> k = 0.2
# As per Lucas et al., 2019, Fertility from t to t + 3 = (total no. offspring in t+3)/(total no. berries in t) * (no. berries for individual in t) = Σrec/Σb * b(z)
# Via two intermediary stages: 1-year-old seedlings (S1) and 2-year-old seedlings (S2)
# Σrec(t+3)/Σb(t) = S1(t+1)/b(t) * S2(t+2)/S1(t+1) * rec(t+3)/S2(t+2) = s0 * s1 * s2
# Assuming constant survival as we cannot measure establishment s0, survival of 1-year-old seedlings s1 or survival of 2-year-old seedlings s2, s0 = s1 = s2
# therefore s0 = (Σrec/Σb)^(1/3)
# Σb is estimated using the berry production constant k, where b = k*exp(z)
# s0 will later be estimated more realistically through iteration to obtain λ=1.1

# When the IPM is built using a particular berry production constant k, s0 will be estimated as follows
get.s0 <- function(k){
  # Filter long df for only fruiting individuals and estimate berry production
  fru_df <- filter(wrangled_by_yr_df, Fruit == 1) %>%
            mutate(Berry_no = k*exp(logArea))
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
  s0 <- mean_f_t_t3 ^ (1/3)
  
  # Return s0
  return(s0)
}

# Examine how s0 varies with geometric sequence of ks
geomSeries <- function(a, r, max) {
  a*r^(0:floor(log(max, r)))
}

ks <- geomSeries(a=10^(-5), r=10, max=10^7)
s0_ls <- list()
for(i in 1:length(ks)){
  s0_ls[[i]] <- get.s0(ks[i])
}

s0s <- unlist(s0_ls)
plot(x = ks, y = s0s)

# k = 0.2 suggests s0 ~ 0.05, but survival difficult to estimate, so s0 will be fixed such that λ=1.1
get.s0(0.2)

# Find distribution of new recruit sizes
# Filter only new recruits
rec_only_df <- filter(wrangled_by_yr_df, 
                      Year >= 17 & Status == "FC" & Stage == "Juvenile")

# Figure S5
# Create new x label
xl <- expression(log ~ "area" ~ cm^2)

# Plot distribution of recruit sizes
Fig_S5a <- ggplot(data = rec_only_df, aes(x = logArea)) +
           geom_histogram(aes(y = ..density..), alpha = 0.3, fill = colours(6)[4], bins = 10) +
           stat_function(fun = function(x) {
             dnorm(x, mean = mean(rec_only_df$logArea), sd = sd(rec_only_df$logArea))
            }, color = colours(6)[4]) +
           labs(x = xl, y = "Density") + 
           theme_bw() +
           theme(text = element_text(size = 14),           
                 axis.title = element_text(size = 16),    
                 axis.text = element_text(size = 14))

# Plot distribution of heights for all individuals to get a better height distribution
Fig_S5b <- ggplot(data = wrangled_by_yr_df, aes(x = Height)) +
           geom_histogram(aes(y = ..density..), alpha = 0.3, fill = colours(6)[4], bins = 10) +
           stat_function(fun = function(x) {
             dnorm(x, mean = mean(wrangled_by_yr_df$Height), sd = sd(wrangled_by_yr_df$Height))
             }, color = colours(6)[4]) +
           labs(x = "Height (m)", y = "Density") + 
           theme_bw() +
           theme(text = element_text(size = 14),           
                  axis.title = element_text(size = 16),    
                  axis.text = element_text(size = 14))

# Plot Figure S5
Fig_S5 <- (Fig_S5a + Fig_S5b) + plot_annotation(tag_levels = "a")
Fig_S5

# Build IPM ---------------------------------------------------------------
## Extract parameters
# Extract parameters from selected models
# Cannot used random effect because there is no way of tracking individual ID in an IPM
# Juvenile survival ~ log(Area)
sur_juv_area_glm <- glm(Survive ~ logArea_t0, 
                        data = wrangled_juv_df, family = binomial)
summary(sur_juv_area_glm)

# Juvenile growth ~ log(Area) - quantify as absolute growth rather than relative rate to project next size
gro_abs_juv_area_lm <- lm(logArea_t1 ~ logArea_t0, 
                      data = wrangled_juv_df)
summary(gro_abs_juv_area_lm)

# Juvenile maturation ~ log(Area) * Height
mat_area_ht_int_glm <- glm(Mature ~ logArea_t0 * Height, 
                    data = wrangled_juv_sur_df, family = binomial)
summary(mat_area_ht_int_glm)

# Adult survival ~ log(Area) * Height
sur_adu_area_ht_int_glm <- glm(Survive ~ logArea_t0 * Height, 
                           data = wrangled_adu_df, family = binomial)
summary(sur_adu_area_ht_int_glm)

# Adult growth ~ log(Area) + log(Area)^2
gro_abs_adu_area_quad_lm <- lm(logArea_t1 ~ logArea_t0 + I(logArea_t0^2), 
                          data = wrangled_adu_df)
summary(gro_abs_adu_area_quad_lm)

# Adult fruiting ~ log(Area) + Height
fru_area_ht_glm <- glm(Fruit ~ logArea + Height, 
                       data = wrangled_by_yr_adu_df, family = binomial)
summary(fru_area_ht_glm)

# Set up model parameters from fitted models
params <- data.frame(
  juv.sur.int    = summary(sur_juv_area_glm)$coefficients[1,1],
  juv.sur.slope  = summary(sur_juv_area_glm)$coefficients[2,1],
  juv.gro.int    = summary(gro_abs_juv_area_lm)$coefficients[1,1],
  juv.gro.slope  = summary(gro_abs_juv_area_lm)$coefficients[2,1],
  juv.gro.sd     = sd(resid(gro_abs_juv_area_lm)),
  mat.int    = summary(mat_area_ht_int_glm)$coefficients[1,1],
  mat.area.slope    = summary(mat_area_ht_int_glm)$coefficients[2,1],
  mat.ht.slope      = summary(mat_area_ht_int_glm)$coefficients[3,1],
  mat.area.ht.slope = summary(mat_area_ht_int_glm)$coefficients[4,1],
  adu.sur.int           = summary(sur_adu_area_ht_int_glm)$coefficients[1,1],
  adu.sur.area.slope    = summary(sur_adu_area_ht_int_glm)$coefficients[2,1],
  adu.sur.ht.slope      = summary(sur_adu_area_ht_int_glm)$coefficients[3,1],
  adu.sur.area.ht.slope = summary(sur_adu_area_ht_int_glm)$coefficients[4,1],
  adu.gro.int    = summary(gro_abs_adu_area_quad_lm)$coefficients[1,1],
  adu.gro.slope  = summary(gro_abs_adu_area_quad_lm)$coefficients[2,1],
  adu.gro.slope2 = summary(gro_abs_adu_area_quad_lm)$coefficients[3,1],
  adu.gro.sd     = sd(resid(gro_abs_adu_area_quad_lm)),
  fru.int        = summary(fru_area_ht_glm)$coefficients[1,1],
  fru.area.slope = summary(fru_area_ht_glm)$coefficients[2,1],
  fru.ht.slope   = summary(fru_area_ht_glm)$coefficients[3,1],
  k.berries      = 0.2, 
  s0             = get.s0(0.2),  
  s1             = get.s0(0.2),  
  s2             = get.s0(0.2),
  rec.area.mean  = mean(rec_only_df$logArea),
  rec.area.sd    = sd(rec_only_df$logArea),
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
  # Juvenile survival function
  sur_juv_fun <- function(z){
    u = exp(params$juv.sur.int + params$juv.sur.slope * z)
    return(u/(1+u))
  }
  
  # Juvenile growth
  gro_juv_fun <- function(z1,z) { 			
    dnorm(z1, mean = params$juv.gro.int + params$juv.gro.slope * z, 
          sd = params$juv.gro.sd)
  }
  
  # Maturation
  mat_fun <- function(z, h) {
    u = exp(params$mat.int + params$mat.area.slope * z + params$mat.ht.slope * h + params$mat.area.ht.slope * z * h)
    return(u/(1+u))
  }
  
  # Adult survival
  sur_adu_fun <- function(z, h) {
    u = exp(params$adu.sur.int + params$adu.sur.area.slope * z + params$adu.sur.ht.slope * h + params$adu.sur.area.ht.slope * z * h)
    return(u/(1+u))
  }
  
  # Adult growth
  gro_adu_fun <- function(z1, z) {
    dnorm(z1, mean = params$adu.gro.int + params$adu.gro.slope * z + params$adu.gro.slope2 * z * z, 
          sd = params$adu.gro.sd)
  }
  
  # Fruiting
  fru_fun <- function(z, h) {
    u = exp(params$fru.int + params$fru.area.slope * z + params$fru.ht.slope * h)
    return(u/(1+u))
  }
  
  # Fecundity (berries produced per fruiting individual)
  fec_fun <- function(z) {
    params$k.berries * exp(z)
  }
  
  # 1-year-old seedling heights
  S1_ht_fun <- function(h1){
    dnorm(h1, mean = params$rec.ht.mean, sd = params$rec.ht.sd)
  }
  
  # Juvenile sizes
  J_z_fun <-function(z1){
    dnorm(z1, mean = params$rec.area.mean, sd = params$rec.area.sd)
  }

  ## Construct sub-kernels
  ### JJ sub-kernel
  # Describes probability of a juvenile of size z and height h becoming a juvenile of size z1 and height h1
  
  # JJ sub-kernel - juv sur * juv gro * (1 - mat)
  JJ_zh <- function(z1, h1, z, h) {
    return(sur_juv_fun(z) * gro_juv_fun(z1, z) * (1 - mat_fun(z, h)) * ht_fun(h1, h))
  }
  
  # Number of meshpoints
  meshz <- mesh
  meshh <- mesh
  
  # Set limits for integration
  Ljz <- 0.99 * min(wrangled_by_yr_df$logArea, na.rm = TRUE)
  Ujz <- min_size_rep
  Lh  <- 0.99 * min(wrangled_by_yr_df$Height, na.rm = TRUE)
  Uh  <- 1.01 * max(wrangled_by_yr_df$Height, na.rm = TRUE)
  
  # Set bin size and midpoints for each trait variable
  
  # Log area bin size for juveniles
  hjz <- (Ujz - Ljz) / meshz
  
  # Log area midpoints for juveniles
  yjz <- Ljz + hjz * ((1:meshz) - 0.5)
  
  # Height bin size
  hh <- (Uh - Lh) / meshh
  
  # Height midpoints
  yh <- Lh + hh * ((1:meshh) - 0.5)
  
  # Create function to fun through mesh points
  eta_ij <- function(i, j, meshz) {(j - 1) * meshz + i}
  
  # Create matrix of indices with which to evaluate mega matrix
  Eta <- outer(1:meshz, 1:meshh, eta_ij, meshz = meshz); 
  
  # Create mega matrix of 0s in which to put transition probabilities at each area-height combination
  JJ <- matrix(0, meshz*meshh, meshz*meshh)
  
  # Create array of 0s in which to transition probabilities at each area-height combination
  JJvals <- array(0, c(meshz, meshh, meshz, meshh))
  
  # Run through each set of meshpoints and evaluate transition probabilities
  for(i in 1:meshz){
    for(j in 1:meshh){
      for(k in 1:meshh){
        jjvals = JJ_zh(yjz, yh[k], yjz[i], yh[j])
        JJ[Eta[,k], Eta[i,j]] = jjvals
        JJvals[,k,i,j] = jjvals
      }
    }
  }
  
  # Multiply matrix by bin sizes to standardise transition probabilities
  JJ <- hjz * hh * JJ
  
  # Plot JJ kernel at median height and area
  image(t(log(JJ[(meshz*meshh/2+1):(meshz*meshh/2+meshz), (meshz*meshh/2+1):(meshz*meshh/2+meshz)] + 0.01)))

  ### JA sub-kernel
  # Describes probability of a juvenile of size z and height h becoming an adult of size z1 and height h1
  
  # JA sub-kernel - juv sur * juv gro * mat
  JA_zh <- function(z1, h1, z, h) {
    return(sur_juv_fun(z) * gro_juv_fun(z1, z) * mat_fun(z, h) * ht_fun(h1, h))
  }
  
  # Same number of meshpoints (meshz, meshh) as before
  
  # Set limits for integration - adults
  Laz <- min_size_rep
  Uaz <- 1.01 * max(wrangled_by_yr_df$logArea, na.rm = TRUE)
  
  # Set bin size and midpoints for each trait variable - juvenile size and height as before
  
  # Log area bin size for adults
  haz <- (Uaz - Laz) / meshz
  
  # Log area midpoints for adults
  yaz <- Laz + haz * ((1:meshz) - 0.5)
  
  # Create matrix of indices with which to evaluate mega matrix
  Eta <- outer(1:meshz, 1:meshh, eta_ij, meshz = meshz); 
  
  # Create mega matrix of 0s in which to put transition probabilities at each area-height combination
  JA <- matrix(0, meshz*meshh, meshz*meshh)
  
  # Create array of 0s in which to transition probabilities at each area-height combination
  JAvals <- array(0, c(meshz, meshh, meshz, meshh))
  
  # Run through each set of meshpoints and evaluate transition probabilities
  for(i in 1:meshz){
    for(j in 1:meshh){
      for(k in 1:meshh){
        javals = JA_zh(yaz, yh[k], yjz[i], yh[j])
        JA[Eta[,k], Eta[i,j]] = javals
        JAvals[,k,i,j] = javals
      }
    }
  }
  
  # Multiply matrix by bin sizes to standardise transition probabilities
  JA <- haz * hh * JA
  
  # Plot JA kernel at median height and area
  image(t(log(JA[(meshz*meshh/2+1):(meshz*meshh/2+meshz), (meshz*meshh/2+1):(meshz*meshh/2+meshz)] + 0.01)))
  
  ### AA sub-kernel
  # Describes probability of an adult of size z and height h becoming an adult of size z1 and height h1
  # AA sub-kernel - adu sur * adu gro
  AA_zh <- function(z1, h1, z, h) {
    return(sur_adu_fun(z, h) * gro_adu_fun(z1, z) * ht_fun(h1, h))
  }
  
  # Same meshpoints, midpoints and bin widths as before for adults
  # Create matrix of indices with which to evaluate mega matrix
  Eta <- outer(1:meshz, 1:meshh, eta_ij, meshz = meshz); 
  
  # Create mega matrix of 0s in which to put transition probabilities at each area-height combination
  AA <- matrix(0, meshz*meshh, meshz*meshh)
  
  # Create array of 0s in which to transition probabilities at each area-height combination
  AAvals <- array(0, c(meshz, meshh, meshz, meshh))
  
  # Run through each set of meshpoints and evaluate transition probabilities
  for(i in 1:meshz){
    for(j in 1:meshh){
      for(k in 1:meshh){
        aavals = AA_zh(yaz, yh[k], yaz[i], yh[j])
        AA[Eta[,k], Eta[i,j]] = aavals
        AAvals[,k,i,j] = aavals
      }
    }
  }
  
  # Multiply matrix by bin sizes to standardise transition probabilities
  AA <- haz * hh * AA
  
  # Plot JA kernel at median height and area
  image(t(log(AA[(meshz*meshh/2+1):(meshz*meshh/2+meshz), (meshz*meshh/2+1):(meshz*meshh/2+meshz)] + 0.01)))

  ### AS1 sub-kernel
  # Describes number of S1 from an adult of size z and height h establishing at height h1
  # AS1 sub-kernel - fru prob * fec * estab prob * new ht
  AS1_zh <- function(h1, z, h) {
    return(fru_fun(z, h) * fec_fun(z) * params$s0 * S1_ht_fun(h1))
  }
  
  # Create function to fun through mesh points
  eta_ij <- function(i, j, meshz) {(j - 1) * meshz + i}
  
  # Modify kernel sizes as z1 is discrete
  # Create mega matrix of 0s in which to put transition probabilities at each area-height combination
  AS1 <- matrix(0, meshh, meshz * meshh)
  
  # Create array of 0s in which to transition probabilities at each area-height combination
  AS1vals <- array(0, c(meshh, meshz, meshh))
  
  # Run through each set of meshpoints and evaluate transition probabilities
  for(i in 1:meshz){
    for(j in 1:meshh){
      for(k in 1:meshh){
        as1vals = AS1_zh(yh[k], yaz[i], yh[j])
        AS1[k, eta_ij(i, j, meshz)] = as1vals
        AS1vals[i,k,j] = as1vals
      }
    }
  }
  
  # Multiply matrix by bin sizes to standardise transition probabilities
  AS1 <- hh * AS1
  
  # Plot AS1 kernel at lowest height
  image(t(log(AS1[, 1:meshz] + 0.01)))

  # S1S2 sub-kernel
  # Describes probability of S1 surviving to become S2
  S1S2_zh <- function(h1, h) {
    return(params$s1 * ht_fun(h1, h))
  }
  
  # Modify kernel sizes as z1 and z are discrete
  # Create mega matrix of 0s in which to put transition probabilities at each area-height combination
  S1S2 <- matrix(0, meshh, meshh)
  
  # Create array of 0s in which to transition probabilities at each area-height combination
  S1S2vals <- array(0, c(meshh, meshh))
  
  # Run through each set of meshpoints and evaluate transition probabilities
  for(j in 1:meshh){
    for(k in 1:meshh){
        s1s2vals = S1S2_zh(yh[k], yh[j])
        S1S2[k, j] = s1s2vals
        S1S2vals[k, j] = s1s2vals
     }
  }
  
  # Multiply matrix by bin sizes to standardise transition probabilities - TRIPLE CHECK THIS IS CORRECT; bin size of z is 1
  S1S2 <- hh * S1S2
  
  # Plot S2S1 kernel
  image(t(log(S1S2 + 0.01)))
  
  ### S2J sub-kernel
  # Describes probability of an S2 establishing as a juvenile of size z1 at height h1
  # S2J sub-kernel - fru prob * fec * estab prob * new ht
  S2J_zh <- function(z1, h1, h) {
    return(params$s2 * J_z_fun(z1) * ht_fun(h1, h))
  }
  
  # Create function to fun through mesh points
  eta_ij <- function(i, j, meshz) {(j - 1) * meshz + i}
  
  # Modify kernel sizes as z is discrete
  # Create mega matrix of 0s in which to put transition probabilities at each area-height combination
  S2J <- matrix(0, meshz * meshh, meshh)
  
  # Create array of 0s in which to transition probabilities at each area-height combination
  S2Jvals <- array(0, c(meshz, meshh, meshh))
  
  # Run through each set of meshpoints and evaluate transition probabilities
  for(i in 1:meshz){
    for(j in 1:meshh){
      for(k in 1:meshh){
        s2jvals = S2J_zh(yjz[i], yh[j], yh[k])
        S2J[eta_ij(i, j, meshz), k] = s2jvals
        S2Jvals[i,k,j] = s2jvals
      }
    }
  }
  
  # Multiply matrix by bin sizes to standardise transition probabilities
  S2J <- hjz * hh * S2J
  
  # Plot S2J kernel at median height
  image(t(log(S2J[1:meshh, ] + 0.01)))
  
  # Remove vals to clear space
  JJvals <- 0
  JAvals <- 0
  AAvals <- 0
  AS1vals <- 0
  S1S2vals <- 0
  S2Jvals <- 0

  # Combine sub-kernels into K kernel
  # Construct P kernel
  # Combine A-row matrices with structural 0 matrices
  PA_row <- cbind(matrix(0, nrow = meshz * meshh, ncol = meshh),
                 matrix(0, nrow = meshz * meshh, ncol = meshh),
                 JA,
                 AA)
  
  # Combine J-row matrices with structural 0 matrices - leave S2J for F kernel
  PJ_row <- cbind(matrix(0, nrow = meshz * meshh, ncol = meshh),
                 S2J,
                 JJ,
                 matrix(0, nrow = meshz * meshh, ncol = meshz * meshh))
  
  # Combine S2-row matrices with structural 0 matrices - leave S2J for F kernel
  PS2_row <- cbind(S1S2, 
                  matrix(0, nrow = meshh, ncol = meshh),
                  matrix(0, nrow = meshh, ncol = meshz * meshh),
                  matrix(0, nrow = meshh, ncol = meshz * meshh))
  
  # Combine S1 row matrices with structural 0 matrices - leave AS1 for F kernel
  PS1_row <- cbind(matrix(0, nrow = meshh, ncol = meshh),
                  matrix(0, nrow = meshh, ncol = meshh),
                  matrix(0, nrow = meshh, ncol = meshz * meshh),
                  0*AS1)
  
  # Combine rows into P kernel
  P_kernel <- rbind(PA_row, PJ_row, PS2_row, PS1_row)
  
  # Construct F kernel
  # Combine A-row matrices with structural 0 matrices
  FA_row <- cbind(matrix(0, nrow = meshz * meshh, ncol = meshh),
                 matrix(0, nrow = meshz * meshh, ncol = meshh),
                 0*JA,
                 0*AA)
  
  # Combine J-row matrices with structural 0 matrices
  FJ_row <- cbind(matrix(0, nrow = meshz * meshh, ncol = meshh),
                 0*S2J,
                 0*JJ,
                 matrix(0, nrow = meshz * meshh, ncol = meshz * meshh))
  
  # Combine S2-row matrices with structural 0 matrices
  FS2_row <- cbind(0*S1S2, 
                  matrix(0, nrow = meshh, ncol = meshh),
                  matrix(0, nrow = meshh, ncol = meshz * meshh),
                  matrix(0, nrow = meshh, ncol = meshz * meshh))
  
  # Combine S1 row matrices with structural 0 matrices
  FS1_row <- cbind(matrix(0, nrow = meshh, ncol = meshh),
                  matrix(0, nrow = meshh, ncol = meshh),
                  matrix(0, nrow = meshh, ncol = meshz * meshh),
                  AS1)
  
  # Combine rows into F kernel
  F_kernel <- rbind(FA_row, FJ_row, FS2_row, FS1_row)

  # Remove rows to clear space
  PA_row <- 0
  PJ_row <- 0
  FJ_row <- 0
  FA_row <- 0

  # Combine P and F into mega matrix kernel K
  K <- P_kernel + F_kernel
  
  ## Plot IPM
  # Plot IPM at reference height
  
  # Find bin of reference height
  ref_ht_bin <- as.integer((ref_ht - Lh)/(Uh - Lh) * meshh)
  
  # Construct K kernel at mean height
  # Divide by maximum value for plotting
  JA_plot <- JA/max(JA)
  AA_plot <- AA/max(AA)
  S2J_plot <- S2J/max(S2J)
  JJ_plot <- JJ/max(JJ)
  S1S2_plot <- S1S2/max(S1S2)
  AS1_plot <- AS1/max(AS1)
  
  # Combine A-row matrices with structural 0 matrices
  plot_A_row <- cbind(matrix(0, nrow = meshz, ncol = 1),
                  matrix(0, nrow = meshz, ncol = 1),
                  JA_plot[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)],
                  AA_plot[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)])
  
  # Combine J-row matrices with structural 0 matrices
  plot_J_row <- cbind(matrix(0, nrow = meshz, ncol = 1),
                  S2J_plot[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ref_ht_bin],
                  JJ_plot[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)],
                  matrix(0, nrow = meshz, ncol = meshz))
  
  # Combine S2-row matrices with structural 0 matrices
  plot_S2_row <- cbind(S1S2_plot[ref_ht_bin, ref_ht_bin], 
                   matrix(0, nrow = 1, ncol = 1),
                   matrix(0, nrow = 1, ncol = meshz),
                   matrix(0, nrow = 1, ncol = meshz))
  
  # Combine S1 row matrices with structural 0 matrices
  plot_S1_row <- cbind(matrix(0, nrow = 1, ncol = 1),
                   matrix(0, nrow = 1, ncol = 1),
                   matrix(0, nrow = 1, ncol = meshz),
                   matrix(AS1_plot[ref_ht_bin, ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)], nrow = 1))
  
  
  # Combine rows into plottable kernel
  plot_kernel <- rbind(plot_S1_row, plot_S2_row, plot_J_row, plot_A_row)
  
  # Plot kernel
  image(t(plot_kernel))
  
  # Create P kernel at reference height
  # Combine A-row matrices with structural 0 matrices
  ref_PA_row <- cbind(matrix(0, nrow = meshz, ncol = 1),
                      matrix(0, nrow = meshz, ncol = 1),
                      JA[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)],
                      AA[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)])
  
  # Combine J-row matrices with structural 0 matrices
  ref_PJ_row <- cbind(matrix(0, nrow = meshz, ncol = 1),
                      S2J[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ref_ht_bin],
                      JJ[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)],
                      matrix(0, nrow = meshz, ncol = meshz))
  
  # Combine S2-row matrices with structural 0 matrices
  ref_PS2_row <- cbind(S1S2[ref_ht_bin, ref_ht_bin], 
                       matrix(0, nrow = 1, ncol = 1),
                       matrix(0, nrow = 1, ncol = meshz),
                       matrix(0, nrow = 1, ncol = meshz))
  
  # Combine S1 row matrices with structural 0 matrices
  ref_PS1_row <- cbind(matrix(0, nrow = 1, ncol = 1),
                       matrix(0, nrow = 1, ncol = 1),
                       matrix(0, nrow = 1, ncol = meshz),
                       0*matrix(AS1[ref_ht_bin, ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)], nrow = 1))
  
  ref_P_kernel <- rbind(ref_PS1_row, ref_PS2_row, ref_PJ_row, ref_PA_row)
  
  # Create F kernel at mean height
  # Combine A-row matrices with structural 0 matrices
  ref_FA_row <- cbind(matrix(0, nrow = meshz, ncol = 1),
                       matrix(0, nrow = meshz, ncol = 1),
                       0*JA[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)],
                       0*AA[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)])
  
  # Combine J-row matrices with structural 0 matrices
  ref_FJ_row <- cbind(matrix(0, nrow = meshz, ncol = 1),
                       0*S2J[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ref_ht_bin],
                       0*JJ[((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh), ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)],
                       matrix(0, nrow = meshz, ncol = meshz))
  
  # Combine S2-row matrices with structural 0 matrices
  ref_FS2_row <- cbind(0*S1S2[ref_ht_bin, ref_ht_bin], 
                        matrix(0, nrow = 1, ncol = 1),
                        matrix(0, nrow = 1, ncol = meshz),
                        matrix(0, nrow = 1, ncol = meshz))
  
  # Combine S1 row matrices with structural 0 matrices
  ref_FS1_row <- cbind(matrix(0, nrow = 1, ncol = 1),
                        matrix(0, nrow = 1, ncol = 1),
                        matrix(0, nrow = 1, ncol = meshz),
                        matrix(AS1[ref_ht_bin, ((ref_ht_bin-1)*meshh + 1):(ref_ht_bin*meshh)], nrow = 1))
  
  ref_F_kernel <- rbind(ref_FS1_row, ref_FS2_row, ref_FJ_row, ref_FA_row)
  
  # Plot full kernel
  ref_kernel <- ref_P_kernel + ref_F_kernel
  image(t(ref_kernel))
  
  # Output P kernel
  output_ls[[1]] <- ref_P_kernel
  # Output F kernel
  output_ls[[2]] <- ref_F_kernel
  # Output full kernel
  output_ls[[3]] <- ref_kernel
  # Output meshes
  output_ls[[4]] <- mesh
  
  # Rename outputs
  names(output_ls) <- c("P", "F", "K", "mesh")
  
  # Return outputs
  return(output_ls)
}

# Build initial IPM with 100 meshpoints and reference height as mean height
mean_ht <- params$rec.ht.mean
outputs <- build.ipm(params = params, mesh = 50, ref_ht = mean_ht)
lambda(outputs[["K"]])

# Extract lambda and life history traits from IPM
get.traits <- function(outputs){
  # Specify first non-propagule (non-seedling) stage
  notProp <- 3
  
  # Calculate lambda
  outputs[[5]] <- eigen.analysis(outputs[["K"]])$lambda1
  names(outputs)[5] <- "Lambda"
  
  # Calculate R0 from P and F sub-kernels, starting at first non-propagule stage
  outputs[[6]] <- net_repro_rate(outputs[["P"]], outputs[["F"]], notProp)
  names(outputs)[6] <- "R0"
  
  # Calculate generation time from R0 and lambda
  outputs[[7]] <- abs(log(outputs[["R0"]])/log(outputs[["Lambda"]]))
  names(outputs)[7] <- "GenTfun"
  
  # Calculate mean life expectancy from mature distribution
  # Identify which stages are reproductive
  mat_dist <- mature_distrib(outputs[["P"]], start = notProp, repro_stages = repro_stages(outputs[["F"]]))
  outputs[[8]] <- life_expect_mean(matU = outputs[["P"]], mixdist = mat_dist, start = NULL)
  names(outputs)[8] <- "Lmean"
  
  # Calculate mean age at maturity, from first non-propagule stage
  outputs[[9]] <- mature_age(matU = outputs[["P"]], matF = outputs[["F"]], start = notProp)
  names(outputs)[9] <- "La"
  
  # Calculate mean reproductive window
  outputs[[10]] <- outputs[["Lmean"]] - outputs[["La"]]
  names(outputs)[10] <- "Lamean"

  # Return outputs
  return(outputs)
}

for(lx_crit in lx_crits){}

# Extract initial life history traits
traits <- get.traits(outputs)

# Fix establishment constant s0 -------------------------------------------
# Create lists to store lambdas and s0 values
lambdas <- list()
s0s <- list()

# Begin with initial lambda
lambdas[[1]] <- traits[["Lambda"]]

# Start with s0 derived from k = 0.2
s0s[[1]] <- get.s0(params$k.berries)

# Initialise loop, set number of iterations and set progress bar
i <- 1
iter <- 10
pb <- txtProgressBar(min = 0, max = iter, style = 3)

# If lambda > 1.11, decrease s0; if lambda < 1.09, increase s0 
while(abs(1.1 - as.numeric(lambdas[[i]])) > 0.01 & i < iter) {
  if(lambdas[[i]] > 1.1){
    params$s0 <- 0.9 * s0s[[i]]
    params$s1 <- params$s0
    params$s2 <- params$s0
  } else {
    params$s0 <- 1.1 * s0s[[i]]
    params$s1 <- params$s0
    params$s2 <- params$s0
  }
  s0s[[i+1]] <- params$s0
  outputs <- build.ipm(params = params, mesh = 50, ref_ht = params$rec.ht.mean)
  names(outputs) <- c("P", "F", "K", "mesh")
  lambdas[[i+1]] <- lambda(outputs[["K"]])
  setTxtProgressBar(pb, i)
  i <- i + 1
}
close(pb)

# # Create lists to store longevities and s0 values
# longs <- list()
# s0s <- list()
# 
# # Begin with initial longevity
# longs[[1]] <- longevity(outputs[["P"]], 
#                         start = mature_distrib(outputs[["P"]], start = 3, repro_stages = repro_stages(outputs[["F"]])))
# 
# # Start with s0 derived from k = 0.2
# s0s[[1]] <- get.s0(params$k.berries)
# 
# # Initialise loop, set number of iterations and set progress bar
# i <- 1
# iter <- 10
# pb <- txtProgressBar(min = 0, max = iter, style = 3)
# # If longevity > 31, decrease s0; if longevity < 29, increase s0 
# while(abs(30 - as.numeric(longs[[i]])) > 1 & i < iter) {
#   if(longs[[i]] > 31){
#     params$s0 <- 0.9 * s0s[[i]]
#     params$s1 <- params$s0
#     params$s2 <- params$s0
#   } else {
#     params$s0 <- 1.1 * s0s[[i]]
#     params$s1 <- params$s0
#     params$s2 <- params$s0
#   }
#   s0s[[i+1]] <- params$s0
#   outputs <- build.ipm(params = params, mesh = 50, ref_ht = params$rec.ht.mean)
#   names(outputs) <- c("P", "F", "K", "mesh")
#   longs[[i+1]] <- longevity(outputs[["P"]], start = mature_distrib(outputs[["P"]], start = 3, repro_stages = repro_stages(outputs[["F"]])))
#   setTxtProgressBar(pb, i)
#   i <- i + 1
# }
# close(pb)

# Set s0 as iteration output
params$s0 <- s0s[[length(s0s)]]

# Build final IPM ---------------------------------------------------------
# Build and plot final IPM with 50 meshpoints
outputs <- build.ipm(params = params, mesh = 50, ref_ht = params$rec.ht.mean)

# Extract final life history traits
traits <- get.traits(outputs)

# Extract output LHTs -----------------------------------------------------
# Create data frame for mistletoe LHTs
Mst_LHTs <- data.frame(SpeciesAccepted = "Viscum album",
                       GenTfun = traits[["GenTfun"]],
                       Lmean = traits[["Lmean"]],
                       La = traits[["La"]],
                       Lamean = traits[["Lamean"]])

# Export .csv for mistletoe LHTs
write.csv(Mst_LHTs, "Mst_LHTs.csv")

# Sensitivity of λ to mesh points ----------------------------------------------
# Set mesh sizes
meshes <- seq(30, 100, by = 10)

# Create list to store lambda values
lambdas <- list()

# Extract lambda for each mesh value
for(mesh in meshes){
  outputs <- build.ipm(params = params, mesh = mesh, ref_ht = params$rec.ht.mean)
  lambdas[[mesh/5 - 1]] <- lambda(outputs[["K"]])
}

# Plot lambda against number of mesh points
lambda_mesh_df <- data.frame(mesh = meshes, lambda = unlist(lambdas))
Fig_S6 <- ggplot(data = lambda_mesh_df, aes(x = (mesh*2+2), y = lambda)) +
  labs(x = "Total meshpoints", y = "λ") +
          geom_line() +
          theme_bw()
Fig_S6

# Sensitivity of λ to parameters (including k) ---------------------------------
lambdas <- list()
for(i in 1:ncol(params)){
  params[,i] <- params[,i] + 0.001
  outputs <- build.ipm(params = params, mesh = mesh, ref_ht = params$rec.mean.ht)
  lambdas[[i]] <- lambda(outputs[["K"]])
  params[,i] <- params[,i] - 0.001
}

# Plot IPM at different heights -------------------------------------------
# Specify heights at which to plot IPM
hts <- c(params$rec.ht.mean - 2 * params$rec.ht.sd,
         params$rec.ht.mean - params$rec.ht.sd,
         params$rec.ht.mean,
         params$rec.ht.mean + params$rec.ht.sd,
         params$rec.ht.mean + 2 * params$rec.ht.sd)

# Plot IPMs
for(ht in hts){
  outputs <- build.ipm(params = params, mesh = mesh, ref_ht = ht)
}
