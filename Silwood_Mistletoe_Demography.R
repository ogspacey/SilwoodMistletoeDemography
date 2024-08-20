### Silwood Mistletoe Demography
### Oliver G. Spacey, Owen R. Jones, Sydne Record, Roberto Salguero-Gómez
### 31.07.2024

# Pre-amble ---------------------------------------------------------------
# Clear environment
rm(list = ls())

# Load packages
library(ggplot2)   # Data visualisation
library(tidyverse) # Data wrangling

# Set working directory as required

# Load data ---------------------------------------------------------------
# Load mistletoe demographic data
mst_df_raw <- read.csv("mistletoe_data_raw.csv")

# Load host metadata
hst_df_raw <- read.csv("host_metadata_raw.csv")

# Wrangle data ------------------------------------------------------------
# Remove bottom rows with NAs?

# Remove hosts from host data frame which mistletoes were not measured in
hst_df_raw <- filter(hst_df_raw, Host_ID %in% unique(mst_df_raw$Host_ID))

# Create genus variable for hosts
hst_df_raw <- mutate(hst_df_raw, Genus = word(Genus_spp))

# Combine mistletoe demographic data and host data
mst_hst_df <- full_join(mst_df_raw, hst_df_raw, by = "Host_ID")

# Create individual ID combining mistletoe and host IDs
mst_hst_df <- mutate(mst_hst_df, Indiv_ID = paste(Mistletoe_ID, Host_ID, sep = "_"))

# Remove columns measuring perimeter as we will only use area as an indicator of size
mst_hst_df <- select(mst_hst_df, -contains('Perim'))

# Remove mistletoes that are never seen
mst_hst_df <- filter(mst_hst_df, Notes != "Never seen")

# DEAL WITH MULTI-CLUMP


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

# Rectify size based on horizontal distance to base of tree (standard)
# If an object is twice as far away relative to the standard, it will appear twice as small in length, and four times as small in area
# Calculate the area adjusting for distance relative to the standard (horizontal distance) as Adjusted area = (Distance to object/Distance to standard)^2 * Measured area

# Select columns of areas, hypotenuses and horizontal distance to standard (host), and Indiv_ID
adjusted_df <- select(mst_hst_df, Indiv_ID, starts_with("Area"), Combined_hypotenuse, Host_horizontal_distance)
adjusted_df[,2:11] <- adjusted_df[,2:11] * (adjusted_df$Combined_hypotenuse/adjusted_df$Host_horizontal_distance)^2

# Stack area columns
stacked_areas_df <- pivot_longer(adjusted_df, cols = starts_with("Area"), names_to = "Year", values_to = "Area")

# Visualise areas in histogram
ggplot(data = stacked_areas_df, aes(x = Area)) +
  geom_histogram() +
  theme_bw()
# Highly skewed distribution

# Log-transform to reduce and re-plot areas, also examining distribution of years
ggplot(data = stacked_areas_df, aes(x = log(Area), fill = Year)) +
  geom_histogram() +
  theme_bw()
# More normally distributed
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
years <- c(14:22) # vector of years to run for loop through
long_df_ls <- list()   # list of long data frames to store each year
for(year in years){ # run through years sequentially
  Status_t0  <- paste("Status", year, sep = "") 
  logArea_t0 <- paste("logArea", year, sep = "")
  Fruit_t0   <- paste("Fruit", year, sep = "")
  Status_t1  <- paste("Status", year+1, sep = "")
  logArea_t1 <- paste("logArea", year+1, sep = "")
  Fruit_t1   <- paste("Fruit", year+1, sep = "")
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
years <- c(14:23) # vector of years to run for loop through
long_df_ls <- list()   # list of long data frames to store each year
for(year in years){ # run through years sequentially
  Status  <- paste("Status", year, sep = "") 
  logArea <- paste("logArea", year, sep = "")
  Fruit   <- paste("Fruit", year, sep = "")
  mst_df_year <- data.frame(Host_ID = mst_hst_df$Host_ID,
                            Mistletoe_ID = mst_hst_df$Mistletoe_ID,
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

## Calculate parasite loads ------------------------------------------------

# Sum parasite load for each host each year
# Create list of host IDs
hosts <- hst_df_raw$Host_ID

# Create empty dataframe
par_lds_df <- data.frame(Host_ID = hst_df_raw$Host_ID)
row.names(par_lds_df) <- hst_df_raw$Host_ID
par_lds_df[c("PL_14", "PL_15", "PL_16", "PL_17", "PL_18",
             "PL_19", "PL_20", "PL_21", "PL_22", "PL_23")] <- NA
years <- as.character(14:23)
for(year in years){
  for(host in hosts){
  filtered_df <- filter(mst_hst_df, Host_ID == host) %>%
                 select(ends_with(year))
  parasite_load <- length(which(filtered_df[1] == "FC" | filtered_df[1] == "Surv"))
  pl_col_name <- paste("PL_", year, sep = "")
  par_lds_df[host, pl_col_name] <- parasite_load
  }
}

# Calculate how many mistletoes seen each year
colSums(par_lds_df[,-1])

# Calculate maximum parasite load observed for each tree
par_lds_df$Max_PL <- pmax(par_lds_df$PL_14,
                          par_lds_df$PL_15,
                          par_lds_df$PL_16,
                          par_lds_df$PL_17,
                          par_lds_df$PL_18,
                          par_lds_df$PL_19,
                          par_lds_df$PL_20,
                          par_lds_df$PL_21,
                          par_lds_df$PL_22,
                          par_lds_df$PL_23)

# Add host metadata to PL data frame
par_lds_df <- full_join(par_lds_df, hst_df_raw)

## Assign sex ----------------------------------------------------
# Assign known females 
sex_df <- mst_hst_df %>%
  select(Indiv_ID, starts_with("Fruit"))
sex_df <- mutate(sex_df, 
                 Known_sex = case_when(rowSums(sex_df[-1], na.rm = TRUE) >= 1 ~ "Female",
                                         TRUE ~ "Unknown"))

# Count number of times each sampled for fruiting
sex_df <- mutate(sex_df, 
                 N_fru_sampled = rowSums(!is.na(sex_df[2:11])))
  
# Join Known_sex to long data frame
long_sex_df <- long_by_yr_df %>% 
                 full_join(sex_df[,-c(2:11)], by = "Indiv_ID")

# Assign individuals below minimum size at reproduction as juveniles
# Calculate minimum size at reproduction
min_size_rep <- min(filter(long_sex_df, Fruit == 1)$logArea, na.rm = TRUE)

# Assign juveniles, keep known females and leave others as unknown
long_sex_df$Known_sex <- case_when(long_sex_df$Known_sex == "Female" ~ "Female",
                                   long_sex_df$logArea < min_size_rep ~ "Juvenile",
                                   .default = "Unknown")

# Create data frame of just adults in which to assign sex
adults_df <- filter(long_sex_df, Known_sex != "Juvenile")

# For each unknown adult, we want to estimate the probability that it is female (Fem) given that it is not fruiting (Fru'), P(Fem|Fru')
# P(Fem|Fru') = P(Fru'|Fem) * P(Fem) / P(Fru'), therefore must calculate P(Fru'|Fem), P(Fem) and P(Fru')
# P(Fru'|Fem) = 1 - P(Fru|Fem), and P(Fru|Fem) may be a function of mistletoe size and height

# Calculate probability that a female fruits in a given year dependent on size
# Filter long data frame for females
known_fem_df <- filter(adults_df, Known_sex == "Female")

# Plot fruiting probability of females based on size
ggplot(data = known_fem_df, aes(x = logArea, y = Fruit)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()
# Larger females more likely to fruit

# Plot fruiting probability of females based on height
ggplot(data = known_fem_df, aes(x = Height, y = Fruit)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()
# Fruits more likely lower down - possibly a biological effect or due to observation probability
  
# Model probability of fruiting for females as a function of size and height
fru_area_ht_glm <- glm(Fruit ~ logArea + Height, data = known_fem_df, family = binomial)
summary(fru_area_ht_glm)
# Female fruiting probability depends on size and height
# Assume that all females fruit with the same probability function of size and height

# Predict probability of individuals of unknown sex fruiting and not fruiting in a giving year
# if they were female, based on area and height, P(Fru|Fem)
long_uk_sex_df <- filter(long_sex_df, Known_sex == "Unknown")
long_uk_sex_df <- mutate(long_uk_sex_df, 
                         P_uk_fru = predict(fru_area_ht_glm, 
                                              newdata = data.frame(logArea = long_uk_sex_df$logArea, 
                                                                Height = long_uk_sex_df$Height), 
                                              type = "response"),
                         P_uk_no_fru = 1 - P_uk_fru) #P(Fru'|Fem) = 1 - P(Fru|Fem)

# P(Fem) = P(Fru' ∩ Fem) / P(Fru'|Fem), and P(Fem) gives the estimated (1 -) sex ratio for the population
# P(Fem) is the probability of a given adult individual being female
# P(Fru' ∩ Fem) is the probability that a given individual is both not fruiting and a female
# The number of non-fruiting known females/all adults in sample
# Filter data frame to observations where fruiting sampled, count number of rows for adult sample size
fru_sample_size <- nrow(filter(adults_df, is.na(Fruit) == FALSE))
# Count number of rows of known females which are not fruiting
no_non_fru_fem <- nrow(filter(known_fem_df, Fruit == 0))
P_no_fru_and_fem <- no_non_fru_fem/fru_sample_size
# P(Fru'|Fem) is the probability that an individual is not fruiting given it is female
# This probability varies depending on size and height
# To calculate P(Fem), we take the P(Fru|Fem) across known females sampled for mean height and size
# back-transform and subtract from 1
P_no_fru_given_fem <- as.numeric(1 - plogis(predict(fru_area_ht_glm, 
                                         newdata = data.frame(logArea = mean(known_fem_df$logArea, na.rm = TRUE), 
                                                              Height = mean(known_fem_df$Height, na.rm = TRUE), 
                                         type = "response"))))
# Calculate P(Fem) = P(Fru' ∩ Fem) / P(Fru'|Fem)
P_fem <- P_no_fru_and_fem / P_no_fru_given_fem
# Sex ratio of 0.45:0.55 M:F

# P(Fru') = probability of an individual not fruiting
# The number of non-fruiting adults/all adults in sample
no_non_fru <- nrow(filter(adults_df, Fruit == 0))
P_no_fru <- no_non_fru/fru_sample_size

# P(Fem|Fru') = P(Fru'|Fem) * P(Fem) / P(Fru'), calculate P(Fem|Fru') for each mistletoe in a given year
long_uk_sex_df <- mutate(long_uk_sex_df,
                         P_uk_fem = P_uk_no_fru * P_fem / P_no_fru)

# As a test, PREDICT PROBABILITY THAT KNOWN FEMALES WOULD BE FEMALE USING THIS METHOD

# Test whether fruiting probability one year is linked to fruiting probability in previous year for a female
# List known females
known_females <- unique(known_fem_df$Indiv_ID)

# Summarise changes in fruit for known females
fruit_t0_t1 <- select(long_df, Indiv_ID, Fruit_t0, Fruit_t1) %>%
               filter(Indiv_ID %in% known_females) %>%
               mutate(fruit_effect = case_when(Fruit_t0 == 0 & Fruit_t1 == 0 ~ "both0",
                                               Fruit_t0 == 1 & Fruit_t1 == 0 ~ "1then0",
                                               Fruit_t0 == 0 & Fruit_t1 == 1 ~ "0then1",
                                               Fruit_t0 == 1 & Fruit_t1 == 1 ~ "both1")) %>%
               filter(!is.na(fruit_effect))

# Create table of observed fruiting transitions in t0 and t1
observed_fru_mat <- matrix(c(table(fruit_t0_t1$fruit_effect)["both0"],  # Observed females not fruiting t0, not fruiting t1
                             table(fruit_t0_t1$fruit_effect)["0then1"], # Observed females not fruiting t0, fruiting t1
                             table(fruit_t0_t1$fruit_effect)["1then0"], # Observed females fruiting t0, not fruiting t1
                             table(fruit_t0_t1$fruit_effect)["both1"]), # Observed females fruiting t0, fruiting t1
                           nrow = 2, 
                           byrow = TRUE)

# Expected transitions based on how many females observed to be fruiting
no_non_fru_fem <- nrow(filter(fruit_t0_t1, Fruit_t0 == 0))
no_fru_fem <- nrow(filter(known_fem_df, Fruit == 1))
no_fem <- no_non_fru_fem + no_fru_fem
expected_fru_mat <- matrix(c(no_non_fru_fem * no_non_fru_fem / (no_fem * no_fem) * nrow(fruit_t0_t1), # Expected females not fruiting t0, not fruiting t1
                             no_non_fru_fem * no_fru_fem / (no_fem * no_fem) * nrow(fruit_t0_t1),     # Expected females not fruiting t0, fruiting t1
                             no_fru_fem * no_non_fru_fem / (no_fem * no_fem) * nrow(fruit_t0_t1),     # Expected females fruiting t0, not fruiting t1
                             no_fru_fem * no_fru_fem / (no_fem * no_fem) * nrow(fruit_t0_t1)),         # Expected females fruiting t0, fruiting t1
                           nrow = 2, 
                           byrow = TRUE)

# Perform Chi-square test for contingency
chisq.test(observed_fru_mat, p = expected_fru_mat / sum(expected_fru_mat))
# Individuals seen fruiting in one year more likely to be fruiting in next, cannot assume independent fruiting probability over years

# Assuming fruiting probability of a female fruiting is independent between years, estimate probability that 
  
# Explore data --------------------------------------------------------
# Visualise distribution of mistletoe heights
ggplot(data = mst_hst_df, aes(x = Combined_height)) +
  geom_histogram() +
  theme_bw()

# Visualise heights in histogram
ggplot(data = mst_hst_df, aes(x = Combined_height)) +
  geom_histogram() +
  theme_bw()

# Calculate heights of mistletoe relative to tree
mst_hst_df <- mutate(mst_hst_df, 
                     Relative_height = Combined_height/Host_height)

# Visualise relative heights in histogram
ggplot(data = mst_hst_df, aes(x = Relative_height)) +
  geom_histogram() +
  theme_bw()
# Some estimated at >1 due to parallax with mistletoes, making them appear higher than the top of the tree
# Use raw heights instead

# Visualise relationship between host and mistletoe heights
ggplot(data = mst_hst_df, aes(x = Host_height, y = Combined_height)) +
  geom_point() +
  theme_bw() +
  geom_smooth() +
  labs(x = "Host height", y = "Mistletoe height")
# Taller the host, higher up mistletoes generally are

# Visualise distribution of host heights and species
ggplot(data = hst_df_raw, aes(x = Host_height, fill = Genus_spp)) +
  geom_histogram()
# Highly colinear - hard to separate effects

# Visualise distribution of mistletoe heights and species
ggplot(data = mst_hst_df, aes(x = Combined_height, fill = Genus_spp)) +
  geom_histogram()
# Highly colinear - hard to separate effects

# Visualise distribution of host height, species and maximum parasite load
ggplot(data = par_lds_df, aes(x = Host_height, y = Max_PL, col = Genus)) +
  geom_point() +
  theme_bw()
# Any size tree can have small PL, but only huge trees can have large PL

# Visualise distribution of host species and maximum parasite load
ggplot(data = par_lds_df, aes(x = Max_PL, fill = Genus)) +
  geom_histogram() +
  theme_bw()
# Parasite load overlaps slightly between genera

# Visualise distribution of mistletoe areas
ggplot(data = long_by_yr_df, aes(x = logArea)) +
  geom_histogram()

# Visualise relationship between areas and height
ggplot(data = long_by_yr_df, aes(x = logArea, y = Height)) +
  geom_point() +
  theme_bw()
# Size not correlated with height


# Construct and select models --------------------------------------------------------

## Model survival ----------------------------------------------------------
# Create binary variable for survival 
long_df <- mutate(long_df, Survive = case_when(
  Status_t1 == "Dead" ~ 0,
  Status_t1 == "Surv" ~ 1
))

# Plot survival against log(Area)
ggplot(data = long_df, aes(x = logArea_t0, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial))

# Initial model for survival
log_reg_sur_1 <- glm(Survive ~ logArea_t0, data = long_df, family = binomial)
summary(log_reg_sur_1)

## Model growth ------------------------------------------------------------
# Plot growth in time t1 against growth in time t0
ggplot(data = long_df, aes(x = logArea_t0, y = logArea_t1)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  lims(x = c(0, max(long_df$logArea_t0)), y = c(0, max(long_df$logArea_t1))) +
  geom_smooth() +
  theme_bw()

# Select model for growth

## Model fruiting ------------------------------------------------------

# Select model for fruiting
ggplot(data = long_by_yr_df, aes(x = logArea, y = Fruit)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial))

# Initial model for 
log_reg_fru_1 <- glm(Fruit ~ logArea, data = long_by_yr_df, family = binomial)
summary(log_reg_fru_1)

# Construct IPM -----------------------------------------------------------


# Analyse IPM -------------------------------------------------------------
# Calculate asymptotic growth rate


# Calculate generation time




