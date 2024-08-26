### Silwood Mistletoe Demography
### Oliver G. Spacey, Owen R. Jones, Sydne Record, Roberto Salguero-Gómez
### 31.07.2024

# Pre-amble ---------------------------------------------------------------
# Clear environment
rm(list = ls())

# Load packages
library(ggplot2)   # Data visualisation
library(tidyverse) # Data wrangling
library(lme4)      # Regression models

# Set working directory as required

# Load data ---------------------------------------------------------------
# Load mistletoe demographic data
mst_df_raw <- read.csv("mistletoe_data_raw.csv")

# Load host metadata
hst_df_raw <- read.csv("host_metadata_raw.csv")

# Wrangle data ------------------------------------------------------------
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

# Remove record from clump containing multiple individual (ID is longer than )
mst_hst_df <- filter(mst_hst_df, nchar(Mistletoe_ID) == 4)

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
colnames(observed_fru_mat) <- c("No_Fruit_t1", "Fruit_t1")
rownames(observed_fru_mat) <- c("No_Fruit_t0", "Fruit_t0")

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
colnames(expected_fru_mat) <- c("No_Fruit_t1", "Fruit_t1")
rownames(expected_fru_mat) <- c("No_Fruit_t0", "Fruit_t0")

# Perform Chi-square test for contingency
chisq.test(observed_fru_mat, p = expected_fru_mat / sum(expected_fru_mat))
# Across females, individuals seen fruiting in one year also more likely to be fruiting in next
# Fruit in t0 not independent of fruit in t1 across females
# Seems that some females fruit a lot and others fruit less

# For a given female, however, assume that probability of fruiting is constant
# Probability that individual is female is the product of the probabilities of being female each year given fruit has not been observed
uk_predict_sex_df <- long_uk_sex_df %>%
                     group_by(Indiv_ID) %>%
                     summarise(cumul_P_uk_fem = prod(P_uk_fem, na.rm = TRUE))

# If no probability of being female is predicted, e.g., because logArea is missing, cumul_P_uk_fem will return 1
# Assign their sex as unknown
# If cumulative probability of being female > 0.5, assign female
# If cumulative probability of being female < 0.5, assign female
uk_predict_sex_df <- mutate(uk_predict_sex_df, 
                            Assigned_sex = case_when(cumul_P_uk_fem == 1 ~ "Unknown",
                                                     cumul_P_uk_fem < 1 & cumul_P_uk_fem > 0.5 ~ "Female",
                                                     cumul_P_uk_fem < 0.5 ~ "Male"))

# Combine unknown sex individuals with known females, create sex variable incorporating known and assigned sexes
# Individuals which grow from juvenile to adult will be assigned Unknown sex as a juvenile
long_sex_df <- full_join(long_sex_df, uk_predict_sex_df, by = "Indiv_ID")
long_sex_df <- mutate(long_sex_df, Sex = case_when(Known_sex == "Female" & logArea >= min_size_rep ~ "Female",
                                                   logArea < min_size_rep ~ "Unknown",
                                      Known_sex == "Unknown" & Assigned_sex == "Female" ~ "Female",
                                      Known_sex == "Unknown" & Assigned_sex == "Male" ~ "Male",
                                      Known_sex == "Unknown" & Assigned_sex == "Unknown" ~ "Unknown"))

# Create stage variable to separate juveniles and adults
long_stage_df <- mutate(long_sex_df, 
                        Stage = case_when(logArea < min_size_rep ~ "Juvenile",
                                          is.na(logArea) == TRUE ~ "Unknown",
                                          .default = "Adult"))

# Remove Known_sex, N_fru_sampled, cumul_P_uk_fem, Assigned_sex variables
wrangled_by_yr_df <- select(long_stage_df, -c(Known_sex, N_fru_sampled, cumul_P_uk_fem, Assigned_sex))

# Create measurement ID to describe individual and year for wrangled data frame and original long data frame
wrangled_by_yr_df$Measure_ID <- paste0(wrangled_by_yr_df$Indiv_ID, "_", wrangled_by_yr_df$Year)
long_df$Measure_ID <- paste0(long_df$Indiv_ID, "_", long_df$Year_t0)

# Add sex and stage to data frame describing transitions from t0 to t1
wrangled_df <- full_join(long_df, wrangled_by_yr_df[, c("Sex", "Stage", "Measure_ID")])

# Sex ratio = total males/all adults, counting each adult only once
each_adult_df <- filter(wrangled_df, Stage == "Adult") %>%
                 filter(Sex != "Unknown") %>%
                 distinct(Indiv_ID, .keep_all = TRUE)
sex_ratio <- nrow(filter(each_adult_df, Sex == "Male")) / nrow(each_adult_df)
# Close to sex ratio predicted earlier - 0.46:0.54 M:F

# Test effect of number of times sampled on probability of being male or female

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
ggplot(data = hst_df_raw, aes(x = Host_height, fill = Genus)) +
  geom_histogram()
# Highly colinear - hard to separate effects

# Visualise distribution of mistletoe heights and species
ggplot(data = mst_hst_df, aes(x = Combined_height, fill = Genus)) +
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
ggplot(data = wrangled_by_yr_df, aes(x = logArea)) +
  geom_histogram() +
  geom_vline(xintercept = mean_logArea - 3 * sd_logArea) +
  geom_vline(xintercept = mean_logArea + 3 * sd_logArea)

# Remove outliers +- 3SD
mean_logArea <- mean(wrangled_by_yr_df$logArea, na.rm = TRUE)
sd_logArea <- sd(wrangled_by_yr_df$logArea, na.rm = TRUE)
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

# Visualise relationship between areas and height
ggplot(data = wrangled_by_yr_df, aes(x = logArea, y = Height)) +
  geom_point() +
  theme_bw()
# Size not correlated with height - CHECK IF WAS TRUE BEFORE DISTANCE CORRECTION

# Compare area and height before distance correction


# Construct and select models --------------------------------------------------------

## Model survival ----------------------------------------------------------
# Create binary variable for survival 
wrangled_df <- mutate(wrangled_df, Survive = case_when(
  Status_t1 == "Dead" ~ 0,
  Status_t1 == "Surv" ~ 1
))

# Plot Survive ~ log(Area)
ggplot(data = wrangled_df, aes(x = logArea_t0, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# HOW TO TREAT RANDOM EFFECTS

# Model Survive ~ log(Area) with Indiv_ID as random effect
sur_area_glmm <- glmer(Survive ~ logArea_t0 + (logArea_t0 | Indiv_ID), 
                     data = wrangled_df, family = binomial)
summary(sur_area_glmm)

# Plot Survival ~ Height
ggplot(data = wrangled_df, aes(x = Height, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "Height above ground (m)", y = "Probability of survival to t1") +
  theme_bw()

# Model Survive ~ Height with Indiv_ID as random effect
sur_ht_glmm <- glmer(Survive ~ Height + (Height | Indiv_ID), 
                        data = wrangled_df, family = binomial)
summary(sur_ht_glmm)

# Plot Survival ~ log(Area) + Height
ggplot(data = wrangled_df, aes(x = logArea_t0, y = Survive, col = Height)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survive ~ log(Area) + Height with Indiv_ID as random effect
sur_area_ht_glmm <- glmer(Survive ~ logArea_t0 + Height + (logArea_t0 | Indiv_ID), 
                     data = wrangled_df, family = binomial)
summary(sur_area_ht_glmm)
# Height has slightly stronger relationship with survival than area

# Plot Survival ~ log(Area) for different stages and sexes - juveniles given "Unknown" sex
ggplot(data = filter(wrangled_df, is.na(Sex) == FALSE), aes(x = logArea_t0, y = Survive, col = Sex)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "log(Area) in t0", y = "Probability of survival to t1") +
  theme_bw()

# Filter wrangled data frame for adults only
wrangled_adults_df <- filter(wrangled_df, Sex == "Female" | Sex == "Male")

# Plot Survival ~ log(Area) + Sex for adults
ggplot(data = wrangled_adults_df, aes(x = logArea_t0, y = Survive, col = Sex)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# For adults, model Survival ~ log(Area) + Sex with Indiv_ID as random effect
sur_area_sex_glmm <- glmer(Survive ~ logArea_t0 + Sex + (logArea_t0 | Indiv_ID), 
                     data = wrangled_adults_df, family = binomial)
summary(sur_area_sex_glmm)

# For adults, model Survival ~ log(Area) * Sex with Indiv_ID as random effect
sur_area_sex_int_glmm <- glmer(Survive ~ logArea_t0 * Sex + (logArea_t0 | Indiv_ID), 
                           data = wrangled_adults_df, family = binomial)
summary(sur_area_sex_int_glmm)

# Plot Survival ~ Height + Sex for adults
ggplot(data = wrangled_adults_df, aes(x = Height, y = Survive, col = Sex)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# For adults, model Survival ~ Height + Sex with Indiv_ID as random effect
sur_ht_sex_glmm <- glmer(Survive ~ Height + Sex + (Height | Indiv_ID), 
                           data = wrangled_adults_df, family = binomial)
summary(sur_ht_sex_glmm)

# For adults, model Survival ~ Height * Sex with Indiv_ID as random effect
sur_ht_sex_int_glmm <- glmer(Survive ~ Height * Sex + (Height | Indiv_ID), 
                               data = wrangled_adults_df, family = binomial)
summary(sur_ht_sex_int_glmm)

# For adults, model Survival ~ logArea + Height + Sex with Indiv_ID as random effect
sur_area_ht_sex_glmm <- glmer(Survive ~ logArea_t0 + Height + Sex + (logArea_t0 | Indiv_ID), 
                             data = wrangled_adults_df, family = binomial)
summary(sur_area_ht_sex_glmm)

# For adults, model Survival ~ logArea * Sex + Height with Indiv_ID as random effect
sur_area_ht_sex_int_glmm <- glmer(Survive ~ logArea_t0 * Sex + Height + (logArea_t0 | Indiv_ID), 
                              data = wrangled_adults_df, family = binomial)
summary(sur_area_ht_sex_int_glmm)

# Survival depends on sex, therefore model juveniles, females and males separately
# Filter wrangled data for juveniles only
wrangled_juv_df <- filter(wrangled_df, 
                          Stage == "Juvenile")

# Filter wrangled data for females only
wrangled_fem_df <- filter(wrangled_df, 
                          Sex == "Female")

# Filter wrangled data for males only
wrangled_mal_df <- filter(wrangled_df, 
                          Sex == "Male")


### Juvenile survival -------------------------------------------------------

# Plot Survival ~ log(Area) for juveniles
ggplot(data = wrangled_juv_df, aes(x = logArea_t0, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ log(Area) for juveniles with Indiv_ID as random effect
sur_juv_area_glmm <- glmer(Survive ~ logArea_t0 + (1 | Indiv_ID), 
                           data = wrangled_juv_df, family = binomial)
summary(sur_juv_area_glmm)

# Plot Survival ~ Height for juveniles
ggplot(data = wrangled_juv_df, aes(x = Height, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ Height for juveniles with Indiv_ID as random effect
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

# Neither log(Area) nor height significantly predict juvenile survival

### Female survival ---------------------------------------------------------
# Plot Survival ~ log(Area) for females
ggplot(data = wrangled_fem_df, aes(x = logArea_t0, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ log(Area) for females with Indiv_ID as random effect
sur_fem_area_glmm <- glmer(Survive ~ logArea_t0 + (logArea_t0 | Indiv_ID), 
                           data = wrangled_fem_df, family = binomial)
summary(sur_fem_area_glmm)

# Plot Survival ~ Height for females
ggplot(data = wrangled_fem_df, aes(x = Height, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ Height for females with Indiv_ID as random effect
sur_fem_ht_glmm <- glmer(Survive ~ Height + (Height | Indiv_ID), 
                         data = wrangled_fem_df, family = binomial)
summary(sur_fem_ht_glmm)

# Plot Survival ~ log(Area) + Height for females
ggplot(data = wrangled_fem_df, aes(x = logArea_t0, y = Survive, col = Height)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ log(Area) + Height for females with Indiv_ID as random effect
sur_fem_area_ht_glmm <- glmer(Survive ~ logArea_t0 + Height + (1 | Indiv_ID), 
                              data = wrangled_fem_df, family = binomial)
summary(sur_fem_area_ht_glmm)

# Survival is predicted by size for females

### Male survival -----------------------------------------------------------
# Plot Survival ~ log(Area) for males
ggplot(data = wrangled_mal_df, aes(x = logArea_t0, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ log(Area) for males with Indiv_ID as random effect
sur_mal_area_glmm <- glmer(Survive ~ logArea_t0 + (logArea_t0 | Indiv_ID), 
                           data = wrangled_mal_df, family = binomial)
summary(sur_mal_area_glmm)

# Plot Survival ~ Height for males
ggplot(data = wrangled_mal_df, aes(x = Height, y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ Height for males with Indiv_ID as random effect
sur_mal_ht_glmm <- glmer(Survive ~ Height + (Height | Indiv_ID), 
                         data = wrangled_mal_df, family = binomial)
summary(sur_mal_ht_glmm)

# Plot Survival ~ log(Area) + Height for males
ggplot(data = wrangled_mal_df, aes(x = logArea_t0, y = Survive, col = Height)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  theme_bw()

# Model Survival ~ log(Area) + Height for males with Indiv_ID as random effect
sur_mal_area_ht_glmm <- glmer(Survive ~ logArea_t0 + Height + (1 | Indiv_ID), 
                              data = wrangled_mal_df, family = binomial)
summary(sur_mal_area_ht_glmm)

# Survival is predicted by log(Area) in males

## Model growth ------------------------------------------------------------
# Plot growth in time t1 against growth in time t0
ggplot(data = wrangled_df, aes(x = logArea_t0, y = logArea_t1)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  lims(x = c(0, max(long_df$logArea_t0)), y = c(0, max(long_df$logArea_t1))) +
  geom_smooth(method = "lm") +
  theme_bw()

# Plot growth in time t1 against growth in time t0 as a function of sex
ggplot(data = wrangled_df, aes(x = logArea_t0, y = logArea_t1, col = Sex)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  lims(x = c(0, max(long_df$logArea_t0)), y = c(0, max(long_df$logArea_t1))) +
  geom_smooth() +
  theme_bw() +
  geom_vline(xintercept = min_size_rep)

# Calculate per-time-step growth rate for each individual in each time step
wrangled_df$Growth_t0_t1 <- (wrangled_df$logArea_t1 - wrangled_df$logArea_t0)/wrangled_df$logArea_t0

# Plot relative growth from t0 to t1 against growth in time t0
ggplot(data = wrangled_df, aes(x = logArea_t0, y = Growth_t0_t1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() 

# Plot relative growth from t0 to t1 against growth in time t0 as a function of sex
ggplot(data = wrangled_df, aes(x = logArea_t0, y = Growth_t0_t1, col = Sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  geom_vline(xintercept = min_size_rep) +
  geom_hline(yintercept = 0)

# Growth depends on sex, therefore model juveniles, males and females separately

# HOW TO MODEL GROWTH - MODEL RELATIVE GROWTH AND RE-CONVERT?

# MODEL SURVIVAL SEPARATELY FOR MATURING AND NON-MATURING

# Split juveniles into those which matured (i.e., of reproductive size in t1) and those which did not
wrangled_juv_df <- mutate(wrangled_juv_df, 
                          Maturation = case_when()




### Juvenile growth ---------------------------------------------------------


### Female growth -----------------------------------------------------------


### Male growth -------------------------------------------------------------




## Model fruiting ------------------------------------------------------
# Filter wrangled data for females only
wrangled_fem_by_yr_df <- filter(wrangled_by_yr_df, 
                          Sex == "Female")

# Plot fruiting ~ logArea for females only 
ggplot(data = wrangled_fem_by_yr_df, aes(x = logArea, y = Fruit)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "log(Area) in t0", y = "Probability of fruiting in t0") +
  theme_bw()

# Model fruiting ~ logArea with Indiv_ID as random effect
fru_area_glmm <- glmer(Fruit ~ logArea + (logArea | Indiv_ID), 
                       data = wrangled_fem_by_yr_df, family = binomial)
summary(fru_area_glmm)

# Plot fruiting ~ Height for females only 
ggplot(data = wrangled_fem_by_yr_df, aes(x = Height, y = Fruit)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "Height above ground (m)", y = "Probability of fruiting in t0") +
  theme_bw()

# Model fruiting ~ Height with Indiv_ID as random effect
fru_ht_glmm <- glmer(Fruit ~ Height + (Height | Indiv_ID), 
                     data = wrangled_fem_by_yr_df, family = binomial)
summary(fru_ht_glmm)

# Plot fruiting ~ log(Area) + Height for females only 
ggplot(data = wrangled_fem_by_yr_df, aes(x = logArea, y = Fruit, col = Height)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial)) +
  labs(x = "log(Area) in t0", y = "Probability of fruiting in t0") +
  theme_bw()

# Model fruiting ~ logArea + Height with Indiv_ID as random effect
fru_area_ht_glmm <- glmer(Fruit ~ logArea + Height + (logArea | Indiv_ID), 
                     data = wrangled_fem_by_yr_df, family = binomial)
summary(fru_area_ht_glmm)

# Fruiting in females is predicted by both height and sex



# Construct IPM -----------------------------------------------------------

## How to incorporate height if selected?!

# Function to calculate meshpoints
calculate_mesh_points <- function(number_of_cells, min_size, max_size) {
  # number of cells in the discretized kernel
  n = number_of_cells
  # Boundary points (the edges of the cells defining the kernel)
  b = min_size + c(0:n) * (max_size-min_size) / n 
  # Mesh points (midpoints of the cells)
  y = 0.5 * (b[1:n] + b[2:(n+1)])
  return(y)
}

# Juvenile meshpoints
min_juv_area <- log(min(host_paras_rep_comb$Area, na.rm = TRUE))
max_juv_area <- log(min_fem_area)
J_meshpoints <- calculate_mesh_points(number_of_cells = 100,
                                      min_size = min_juv_area,
                                      max_size = max_juv_area)

# Adult meshpoints
min_adu_area <- log(min_fem_area)
max_adu_area <- log(max(host_paras_rep_comb$Area, na.rm = TRUE))
A_meshpoints <- calculate_mesh_points(number_of_cells = 100,
                                      min_size = min_adu_area,
                                      max_size = max_adu_area)


# Function to calculate remain P kernels: JJ
build_p_kernel_juv_remain <- function(y, growth_function, survival_function, mature_function, params) {
  # Width of the cells
  h = y[2] - y[1] 
  G = h * outer(y, y, growth_function, params) # growth sub-kernel (which does not yet have survival in it)
  S = survival_function(y, params) # survival 
  M = mature_function(y, params) # probability of not remaining in the same stage (e.g, juv to juv or adu to adu)
  kernel = G # placeholder; we're about to redefine P on the next line
  n = length(y)
  try(for(i in 1:n) kernel[,i] <- G[,i] * S[i] * (1-M[i]))  # growth/survival/remain sub-kernel
  return(kernel)
}

# Function to calculate remain P kernels: FF, MM
build_p_kernel_adult <- function(y, growth_function, survival_function, params) {
  # Width of the cells
  h = y[2] - y[1] 
  G = h * outer(y, y, growth_function, params) # growth sub-kernel (which does not yet have survival in it)
  S = survival_function(y, params) # survival 
  kernel = G # placeholder; we're about to redefine P on the next line
  n = length(y)
  try(for(i in 1:n) kernel[,i] <- G[,i] * S[i])  # growth/survival sub-kernel
  return(kernel)
}

# Function to calculate leave P kernels: FJ, MJ
build_p_kernel_mature <- function(y0, y1, growth_function, survival_function, mature_function, sex_ratio, params) {
  # Width of the cells
  h = y0[2] - y0[1] 
  G = h * outer(y1, y0, growth_function, params) # growth sub-kernel (which does not yet have survival in it)
  S = survival_function(y0, params) # survival 
  M = mature_function(y0, params) # probability of moving to a different stage (e.g, juv to adu or adu to juv)
  n = length(y0)
  kernel = G # placeholder; we're about to redefine P on the next line
  try(for(i in 1:n) kernel[,i] <- G[,i] * S[i] * M[i] * sex_ratio)  # growth/survival/leave sub-kernel
  return(kernel)
}

# Function to calculate F kernels: JF, JM
build_f_kernel <- function(y0, y1, reprod_function, growth_function, mean_fec, params) {
  # Width of the cells
  h = y0[2] - y0[1] 
  F = h * outer(y0, y1, reprod_function, params)
  n = length(y1)
  kernel = h * outer(y1, y0, growth_function, params) # placeholder; we're about to redefine F on the next line
  try(for(i in 1:n) kernel[,i] <- 0.5 * F[,i] * mean_fec)  # growth/survival/leave sub-kernel
  return(kernel)
}

# Function to build empty kernel: MF, FM
build_empty_kernel <- function(y0, y1, growth_function, params) {
  h = y0[2] - y0[1] 
  n = length(y1)
  kernel = h * outer(y1, y0, growth_function, params) # placeholder; we're about to redefine F on the next line
  try(for(i in 1:n) kernel[,i] <- 0)  # growth/survival/leave sub-kernel
  return(kernel)
}


#JJ: P Kernel for transitions from juvenile in time t to juvenile in time t+1
#   JJ = f(juvenile survival x prob remaining juvenile x juvenile growth to juvenile)

#FJ: P Kernel for transitions from juvenile in time t to female in time t+1
#   FJ = f(juvenile survival x (1 - prob remaining juvenile) x sex ratio x juvenile growth to adult)

#MJ: P Kernel for transitions from juvenile in time t to male in time t+1
#   MJ = f(juvenile survival x (1 - prob remaining juvenile) x (1 - sex ratio) x juvenile growth to adult)

#FF: P Kernel for transitions from female in time t to female in time t+1
#   FF = f(female survival x adult growth to adult)

#MM: P Kernel for transitions from male in time t to male in time t+1
#   MM = f(male survival x adult growth to adult)

#JF: F Kernel for transitions from female in time t to juvenile in time t+1
#   JF = f(0.5 x female fruiting x mean number of offspring x offspring size distribution x offspring relative height distribution)

#JM: F Kernel for transitions from female in time t to juvenile in time t+1
#   JM = f(0.5 x female fruiting x mean number of offspring x offspring size distribution x offspring relative height distribution)


## JJ: J -> J
## JJ = f(juvenile survival x prob remaining juvenile x juvenile growth to juvenile)

# Define functions for JJ p-kernel
# 1. Probability of surviving
survival_jj_function <- function(x,params) {
  u = exp(params$juv.surv.int + params$juv.surv.slope * x)
  return(u/(1+u))
}

# 2. Probability of maturing
mature_jj_function <- function(x,params) {
  u = exp(params$juv.mat.int + params$juv.mat.slope * x)
  return(u/(1+u))
}

# 3. Growth function
growth_jj_function <- function(xp,x,params) { 			
  dnorm(xp, mean = params$juv.growth.int + params$juv.growth.slope * x, 
        sd = params$juv.growth.sd)
}

# Build the JJ P kernel
JJ_kernel <- build_p_kernel_juv_remain(y = J_meshpoints,
                                       growth_function = growth_jj_function,
                                       survival_function = survival_jj_function,
                                       mature_function = mature_jj_function,
                                       params = params)

# Plot the kernel
cols <- brewer.pal(9, "Blues")
cols
pal <- colorRampPalette(cols)
image(JJ_meshpoints,JJ_meshpoints,t(log(JJ_kernel+1)), xlab = "Log area (t)",ylab = "Log area (t+1)", col = pal(100))
abline(0,1, col = "black", lty=2)


## FJ: J -> F
## FJ = f(juvenile survival x (1 - prob remaining juvenile) x sex ratio x juvenile growth to adult)

# Define functions for FJ p-kernel
# 1. Probability of surviving same as for all juveniles
survival_fj_function <- function(x,params) {
  u = exp(params$juv.surv.int + params$juv.surv.slope * x)
  return(u/(1+u))
}

# 2. Probability of maturing
mature_fj_function <- function(x,params) {
  u = exp(params$juv.mat.int + params$juv.mat.slope * x)
  return(u/(1+u))
}

# 3. Growth function same as for all juveniles
growth_fj_function <- function(xp,x,params) { 			
  dnorm(xp, mean = params$juv.growth.int + params$juv.growth.slope * x, 
        sd = params$juv.growth.sd)
}

# Build the FJ P kernel
FJ_kernel <- build_p_kernel_mature(y0 = J_meshpoints,
                                   y1 = A_meshpoints,
                                   growth_function = growth_fj_function,
                                   survival_function = survival_fj_function,
                                   mature_function = mature_fj_function,
                                   sex_ratio = params$sex.ratio,
                                   params = params)

# Plot the kernel
cols <- brewer.pal(9, "Blues")
cols
pal <- colorRampPalette(cols)
image(J_meshpoints,A_meshpoints,t(log(FJ_kernel+1)), xlab = "Log area (t)",ylab = "Log area (t+1)", col = pal(100))
abline(0,1, col = "black", lty=2)

## MJ: J -> M
## MJ = f(juvenile survival x (1 - prob remaining juvenile) x (1 - sex ratio) x juvenile growth to adult)

# Define functions for MJ p-kernel
# 1. Probability of surviving same as for all juveniles
survival_mj_function <- function(x,params) {
  u = exp(params$juv.surv.int + params$juv.surv.slope * x)
  return(u/(1+u))
}

# 2. Probability of maturing
mature_mj_function <- function(x,params) {
  u = exp(params$juv.mat.int + params$juv.mat.slope * x)
  return(u/(1+u))
}

# 3. Growth function same as for all juveniles
growth_mj_function <- function(xp,x,params) { 			
  dnorm(xp, mean = params$juv.growth.int + params$juv.growth.slope * x, 
        sd = params$juv.growth.sd)
}

# Build the MJ P kernel
MJ_kernel <- build_p_kernel_mature(y0 = J_meshpoints,
                                   y1 = A_meshpoints,
                                   growth_function = growth_fj_function,
                                   survival_function = survival_fj_function,
                                   mature_function = mature_fj_function,
                                   sex_ratio = 1 - params$sex.ratio,
                                   params = params)

# Plot the kernel
cols <- brewer.pal(9, "Blues")
cols
pal <- colorRampPalette(cols)
image(J_meshpoints,A_meshpoints,t(log(MJ_kernel+1)), xlab = "Log area (t)",ylab = "Log area (t+1)", col = pal(100))
abline(0,1, col = "black", lty=2)


## FF: F -> F
## FF = f(female survival x adult growth to adult)
# Define functions for FF p-kernel
# 1. Probability of surviving
survival_ff_function <- function(x,params) {
  u = exp(params$fem.surv.int + params$adu.surv.slope.area * x) # ignore rht for now
  return(u/(1+u))
}

# 2. Growth function
growth_ff_function <- function(xp,x,params) { 			
  dnorm(xp, mean = params$adu.growth.int + params$adu.growth.slope * x, 
        sd = params$adu.growth.sd)
}

# Build the FF P kernel
FF_kernel <- build_p_kernel_adult(y = A_meshpoints,
                                  growth_function = growth_ff_function,
                                  survival_function = survival_ff_function,
                                  params = params)

# Plot the kernel
cols <- brewer.pal(9, "Blues")
cols
pal <- colorRampPalette(cols)
image(A_meshpoints,A_meshpoints,t(log(FF_kernel+1)), xlab = "Log area (t)",ylab = "Log area (t+1)", col = pal(100))
abline(0,1, col = "black", lty=2)

## MM: M -> M
## MM = f(male survival x adult growth to adult)
# Define functions for MM p-kernel
# 1. Probability of surviving
survival_mm_function <- function(x,params) {
  u = exp(params$male.surv.int + params$adu.surv.slope.area * x) # ignore rht for now
  return(u/(1+u))
}

# 2. Growth function
growth_mm_function <- function(xp,x,params) { 			
  dnorm(xp, mean = params$adu.growth.int + params$adu.growth.slope * x, 
        sd = params$adu.growth.sd)
}

# Build the MM P kernel
MM_kernel <- build_p_kernel_adult(y = A_meshpoints,
                                  growth_function = growth_mm_function,
                                  survival_function = survival_mm_function,
                                  params = params)

# Plot the kernel
cols <- brewer.pal(9, "Blues")
cols
pal <- colorRampPalette(cols)
image(A_meshpoints,A_meshpoints,t(log(MM_kernel+1)), xlab = "Log area (t)",ylab = "Log area (t+1)", col = pal(100))
abline(0,1, col = "black", lty=2)

## JF: F -> J
##  JF = f(0.5 x female fruiting x mean number of offspring x offspring size distribution x offspring relative height distribution)
# Define functions for JF F kernel
# 1. Reproduction function - Probability of Fruiting multiplied recruitment distribution
reprod_jf_function <- function(xp,x,params) {
  u = exp(params$fruit.int + params$fruit.slope * x)
  (u/(1+u))*dnorm(xp, mean = params$recruit.size.mean, sd = params$recruit.size.sd)
}

# 2. Null growth distribution for growth - used as placeholder for kernel
growth_jf_function <- function(xp,x,params) { 			
  dnorm(xp, mean = params$adu.growth.int + params$adu.growth.slope * x, 
        sd = params$adu.growth.sd)
}


# Build the JF F kernel
JF_kernel <- build_f_kernel(y0 = J_meshpoints,
                            y1 = A_meshpoints,
                            reprod_function = reprod_jf_function,
                            growth_function = growth_jf_function,
                            mean_fec = params$fec.mean,
                            params = params)

# Plot the kernel
cols <- brewer.pal(9, "Blues")
cols
pal <- colorRampPalette(cols)
image(A_meshpoints,J_meshpoints,t(log(JF_kernel+1)), xlab = "Log area (t)",ylab = "Log area (t+1)", col = pal(100))
abline(0,1, col = "black", lty=2)

## JM: M -> J
## JM = f(0.5 x female fruiting x mean number of offspring x offspring size distribution x offspring relative height distribution)
# Define functions for JM F kernel
# 1. Reproduction function - Probability of Fruiting multiplied recruitment distribution
reprod_jm_function <- function(xp,x,params) {
  u = exp(params$fruit.int + params$fruit.slope * x)
  (u/(1+u))*dnorm(xp, mean = params$recruit.size.mean, sd = params$recruit.size.sd)
}

# 2. Null growth distribution for growth - used as placeholder for kernel
growth_jm_function <- function(xp,x,params) { 			
  dnorm(xp, mean = params$adu.growth.int + params$adu.growth.slope * x, 
        sd = params$adu.growth.sd)
}

# Build the JF F kernel
JM_kernel <- build_f_kernel(y0 = J_meshpoints,
                            y1 = A_meshpoints,
                            reprod_function = reprod_jm_function,
                            growth_function = growth_jm_function,
                            mean_fec = params$fec.mean,
                            params = params)

# Plot the kernel
cols <- brewer.pal(9, "Blues")
pal <- colorRampPalette(cols)
image(A_meshpoints,J_meshpoints,t(log(JF_kernel+1)), xlab = "Log area (t)",ylab = "Log area (t+1)", col = pal(100))

## MF: empty kernel
# 1. Null growth distribution for growth - used as placeholder for kernel
growth_mf_function <- function(xp,x,params) { 			
  dnorm(xp, mean = params$adu.growth.int + params$adu.growth.slope * x, 
        sd = params$adu.growth.sd)
}

# Build empty MF kernel
MF_kernel <- build_empty_kernel(y0 = A_meshpoints,
                                y1 = A_meshpoints,
                                growth_function = growth_mf_function,
                                params = params)
# Plot the kernel
cols <- brewer.pal(9, "Blues")
pal <- colorRampPalette(cols)
image(A_meshpoints,A_meshpoints,t(log(MF_kernel+1)), xlab = "Log area (t)",ylab = "Log area (t+1)", col = pal(100))

## FM: empty kernel
# 1. Null growth distribution for growth - used as placeholder for kernel
growth_fm_function <- function(xp,x,params) { 			
  dnorm(xp, mean = params$adu.growth.int + params$adu.growth.slope * x, 
        sd = params$adu.growth.sd)
}

# Build empty MF kernel
FM_kernel <- build_empty_kernel(y0 = A_meshpoints,
                                y1 = A_meshpoints,
                                growth_function = growth_fm_function,
                                params = params)
# Plot the kernel
cols <- brewer.pal(9, "Blues")
pal <- colorRampPalette(cols)
image(A_meshpoints,A_meshpoints,t(log(FM_kernel+1)), xlab = "Log area (t)",ylab = "Log area (t+1)", col = pal(100))


#### Putting it all together using cbind and rbind ####
MJ_MF_MM <- cbind(MJ_kernel, MF_kernel, MM_kernel)
FJ_FF_FM <- cbind(FJ_kernel, FF_kernel, FM_kernel)
JJ_JF_JM <- cbind(JJ_kernel, JF_kernel, JM_kernel)
IPM <- rbind(JJ_JF_JM, FJ_FF_FM, MJ_MF_MM)


# Analyse IPM -------------------------------------------------------------
# Calculate asymptotic growth rate


# Calculate generation time




