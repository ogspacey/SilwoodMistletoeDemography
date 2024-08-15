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
# Combine mistletoe demographic data and host data
mst_hst_df <- full_join(mst_df_raw, hst_df_raw, by = "Host_ID")

# Create individual ID combining mistletoe and host IDs
mst_hst_df <- mutate(mst_hst_df, Indiv_ID = paste(Mistletoe_ID, Host_ID, sep = "_"))

# Remove columns measuring perimeter as we will only use area as an indicator of size
mst_hst_df <- select(mst_hst_df, -contains('Perim'))

# Remove mistletoes that are never seen
mst_hst_df <- filter(mst_hst_df, Notes != "Never seen")

# Estimate distance to mistletoes (hypotenuse) for mistletoes with imputed height via Pythagoras, to rectify size perception
mst_hst_df <- mutate(mst_hst_df, 
                     Estimated_hypotenuse = sqrt((mst_hst_df$Imputed_height)^2 + (mst_hst_df$Host_horizontal_distance)^2))

# Convert table to long format to show status, area and fruit in time t0 and status, area and fruit in time t1
years <- c(14:22) # vector of years to run for loop through
long_df_ls <- list()   # list of long data frames to store each year
for(year in years){ # run through years sequentially
  Status_t0   <- paste("Status", year,   sep = "") 
  Area_t0     <- paste("Area",   year,   sep = "")
  Fruit_t0    <- paste("Fruit",  year,   sep = "")
  Status_t1   <- paste("Status", year+1, sep = "")
  Area_t1     <- paste("Area",   year+1, sep = "")
  Fruit_t1    <- paste("Fruit",  year+1, sep = "")
  mst_df_year <- data.frame(Indiv_ID = mst_hst_df$Indiv_ID,  # Individual IDs
                            Year_t0 = year,
                            Status_t0 = mst_hst_df[[Status_t0]],
                            Area_t0 = mst_hst_df[[Area_t0]],
                            Fruit_t0 = mst_hst_df[[Fruit_t0]],
                            Year_t1 = year + 1,
                            Status_t1 = mst_hst_df[[Status_t1]],
                            Area_t1 = mst_hst_df[[Area_t1]],
                            Fruit_t1 = mst_hst_df[[Fruit_t1]]
                            )
  long_df_ls[[year]] <- mst_df_year
}

# Merge rows to create long data frame
long_df <- bind_rows(long_df_ls)

# Rectify size based on horizontal distance to base of tree (standard)
# If an object is twice as far away relative to the standard, it will appear twice as small in length, and four times as small in area
# Calculate the area adjusting for distance relative to the standard (horizontal distance) as Adjusted area = (Distance to object/Distance to standard)^2 * Measured area
Adjusted_Area <- (Hypotenuse/Horizontal_distance)^2 * Area

# Explore data --------------------------------------------------------
# Visualise distribution of mistletoe sizes

# Convert 2D area to log(area)

# Visualise distribution of mistletoe heights

# Visualise distribution of host

# Effect of density dependence


# Assign stage and sex ----------------------------------------------------

# Set cut-off for juveniles and adults

# Assign known females 


# Construct and select linear models --------------------------------------------------------

## Model survival ----------------------------------------------------------

long_df <- mutate(long_df, Survive = case_when(
  Status_t1 == "Dead" ~ 0,
  Status_t1 == "Surv" ~ 1
))

ggplot(data = long_df, aes(x = log(Area_t0), y = Survive)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial))

# Select model for growth

log_reg_sur_1 <- glm(Survive ~ log(Area_t0), data = long_df, family = binomial)
summary(log_reg_sur_1)

## Model growth ------------------------------------------------------------
# Plot growth in time t1 against growth in time t0
ggplot(data = long_df, aes(x = log(Area_t0), y = log(Area_t1), col = Year_t0)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1)) +
  lims(x = c(0, log(max(long_df$Area_t0))), y = c(0, log(max(long_df$Area_t1)))) +
  geom_smooth() +
  theme_bw()

# Select model for growth

## Model fruiting ------------------------------------------------------

# Select model for fruiting
ggplot(data = long_df, aes(x = log(Area_t0), y = Fruit_t0)) +
  geom_point() +
  stat_smooth(method = "glm", 
              method.args = list(family = binomial))

log_reg_fru_1 <- glm(Fruit_t0 ~ log(Area_t0), data = long_df, family = binomial)
summary(log_reg_fru_1)

# Construct IPM -----------------------------------------------------------


# Analyse IPM -------------------------------------------------------------
# Calculate asymptotic growth rate


# Calculate generation time




