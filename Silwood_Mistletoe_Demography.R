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
mst_hst_df <- full_join(mst_df_raw, hst_df_raw, by = Host_ID)

# Create individual ID combining mistletoe and host IDs
indiv_id_df <- mutate(mst_hst_df, Indiv_ID = paste(Mistletoe_ID, Host_ID, sep = "_"))

# Remove columns measuring perimeter as we will only use area as an indicator of size
no_perim_df <- select(indiv_df, -contains('Perim'))

# Estimate distance to mistletoes (hypotenuse) for mistletoes with imputed height via Pythagoras, to rectify size perception
est_hyp_df <- mutate(no_perim_df, 
                     Estimated_hypotenuse = sqrt((no_perim_df$Imputed_height)^2 + (no_perim_df$Host_horizontal_distance)^2))

# Convert table to long format to show status, area and fruit in time t0 and status, area and fruit in time t1
years <- c(14:22) # vector of years to run for loop through
df_ls <- list()   # list of data frames to store each year
for(year in years){ # run through years sequentially
  Status_t0   <- paste("Status", year, sep = "") 
  Area_t0     <- paste("Area", year, sep = "")
  Fruit_t0    <- paste("Fruit", year, sep = "")
  Status_t1 <- paste("Status", year+1, sep = "")
  Area_t1   <- paste("Area", year+1, sep = "")
  Fruit_t1  <- paste("Fruit", year+1, sep = "")
  mst_df_year <- data.frame(Indiv_ID = previous_df$Indiv_ID,
                            Year_t0 = year,
                            Status_t0 = previous_df[[Status_t0]],
                            Area_t0 = previous_df[[Area_t0]],
                            Fruit_t0 = previous_df[[Fruit_t0]],
                            Year_t1 = year + 1,
                            Status_t1 = previous_df[[Status_t1]],
                            Area_t1 = previous_df[[Area_t1]],
                            Fruit_t1 = previous_df[[Fruit_t1]]
                            )
  df_ls[as.character(year)]<- mst_df_year
}


stat_long_df <- pivot_longer(no_perim_df, cols = Status15, names_to = "Year", values_to = "Status_t")

# Rectify size based on horizontal distance to base of tree (standard)
# If an object is twice as far away relative to the standard, it will appear twice as small in length, and four times as small in area
# Calculate the area adjusting for distance relative to the standard (horizontal distance) as Adjusted area = (Distance to object/Distance to standard)^2 * Measured area
Adjusted_Area <- (Hypotenuse/Horizontal_distance)^2 * Area

# Explore data --------------------------------------------------------
# Visualise distribution of mistletoe sizes

# Convert 2D area to log(area)

# Visualise distribution of mistletoe heights

# Visualise distribution of host


# Assign stage and sex ----------------------------------------------------

# Set cut-off for juveniles and adults

# Assign known females 


# Construct and select linear models --------------------------------------------------------

## Model survival ----------------------------------------------------------


## Model growth ------------------------------------------------------------



## Model reproduction ------------------------------------------------------


# Construct IPM -----------------------------------------------------------


# Analyse IPM -------------------------------------------------------------
# Calculate asymptotic growth rate


# Calculate generation time




