### Silwood Park Mistletoe - PCA
### Oliver G. Spacey, Owen R. Jones, Sydne Record, Arya Y. Yue, Wenyi Liu, Alice Rosen, Michael Crawley, Chris J. Thorogood, Roberto Salguero-Gómez
### Adapted from Salguero-Gómez, 2024

# Pre-amble -----------------------------------------------------------
# Clear environment
rm(list=ls())

# Load necessary packages
library(scales)
library(mice)
library(popdemo)   # Extract life history traits
library(popbio)    # Extract life history traits
library(remote)
#remotes::install_github("jonesor/Rcompadre", build_opts = NULL)
#remotes::install_github("jonesor/Rage", build_opts = NULL)
library(Rcompadre) # Extract and manipulate data from COMPADRE
library(Rage)      # Extract life history traits
library(naniar)
library(rotl)      # Query Open Tree of Life
library(ape)
library(geiger)
library(nlme)
library(phytools)
library(caper)
library(agricolae)
library(vioplot)
library(stringr)
library(Ternary)
library(plotrix)
library(gclus)
library(ggplot2) # Data visualisation
library(dplyr)   # Data wrangling
library(rworldmap)
library(scales) # Data visualisation
library(RColorBrewer) # Data visualisation
library(khroma)
library(rphylopic)
library(png)       # for reading PNG images
library(grid)      # Data visualisation

# Set working directory as required

# Set colours for visualisation
colours <- color("vibrant")

# Obtain life history traits -----------------------------------------------------

# Call in COMPADRE
compadre <- cdb_fetch("compadre")

# Check dimensions of initial dataset
dim(compadre)

# Subset MPMs for unmanipulated, wild populations, dim >= 3, 
# with survival and fecundity appropriately measured
compadre <- subset(compadre, MatrixTreatment == "Unmanipulated" &
                     MatrixDimension >= 3 &
                     SurvivalIssue < 1.001 &
                     MatrixSplit == "Divided" &
                     MatrixFec =="Yes" &
                     MatrixCaptivity == "W")

# Check dimensions of new dataset
dim(compadre)

# Obtain unique matrix IDs
uniqueStudiesFirstID <- compadre$MatrixID[which(duplicated(compadre$SpeciesAuthor)==F)]

# Subset accordingly to remove duplicates 
compadre <- subset(compadre, compadre$MatrixID %in% uniqueStudiesFirstID)
# Check dimensions of new dataset
dim(compadre)

# Check how many species have more than a single study
tableSpeciesAccepted <- table(compadre$SpeciesAccepted)
matrixSpeciesAccepted <- as.matrix(tableSpeciesAccepted)
speciesDuplicate <- rownames(matrixSpeciesAccepted)[which(matrixSpeciesAccepted>1)]
speciesDuplicate
#[1] "Actaea elata"                   
#[2] "Alliaria petiolata"             
#[3] "Anthyllis vulneraria"           
#[4] "Araucaria hunsteinii"           
#[5] "Ascophyllum nodosum"            
#[6] "Astragalus scaphoides"          
#[7] "Astrophytum capricorne"         
#[8] "Astrophytum ornatum"            
#[9] "Calochortus lyallii"            
#[10] "Carduus nutans"                 
#[11] "Carlina vulgaris"               
#[12] "Chamaedorea radicalis"          
#[13] "Cirsium perplexans"             
#[14] "Cirsium pitcheri"               
#[15] "Cirsium tracyi"                 
#[16] "Cypripedium calceolus"          
#[17] "Dactylorhiza lapponica"         
#[18] "Digitalis purpurea"             
#[19] "Dracocephalum austriacum"       
#[20] "Erythronium japonicum"          
#[21] "Euterpe precatoria"             
#[22] "Hylocomium splendens"           
#[23] "Lepanthes rupestris"            
#[24] "Lomatium bradshawii"            
#[25] "Mammillaria huitzilopochtli"    
#[26] "Mammillaria pectinifera"        
#[27] "Neobuxbaumia macrocephala"      
#[28] "Neobuxbaumia mezcalaensis"      
#[29] "Neobuxbaumia tetetzo"           
#[30] "Panax quinquefolius"            
#[31] "Phyllanthus emblica"            
#[32] "Pinus lambertiana"              
#[33] "Plantago coronopus"             
#[34] "Primula veris"                  
#[35] "Primula vulgaris"               
#[36] "Pseudophoenix sargentii"        
#[37] "Pyrrocoma radiata"              
#[38] "Rhododendron maximum"           
#[39] "Rhododendron ponticum"          
#[40] "Salix arctica"                  
#[41] "Sarracenia purpurea"            
#[42] "Shorea leprosula"               
#[43] "Solidago altissima"             
#[44] "Succisa pratensis"              
#[45] "Taraxacum campylodes"           
#[46] "Tillandsia multicaulis"         
#[47] "Tillandsia punctulata"          
#[48] "Trillium grandiflorum"          
#[49] "Trollius laxus"                 
#[50] "Vella pseudocytisus subsp. paui"

# Run through duplicated species
for (spp in speciesDuplicate){
  try(subset <- compadre[which(compadre$SpeciesAccepted==spp), c("MatrixID","SpeciesAuthor","SpeciesAccepted","CommonName","StudyDuration","MatrixComposite","MatrixDimension")])
  print(subset)
  matrixID <- compadre$MatrixID[which(compadre$SpeciesAccepted==spp)]
  readline(prompt="Press [enter] to continue")
  for (MPM in matrixID){
    print(compadre$mat[[which(compadre$MatrixID == MPM)]])
    readline(prompt="Press [enter] to continue")
  }
}

# Delete studies that are replicated by given species, prioritising those that have longer duration, higher dimensionality, range of sizes examined, consistency in state variables, inclusion of reproduction
delete_Actaea_elata <- c(241325) # Retain 242308 (longer duration)            
delete_Alliaria_petiolata <- c(241453, 241454) # Retain 241455 (longer duration)          
delete_Anthyllis_vulneraria <- c(241556) # Retain 239132 (longer duration)         
delete_Araucaria_hunsteinii <- c(247391) # Retain 247395 (larger matrix)         
delete_Ascophyllum_nodosum <- c(238271) # Retain 240684 (longer duration)         
delete_Astragalus_scaphoides <- c(241822) # Retain 241833 (longer duration)        
delete_Astrophytum_capricorne <- c(246894) # Retain 246895 (longer duration)       
delete_Astrophytum_ornatum <- c(246904) # Retain 246901 (longer duration)         
delete_Calochortus_lyallii <- c(242128) # Retain 242069 (larger matrix)          
delete_Carduus_nutans <- c(242147) # Retain 242153 (longer duration)               
delete_Carlina_vulgaris <- c(242178) # Retain 242187 (longer duration)             
delete_Chamaedorea_radicalis <- c(246265) # Retain 246232 (longer duration)       
delete_Cirsium_perplexans <- c(242360) # Retain 242359 (larger matrix)         
delete_Cirsium_pitcheri <- c(242367, 242373, 242413) # Retain 242422 (longer duration)             
delete_Cirsium_tracyi <- c(242472) # Retain 242471 (larger matrix)               
delete_Cypripedium_calceolus <- c(242609, 242624) # Retain 242581 (longer duration)       
delete_Dactylorhiza_lapponica <- c(242656) # Retain 242663 (longer duration)      
delete_Digitalis_purpurea <- c(242803) # Retain 242802 (longer duration)         
delete_Dracocephalum_austriacum <- c(238541) # Retain 242856 (longer duration)   
delete_Erythronium_japonicum <- c(238585, 243161) # Retain 243162 (larger matrix)       
delete_Euterpe_precatoria <- c(246326) # Retain 246328 (longer duration)          
delete_Hylocomium_splendens <- c(241108, 241109) # Retain 238673 (longer duration and larger matrix)        
delete_Lepanthes_rupestris <- c(244030) # Retain 244037 (longer duration)        
delete_Lomatium_bradshawii <- c(244255) # Retain 244207 (longer duration)         
delete_Mammillaria_huitzilopochtli <- c(247028) # Retain 247040 (larger matrix) 
delete_Mammillaria_pectinifera <- c(247064, 247065) # Retain 247066 (longer duration and larger matrix)     
delete_Neobuxbaumia_macrocephala <- c(247078) # Retain 247073 (same publication, equal in study duration and dimension)  
delete_Neobuxbaumia_mezcalaensis <- c(247082) # Retain 247086 (broader range of sizes examined)    
delete_Neobuxbaumia_tetetzo <- c(247093, 247097, 247098) # Retain 247102 (longer duration)       
delete_Panax_quinquefolius <- c(244577, 244585) # Retain 238925 (longer durationa and larger matrix)         
delete_Phyllanthus_emblica <- c(247740) # Retain 247692 (same publication, equal in study duration and dimension)         
delete_Pinus_lambertiana <- c(247799) # Retain 247789 (longer duration)          
delete_Plantago_coronopus <- c(244776) # Retain 244780 (longer duration)         
delete_Primula_veris <- c(245010, 245105) # Retain 245031 (longer duration)               
delete_Primula_vulgaris <- c(245197) # Retain 245132 (larger matrix)            
delete_Pseudophoenix_sargentii <- c(246364) # Retain 246371 (longer duration)    
delete_Pyrrocoma_radiata <- c(243393) # Retain 245290 (longer duration)          
delete_Rhododendron_maximum <- c(246735) # Retain 246736 (size-based, more accurate to measure, equal in study duration and dimension)        
delete_Rhododendron_ponticum <- c(247872) # Retain 247902 (larger matrix)      
delete_Salix_arctica <- c(246787) # Retain 246785 (larger matrix)              
delete_Sarracenia_purpurea <- c(245491) # Retain 245494 (longer duration)        
delete_Shorea_leprosula <- c(247953) # Retain 247949 (longer duration)           
delete_Solidago_altissima <- c(245613) # Retain 245611 (contains reproduction, equal in study duration and dimension)          
delete_Succisa_pratensis <- c(245715, 245727) # Retain 245693 (longer duration)           
delete_Taraxacum_campylodes <- c(245818) # Retain 245814 (larger matrix)        
delete_Tillandsia_multicaulis <- c(241194) # Retain 241198 (consistent state variable used, equal in study duration and dimension)     
delete_Tillandsia_punctulata <- c(241200) # Retain 241203 (consistent state variable used, equal in study duration and dimension)      
delete_Trillium_grandiflorum <- c(240174) # Retain 245986 (larger matrix)       
delete_Trollius_laxus <- c(246097) # Retain 246109 (longer duration)              
delete_Vella_pseudocytisus_subsp._paui <- c(246816, 246839) # Retain 246829 (longer duration)

# Delete selected duplicates
deleteDuplicateStudies <- c(
  delete_Actaea_elata,           
  delete_Alliaria_petiolata,          
  delete_Anthyllis_vulneraria,          
  delete_Araucaria_hunsteinii,          
  delete_Ascophyllum_nodosum,           
  delete_Astragalus_scaphoides,         
  delete_Astrophytum_capricorne,        
  delete_Astrophytum_ornatum,           
  delete_Calochortus_lyallii,           
  delete_Carduus_nutans,                
  delete_Carlina_vulgaris,              
  delete_Chamaedorea_radicalis,         
  delete_Cirsium_perplexans,            
  delete_Cirsium_pitcheri,              
  delete_Cirsium_tracyi,                
  delete_Cypripedium_calceolus,         
  delete_Dactylorhiza_lapponica,        
  delete_Digitalis_purpurea,            
  delete_Dracocephalum_austriacum,      
  delete_Erythronium_japonicum,         
  delete_Euterpe_precatoria,            
  delete_Hylocomium_splendens,          
  delete_Lepanthes_rupestris,           
  delete_Lomatium_bradshawii,           
  delete_Mammillaria_huitzilopochtli,   
  delete_Mammillaria_pectinifera,       
  delete_Neobuxbaumia_macrocephala,     
  delete_Neobuxbaumia_mezcalaensis,     
  delete_Neobuxbaumia_tetetzo,          
  delete_Panax_quinquefolius,           
  delete_Phyllanthus_emblica,           
  delete_Pinus_lambertiana,             
  delete_Plantago_coronopus,            
  delete_Primula_veris,                 
  delete_Primula_vulgaris,              
  delete_Pseudophoenix_sargentii,       
  delete_Pyrrocoma_radiata,             
  delete_Rhododendron_maximum,          
  delete_Rhododendron_ponticum,        
  delete_Salix_arctica,                 
  delete_Sarracenia_purpurea,           
  delete_Shorea_leprosula,              
  delete_Solidago_altissima,            
  delete_Succisa_pratensis,             
  delete_Taraxacum_campylodes,          
  delete_Tillandsia_multicaulis,        
  delete_Tillandsia_punctulata,         
  delete_Trillium_grandiflorum,         
  delete_Trollius_laxus,                
  delete_Vella_pseudocytisus_subsp._paui)

# Delete duplicate studies of suboptimal quality
dim(compadre)
compadre <- subset(compadre, !compadre$MatrixID %in% deleteDuplicateStudies)
dim(compadre)

# Make sure there are no duplicates left in the dataset
tableSpeciesAccepted <- table(compadre$SpeciesAccepted)
matrixSpeciesAccepted <- as.matrix(tableSpeciesAccepted)
speciesDuplicate <- rownames(matrixSpeciesAccepted)[which(matrixSpeciesAccepted>1)]
speciesDuplicate

# Make sure all single species MPMs are sensible for this study via visual inspection
for (i in 1:dim(compadre)[1]){
  try(subset <- compadre[i, c("MatrixID","SpeciesAccepted","CommonName","StudyDuration","MatrixComposite","MatrixDimension")])
  print(subset)
  print(compadre$mat[[i]]@matrixClass)
  print(compadre$mat[[i]]@matU)
  print(compadre$mat[[i]]@matF)
  print(compadre$mat[[i]]@matC)
  readline(prompt="Press [enter] to continue")
}

# Check models have reproduction, are not two-sex models and do not have any stages which do not make sense
# 238578 flowering individuals don't reproduce
# 238586 no reproduction
# 238846 flowering shoots don't reproduce
# 238856 no reproduction
# 240738 two-sex model
# 240808 no reproduction
# 240842 just seeds
# 240904 no reproduction
# 240956 no reproduction
# 241606 two-sex model
# 241949 only diagonal matrix (no movement between classes)
# 242040 reproductive stages not actually reproductive
# 243176 seedlings reproducing...
# 243309 multi-sex model
# 243357 multi-sex model
# 244284 no reproduction
# 245418 flowering individuals don't reproduce
# 245646 plants measured in agglomerations, reproduction and survival equated
# 245893 no reproduction
# 246049 no reproduction
# 246393 two-sex model
# 246613 two-sex model
# 247442 two-sex model
# 247551 no reproduction

# Remove models which do not fulfil these criteria

deleteNAStudies <- c(238578, 238586, 238846, 238856, 240738,
                     240808, 240842, 240904, 240956, 241606,
                     241949, 242040, 243176, 243309, 243357,
                     244284, 245418, 245646, 245893, 246049,
                     246393, 246613, 247442, 247551)

#Delete non-applicable studies
dim(compadre)
compadre <- subset(compadre, !compadre$MatrixID %in% deleteNAStudies)
dim(compadre)

# Distributions of matrix dimensions
table(compadre$dim)

# Search for parasitic plant families from Nickrent et al., 2020
parasitic_fams <- c("Convolvulaceae",
                    "Orobanchaceae",
                    "Lennoaceae",
                    "Mitrastemonaceae",
                    "Cytinaceae",
                    "Apodanthaceae",
                    "Rafflesiaceae",
                    "Krameriaceae",
                    "Cynomoriaceae",
                    "Erythropalaceae",
                    "Strombosiaceae",
                    "Octoknemaceae",
                    "Coulaceae",
                    "Ximeniaceae",
                    "Aptandraceae",
                    "Olacaceae",
                    "Balanophoraceae",
                    "Misodendraceae",
                    "Schoepfiaceae",
                    "Mystropetalaceae",
                    "Loranthaceae",
                    "Opiliaceae",
                    "Comandraceae",
                    "Thesiaceae",
                    "Cervantesiaceae",
                    "Nanodeaceae",
                    "Santalaceae",
                    "Amphorogynaceae",
                    "Viscaceae",
                    "Hydnoraceae",
                    "Lauraceae")

# Find species which are in these families
compadre[which(compadre$Family %in% parasitic_fams)]$SpeciesAccepted

# Kunkeliella subsucculenta is also a hemiparasite
# Pedicularis furbishiae is also a hemiparasite
# Others are non-parasitic

# Number of MPMs for which to calculate life history traits
long <- dim(compadre)[1]

# Dataframe to save life history trait outputs
output <- data.frame("unique.number"=rep(NA,long),
                     "SpeciesAuthor"=rep(NA,long),
                     "SpeciesAccepted"=rep(NA,long),
                     "CommonName"=rep(NA,long),
                     "Family"=rep(NA,long),
                     "Class"=rep(NA,long),
                     "Kingdom"=rep(NA,long),
                     "Authors"=rep(NA,long),
                     "Journal"=rep(NA,long),
                     "SourceType"=rep(NA,long),
                     "YearPublication"=rep(NA,long),
                     "DOI_ISBN"=rep(NA,long),
                     "Country"=rep(NA,long),
                     "Continent"=rep(NA,long),
                     "Ecoregion"=rep(NA,long),
                     "Population"=rep(NA,long),
                     "StartYear"=rep(NA,long),
                     "StartMonth"=rep(NA,long),
                     "StartSeason"=rep(NA,long),
                     "EndYear"=rep(NA,long),
                     "EndMonth"=rep(NA,long),
                     "EndSeason"=rep(NA,long),
                     "Treatment"=rep(NA,long),
                     "MatrixComposite"=rep(NA,long),
                     "Lat"=rep(NA,long),
                     "Lon"=rep(NA,long),
                     
                     "MatrixSplit"=rep(NA,long),
                     "MatrixFec"=rep(NA,long),
                     "MatrixDimension"=rep(NA,long),
                     
                     "Ergodic"=rep(NA,long),
                     "Irreducible"=rep(NA,long),
                     "Primitive"=rep(NA,long),
                     "notProp"=rep(NA,long),
                     
                     "Lambda"=rep(NA,long),
                     "GenTfun"=rep(NA,long),
                     "R0"=rep(NA,long),
                     "Lmean"=rep(NA,long),
                     "La"=rep(NA,long)
)


# Convert projection interval to numeric and assuming an annual basis when not stated
compadre$ProjectionInterval <- as.numeric(compadre$ProjectionInterval)
compadre$ProjectionInterval[which(is.na(compadre$ProjectionInterval))] <- 1

# Plot the frequency of projection interval
hist((as.numeric(compadre$ProjectionInterval)), xlab= "Projection interval (years)", col="blue", main ="", cex.lab=1.2, breaks=200, xlim=c(0,max(compadre$ProjectionInterval)), ylim=c(0,140))

# Run through all models and calculate life history traits
d <- compadre
count <- 0

for (i in 1:long){
  count <- count + 1
  
  # Transfer metadata over from the metadata file
  output[count,c("MatrixID", "SpeciesAuthor","SpeciesAccepted","CommonName","Family","Class","Kingdom","Authors","Journal","SourceType","YearPublication","DOI_ISBN","Country","Continent","Ecoregion","Lat","Lon","StartYear","StartSeason","StartMonth","EndYear","EndSeason","EndMonth","Population","Treatment","MatrixComposite","MatrixSplit","MatrixFec")]=
    d@data[i,c("MatrixID", "SpeciesAuthor","SpeciesAccepted","CommonName","Family","Class","Kingdom","Authors","Journal","SourceType","YearPublication","DOI_ISBN","Country","Continent","Ecoregion","Lat","Lon","MatrixStartYear","MatrixStartSeason","MatrixStartMonth","MatrixEndYear","MatrixEndSeason","MatrixEndMonth","MatrixPopulation","MatrixTreatment","MatrixComposite","MatrixSplit","MatrixFec")]
  
  # Define the beginning of life when an individual become established. Thus, we do not consider transitions from the "prop" stages
  try(lifeStages <- d$mat[[i]]@matrixClass$MatrixClassOrganized)
  try(prop <- (which(lifeStages == "prop")))
  try(active <- (which(lifeStages == "active")))
  try(dorm <- (which(lifeStages == "dorm")))
  try(matU <- as.matrix(d$mat[[i]]@matU)^(1/d$ProjectionInterval[i]))
  try(matU[is.na(matU)] <- 0)
  try(matF <- as.matrix(d$mat[[i]]@matF)^(1/d$ProjectionInterval[i]))
  try(matF[is.na(matF)] <- 0)
  try(matC <- as.matrix(d$mat[[i]]@matC)^(1/d$ProjectionInterval[i]))
  try(matC[is.na(matC)] <- 0)
  try(matA <- matU + matF + matC)
  
  #Re-arrange stages if necessary
  try(reArrangeVector <- c(prop, active, dorm))
  try(matU <- as.matrix(matU[reArrangeVector, reArrangeVector]))
  rownames(matU) <- colnames(matU)
  try(matF <- as.matrix(matF[reArrangeVector, reArrangeVector]))
  rownames(matF) <- colnames(matF)
  try(matC <- as.matrix(matC[reArrangeVector, reArrangeVector]))
  rownames(matC) <- colnames(matC)
  try(matA <- as.matrix(matA[reArrangeVector, reArrangeVector]))
  rownames(matA) <- colnames(matA)
  
  # Designate non-propagule stage
  try(output$notProp[count] <- notProp <- min(which(lifeStages != "prop")))
  try(if (is.infinite(notProp)) {notProp <- 1})
  
  #Matrix Dimension
  try(matDim <- dim(matU)[1])
  try(output$MatrixDimension[count] <- matDim)
  
  if(matDim>1){
    
    #Primitivity
    try(output$Primitive[count] <- isPrimitive(matA))
    
    #Ergodicity
    try(output$Ergodic[count] <- isErgodic(matA))
    
    #Irreducibility
    try(output$Irreducible[count] <- isIrreducible(matA))
    
    #Net reproductive rate (R0; Caswell 2001, p 126)"The mean numbe of offspring by which an individual will be replaced by the end of its live.
    try(output$R0[count] <- R0 <- net_repro_rate(matU, matF, notProp))
    
    #lambda
    try(output$Lambda[count] <- Lambda <- eigen.analysis(matA)$lambda1)
    
    #Generation time (T; Caswell 2001, p 129)
    try(output$GenTfun[count] <- abs(log(R0)/log(Lambda)))
    
    # Calculate mean life expectancy from mature distribution
    try(mat_dist <- mature_distrib(matU, start = notProp, repro_stages = repro_stages(matF)))
    try(output[count,"Lmean"] <- life_expect_mean(matU, mixdist = mat_dist, start = NULL))
    
    # Calculate mean age at maturity
    try(output[count,"La"] <- mature_age(matU, matF, notProp))
    
  }
  
  print(paste(i, "-", long, "-", d$SpeciesAuthor[i]))
}

# Estimate reproductive window from mean life expectancy and age at first reproduction
output$Lamean <- output$Lmean - output$La
hist(output$Lamean)
range(output$Lamean, na.rm=T)

# Outlier removal will be important

# Replace all infinite and NaN values with NA
output <- do.call(data.frame, lapply(output, function(x) {
  replace(x, is.infinite(x) | is.nan(x), NA)
}))

# Define all LHTs and demographic properties
allLHTs <- c("Lambda", "GenTfun", "R0", "Lmean", "La", "Lamean")

# Add mistletoe life history traits

# Read in mistletoe LHTs
Mst_LHTs <- read.csv("Mst_LHTs.csv")

# Add to output data frame
output <- full_join(output, Mst_LHTs)

# Create function to remove outliers - use 1.5 interquartlie range as not approximated as normal
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Remove all outliers
for (i in allLHTs){
  output[,i] <- remove_outliers(output[,i])
}

# Run through each LHT and check for normality
for (i in allLHTs){
  par(mfrow=c(2,2))
  hist(output[,i], main=i)
  try(hist(log(output[,i]), main="log"))
  try(hist(log(output[,i]) + 1, main="log+1"))
  try(hist(sqrt(output[,i]), main="sqrt"))
  readline(prompt = "Press [enter] to continue")
}

# Traits left unaltered:
# Lambda, R0

# Figure S10
# Traits to be log-transformed:
# GenTfun
output$GenTfun <- log(output$GenTfun)
hist(output$GenTfun, breaks = 10, main = "Generation time", xlab = "Log-transformed generation time (years)")

#Lmean
output$Lmean <- log(output$Lmean)

# Remove non-positive values and -Inf as due to calculation errors
output$Lmean <- replace(output$Lmean, is.infinite(output$Lmean), NA)
output$Lmean <- replace(output$Lmean, output$Lmean < 0, NA)
hist(output$Lmean, breaks = 10, main = "Mean life expectancy", xlab = "Log-transformed mean life expectancy (years)")

# La
output$La <- log(output$La)
output$La <- replace(output$La, output$La < 0, NA)
hist(output$La, breaks = 10, main = "Age at maturity", xlab = "Log-transformed age at maturity (years)")

# To log-transform Lamean, add minimum value + 1 first
output$Lamean <- log(output$Lamean + abs(min(output$Lamean, na.rm = TRUE)) + 1)
hist(output$Lamean, breaks = 10, main = "Reproductive window", xlab = "Transformed reproductive window (years)")

# Examine distribution of first non-propagule stages
hist(output$notProp)

# Designate parasites -----------------------------------------------------------

# Designate parasitic plants
output <- mutate(output, Parasite = case_when(SpeciesAccepted %in% c("Viscum album", 
                                                                     "Kunkeliella subsucculenta", 
                                                                     "Pedicularis furbishiae") ~ "Parasite", 
                                              .default = "Free-living"))

# Find papers for parasitic plants
subset(output, Parasite == "Parasite")[, c("Authors", "Journal", "SourceType", "YearPublication", "DOI_ISBN")]

# Phylogenetic data ----------------------------------------------------------

# Match names in output to Open Tree Taxonomy
output_species_resolved <- tnrs_match_names(output$SpeciesAccepted)

# Correct subsp. and var. to only species names
output$SpeciesAccepted[which(output$SpeciesAccepted == "Arenaria grandiflora subsp. bolosii")]        <- "Arenaria grandiflora"
output$SpeciesAccepted[which(output$SpeciesAccepted == "Chamaecrista lineata var. keyensis")]         <- "Chamaecrista lineata"
output$SpeciesAccepted[which(output$SpeciesAccepted == "Echinospartum ibericum subsp. algibicum")]    <- "Echinospartum ibericum"
output$SpeciesAccepted[which(output$SpeciesAccepted == "Eriogonum longifolium var. gnaphalifolium")]  <- "Eriogonum longifolium"
output$SpeciesAccepted[which(output$SpeciesAccepted == "Oenothera coloradensis subsp. coloradensis")] <- "Oenothera coloradensis"
output$SpeciesAccepted[which(output$SpeciesAccepted == "Lespedeza juncea var. sericea")]              <- "Lespedeza juncea"
output$SpeciesAccepted[which(output$SpeciesAccepted == "Petrocoptis pyrenaica subsp. pseudoviscosa")] <- "Petrocoptis pyrenaica"

# Re-attempt to match names
output_species_resolved <- tnrs_match_names(output$SpeciesAccepted)

# Create tree using names from Open Tree Taxonomy
try(tree <- tol_induced_subtree(ott_ids = output_species_resolved$ott_id))

# List ott_ids not in TOL
invalid_ott_ids <- c(16952, 3955000, 5144270, 5525303, 5738006, 807307)

# Remove corresponding species from output
invalid_sp <- subset(output_species_resolved, output_species_resolved$ott_id %in% invalid_ott_ids)$search_string %>%
  str_to_sentence()
output <- subset(output, !output$SpeciesAccepted %in% invalid_sp)

# Remove ott_ids not in list
output_species_resolved <- subset(output_species_resolved, !output_species_resolved$ott_id %in% invalid_ott_ids)

# Rename SpeciesAccepted names to those in phylogeny
# Capitalise search string in resolved list
output_species_resolved$search_string <- str_to_sentence(output_species_resolved$search_string)

# Rename row names to search string
rownames(output_species_resolved) <- output_species_resolved$search_string

# Rename rows to accepted names
rownames(output) <- output$SpeciesAccepted

# Select species which have update names
mismatched_sp <- subset(output, output$SpeciesAccepted %in% output_species_resolved$unique_name == FALSE)$SpeciesAccepted

# Rename to accepted names
for(species in mismatched_sp){
  output[species,]$SpeciesAccepted <- output_species_resolved[species, ]$unique_name
}

# Check mismatched species have been resolved
subset(output, output$SpeciesAccepted %in% output_species_resolved$unique_name == FALSE)$SpeciesAccepted

# Re-attempt to create tree using names from Open Tree Taxonomy
tree <- tol_induced_subtree(ott_ids = output_species_resolved$ott_id)

# Check for missing labels
output[-which(output$SpeciesAccepted %in% tree$tip.label),]$SpeciesAccepted

# Strip OTT IDs
tree$tip.label <- strip_ott_ids(tree$tip.label, remove_underscores = T)
tree$node.label <- NULL

# Remove extra brackets
tree$tip.label <- gsub(" \\(species in kingdom Archaeplastida\\)", "", tree$tip.label)

# Add edge lengths
tree <- compute.brlen(tree)
tree$edge.length

# Check if tree is ultrametric
is.ultrametric(tree)

# Plot tree
par(mfrow = c(1, 1))
plot(tree, type = 'radial', cex = 0.4)
plot(tree, cex = 0.4)

# Form phylogenetic object
row.names(output) <- output$SpeciesAccepted
obj <- name.check(tree, output)
obj

tree2 <- tree
tree2$tip.label <- paste("    ", tree2$tip.label, "    ")

# Create .pdf of tree
pdf("Tree.pdf", h = 8, w = 15)
plotTree(ladderize(tree2),type="arc",fsize=.7, depth=2,
         lwd=1,ftype="i",arc_height=.2)
dev.off()

# Assessing "missingness" of data -------------------------------------------

# Perform pair-wise correlations of LHTs to evaluate co-linearities

# Define LHTs to test correlations between
LHTs <- c("GenTfun","Lmean","La","Lamean")

LHTSymbols <-c(expression("T"),
               expression(eta["e"]),
               expression("L"[alpha]),
               expression(paste("L"[paste(alpha,"-",omega)])))


# Function to add correlation coefficients
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y, use = "complete.obs", method = "spearman")) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

# Define diagonal panel of plot
diag.panel <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  
  i <- which(sapply(output[, LHTs], identical, x))
  expr <- LHTSymbols[i]
  
  # Draw only the expression
  text(0.5, 0.5, labels = expr, cex = 6, font = 1)
}


# Plotting the correlation matrix
Fig_S14 <- pairs(output[, LHTs],
                 labels = rep("", length(LHTs)),
                 upper.panel = panel.cor,
                 lower.panel = panel.smooth,
                 diag.panel = diag.panel)
Fig_S14

# Correlation in absolute terms
corr <- abs(cor(output[, LHTs], use = "complete.obs"))

colors <- dmat.color(corr)
#order <- order.single(corr)

cpairs(output[, LHTs],                    # Data frame of variables
       #       order,                   # Order of the variables
       panel.colors = colors,   # Matrix of panel colors
       border.color = alpha("grey70",.2), # Borders color
       gap = 0.45,              # Distance between subplots
       main = "Ordered variables colored by correlation", # Main title
       show.points = TRUE,      # If FALSE, removes all the points
       pch = 21,                # pch symbol
       bg = alpha('blue',.2)) # Colors by group

# Plotting the correlation matrix without imputations
pdf("Pair-wise LHT correlations no imputations.pdf")
pairs(output[, LHTs],
      labels = rep("", length(LHTs)),
      upper.panel = panel.cor,
      lower.panel = panel.smooth,
      diag.panel = diag.panel)
dev.off()

# Calculate percentage of missing data
pct_miss(output[, LHTs])

pdf("Missing LHTs.pdf")
gg_miss_var(output[, LHTs], show_pct = TRUE)
dev.off()

pdf("Missing LHTs visualised.pdf")
vis_miss(output[, LHTs])
dev.off()

# Impute missing data -----------------------------------------------------

unimputedOutput <- output

# Imputing missing data
imputedOutput <- imputedOutput1 <- imputedOutput2 <- imputedOutput3 <- imputedOutput4 <- imputedOutput5 <- 
  imputedOutput6 <- imputedOutput7 <- imputedOutput8 <- imputedOutput9 <- imputedOutput10 <- output

set.seed(120)
imputedLHTs1 <- mice(output[,LHTs], m=1)
set.seed(121)
imputedLHTs2 <- mice(output[,LHTs], m=1)
set.seed(122)
imputedLHTs3 <- mice(output[,LHTs], m=1)
set.seed(123)
imputedLHTs4 <- mice(output[,LHTs], m=1)
set.seed(124)
imputedLHTs5 <- mice(output[,LHTs], m=1)
set.seed(125)
imputedLHTs6 <- mice(output[,LHTs], m=1)
set.seed(126)
imputedLHTs7 <- mice(output[,LHTs], m=1)
set.seed(127)
imputedLHTs8 <- mice(output[,LHTs], m=1)
set.seed(128)
imputedLHTs9 <- mice(output[,LHTs], m=1)
set.seed(129)
imputedLHTs10 <- mice(output[,LHTs], m=1)

#imputedLHTs <- mice(output[,LHTs], m=1)
imputedOutput1[,LHTs] <- complete(imputedLHTs1, 1)
imputedOutput2[,LHTs] <- complete(imputedLHTs2, 1)
imputedOutput3[,LHTs] <- complete(imputedLHTs3, 1)
imputedOutput4[,LHTs] <- complete(imputedLHTs4, 1)
imputedOutput5[,LHTs] <- complete(imputedLHTs5, 1)
imputedOutput6[,LHTs] <- complete(imputedLHTs6, 1)
imputedOutput7[,LHTs] <- complete(imputedLHTs7, 1)
imputedOutput8[,LHTs] <- complete(imputedLHTs8, 1)
imputedOutput9[,LHTs] <- complete(imputedLHTs9, 1)
imputedOutput10[,LHTs] <- complete(imputedLHTs10, 1)

for (i in LHTs){
  for (j in 1:dim(imputedOutput)[1]){
    imputedOutput[j,i] <- mean(imputedOutput1[j,i],
                               imputedOutput2[j,i],
                               imputedOutput3[j,i],
                               imputedOutput4[j,i],
                               imputedOutput5[j,i],
                               imputedOutput6[j,i],
                               imputedOutput7[j,i],
                               imputedOutput8[j,i],
                               imputedOutput9[j,i],
                               imputedOutput10[j,i],
                               na.rm=T)
  }
}

#imputedOutput[,LHTs] <- complete(imputedLHTs, 1)

# Check whether NAs returned
imputedOutput[,LHTs] %>% 
  dplyr::summarise(across(everything(), ~ sum(is.na(.x))))

# No NAs returned

# Plotting the correlation matrix with imputations
pdf("Pair-wise LHT correlations with imputations.pdf")
pairs(imputedOutput[, LHTs],
      labels = LHTSymbols,
      upper.panel = panel.cor,    # Correlation panel
      lower.panel = panel.smooth) # Smoothed regression lines
dev.off()

# PCA --------------------------------------------------------

# Check how many species will be used in the PCA
nrow(imputedOutput)

# Check how many parasites and free-living species
table(imputedOutput$Parasite)

#Traits used
LHTs <- c("GenTfun", "Lmean", "La", "Lamean")
# LHTSymbols <- c("T",
#                 expression(eta["e"]),
#                 expression("L"[alpha]),
#                 expression(paste("L"[paste(alpha,"-",omega)])))
# LHTcols <- c("black", "darkgreen", "blue", "purple")
# 
# 
#PCA without phylogeny
pca <- prcomp(imputedOutput[, LHTs], center = T, scale. = T)
row.names(pca$x) <- imputedOutput$SpeciesAccepted
variancePCA <- summary(pca)$importance[2,]
variancePCA

#Invert order on PC1 for interpretability
pca$x[,"PC1"] <- - pca$x[,"PC1"]
pca$rotation[,"PC1"] <- - pca$rotation[,"PC1"]

# par(mfrow = c(1,1))
# # Plot PCA
#   plot(pca$x[,"PC1"], pca$x[,"PC2"], 
#        col=alpha("black",0.5), pch=21, ylim=c(-4, 5), xlim=c(-4, 4),
#        xlab=paste0("PC1 (",round(variancePCA[1]*100,2),"%)"),ylab=paste0("PC2 (",round(variancePCA[2]*100,2),"%)"))
#   arrows(0,0, pca$rotation[,"PC1"]*5, pca$rotation[,"PC2"]*5, lwd=6, col=LHTcols)
#   arrows(0,0, pca$rotation[,"PC1"]*5, pca$rotation[,"PC2"]*5, lwd=3, col="white")
#   points(pca$x["Viscum album","PC1"], pca$x["Viscum album","PC2"], col = "red", pch = 16)
#   points(pca$x["Thesium subsucculentum","PC1"], pca$x["Thesium subsucculentum","PC2"], col = "darkorange", pch = 16)
#   points(pca$x["Pedicularis furbishiae","PC1"], pca$x["Pedicularis furbishiae","PC2"], col = "brown", pch = 16)
#   text(pca$rotation[,"PC1"]*6, pca$rotation[,"PC2"]*6, LHTSymbols, col=LHTcols,cex=1.2)

# Plot PCA with species names
# Assuming pca is the result of `prcomp()`
# Create a data frame for the PCA results (with x and rotation components)
pca_data <- as.data.frame(pca$x)
pca_data$Species <- row.names(pca_data)
pca_rotations <- as.data.frame(pca$rotation)

# Get UUIDs from scientific names for particular species
viscum_uuid <- "e354b2f8-cda4-4fdb-9040-4b01cdb9eaeb"
pedicularis_uuid <- get_uuid("Pedicularis")
helianthus_uuid <- get_uuid("Helianthus")
brassica_uuid <- get_uuid("Brassica rapa")
abies_uuid <- get_uuid("Abies magnifica")
sequoia_uuid <- "51270d2b-6d63-43fc-b899-a9ff7a338d6c"
escontria_uuid <- "b95ba2c8-8b80-4346-acd9-b27ea03622f9"
calocedrus_uuid <- "db046be5-7aad-4889-b55c-6c0052dea91a"
quercus_uuid <- "b34746dd-ec23-4fa0-83ed-013f1d9dd6be"

# Create silhouette data frame
sil_df <- data.frame(
  PC1 = c(subset(pca_data, rownames(pca_data) == "Viscum album")$PC1,
          subset(pca_data, rownames(pca_data) == "Pedicularis furbishiae")$PC1,
          subset(pca_data, rownames(pca_data) == "Helianthus divaricatus")$PC1,
          subset(pca_data, rownames(pca_data) == "Brassica napus")$PC1,
          subset(pca_data, rownames(pca_data) == "Abies magnifica")$PC1,
          subset(pca_data, rownames(pca_data) == "Sequoia sempervirens")$PC1,
          subset(pca_data, rownames(pca_data) == "Escontria chiotilla")$PC1*1.2,
          subset(pca_data, rownames(pca_data) == "Calocedrus macrolepis")$PC1),
  PC2 = c(subset(pca_data, rownames(pca_data) == "Viscum album")$PC2*5,
          subset(pca_data, rownames(pca_data) == "Pedicularis furbishiae")$PC2*6,
          subset(pca_data, rownames(pca_data) == "Helianthus divaricatus")$PC2*3.9,
          subset(pca_data, rownames(pca_data) == "Brassica napus")$PC2,
          subset(pca_data, rownames(pca_data) == "Abies magnifica")$PC2,
          subset(pca_data, rownames(pca_data) == "Sequoia sempervirens")$PC2*1.2,
          subset(pca_data, rownames(pca_data) == "Escontria chiotilla")$PC2*0.8,
          subset(pca_data, rownames(pca_data) == "Calocedrus macrolepis")$PC1),
  uuid = c(viscum_uuid,
           pedicularis_uuid,
           helianthus_uuid,
           brassica_uuid,
           abies_uuid,
           sequoia_uuid,
           escontria_uuid,
           calocedrus_uuid),
  col = c(colours(6)[5],
          colours(6)[5],
          "lightgrey",
          "lightgrey",
          "lightgrey",
          "lightgrey",
          "lightgrey",
          "lightgrey"))

# Create an expression vector
LHTSymbols <- list(
  expression(T),
  expression(eta["e"]),
  expression("L"[alpha]),
  expression(paste("L"[paste(alpha, "-", omega)]))
)

# Add symbols to the PCA rotation data frame
pca_rotations$Symbols <- LHTSymbols

# Extract variance explained for the first two components (assuming this is calculated before)
variancePCA <- (pca$sdev^2) / sum(pca$sdev^2)

# Plot PCA
PCA_plot <- ggplot(data = pca_data, aes(x = PC1, y = PC2)) +
  # Points for the PCA plot
  geom_point(color = alpha("black", 0.3), shape = 21, size = 2.3) +
  
  # Overlay images
  add_phylopic(uuid = sil_df$uuid,
               x = sil_df$PC1*1.2, 
               y = sil_df$PC2*1.2,
               alpha = 0.8, 
               height = 1.5,
               fill = sil_df$col) +
  
  # Arrows for PCA loadings (rotation)
  geom_segment(data = pca_rotations, 
               aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")), 
               size = 1.5, color = "black") +
  
  # Arrows for PCA loadings (white outline)
  geom_segment(data = pca_rotations, 
               aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")), 
               size = 0.8, color = "white") +
  
  # Points for specific species and their labels
  geom_point(data = subset(pca_data, rownames(pca_data) == "Viscum album"), 
             aes(x = PC1, y = PC2), color = colours(6)[5], size = 3) +
  geom_text(data = subset(pca_data, rownames(pca_data) == "Viscum album"), 
            aes(x = PC1 * 1.4, y = PC2 * 1.4, label = "V. album"), 
            color = colours(6)[5], size = 5) +
  
  geom_point(data = subset(pca_data, rownames(pca_data) == "Thesium subsucculentum"), 
             aes(x = PC1, y = PC2), color = colours(6)[5], size = 3) +
  geom_text(data = subset(pca_data, rownames(pca_data) == "Thesium subsucculentum"), 
            aes(x = PC1 * 2.5, y = PC2 * -1.5, label = "T. subsucculentum"), 
            color = colours(6)[5], size = 5) +
  
  geom_point(data = subset(pca_data, rownames(pca_data) == "Pedicularis furbishiae"), 
             aes(x = PC1, y = PC2), color = colours(6)[5], size = 3) +
  geom_text(data = subset(pca_data, rownames(pca_data) == "Pedicularis furbishiae"), 
            aes(x = PC1 * 1.4, y = PC2 * 2.3, label = "P. furbishiae"), 
            color = colours(6)[5], size = 5) +
  
  
  # Axis labels and plot title with variance explained
  labs(x = paste0("PC1 - Fast-slow continuum (", round(variancePCA[1] * 100, 2), "%)"),
       y = paste0("PC2 - Reproductive strategy (", round(variancePCA[2] * 100, 2), "%)")) +
  
  # Customize theme for cleaner visuals
  theme_bw() +
  theme(
    legend.position = "none",  # Remove legend (not needed in this case)
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )

# Manually add each expression using annotate()
for (i in seq_along(LHTSymbols)) {
  PCA_plot <- PCA_plot + annotate(
    "text",
    x = pca_rotations$PC1[i] * 5.5,
    y = pca_rotations$PC2[i] * 5.5,
    label = LHTSymbols[[i]],
    color = "black",
    size = 6,
    parse = TRUE
  )
}


# PCA with phylogeny
# Scale values because phyl.pca does not do it, and otherwise PCA is very stretched out on PC1
outputPhyl <- imputedOutput
for (i in LHTs){
  outputPhyl[,i]  <- scale(outputPhyl[,i], center= T, scale = T)
}

rownames(outputPhyl) <- outputPhyl$SpeciesAccepted

# Remove rows not in tree
outputPhyl <- subset(outputPhyl, SpeciesAccepted %in% tree$tip.label)

# Perform pPCA
phyloPCA <- phyl.pca(tree, outputPhyl[, LHTs], method = "lambda")
phyloPCA$lambda
# Pagel's λ = 0.101

# Percentage variance explained
variancePhyloPCA <- diag(phyloPCA$Eval)/sum(phyloPCA$Eval)
variancePhyloPCA

#Invert order on PC1 for interpretability
phyloPCA$S[,"PC1"] <- -phyloPCA$S[,"PC1"]
phyloPCA$L[,"PC1"] <- -phyloPCA$L[,"PC1"]

phylo_pca_data <- data.frame(PC1 = phyloPCA$S[,"PC1"], 
                             PC2 = phyloPCA$S[,"PC2"])
phylo_pca_data$Species <- row.names(phylo_pca_data)

# Create rotations data frame
phylo_pca_rotations <- data.frame(PC1 = phyloPCA$L[,"PC1"],
                                  PC2 = phyloPCA$L[,"PC2"])

# Get UUIDs from scientific names for particular species
viscum_uuid <- "e354b2f8-cda4-4fdb-9040-4b01cdb9eaeb"
pedicularis_uuid <- get_uuid("Pedicularis")
helianthus_uuid <- get_uuid("Helianthus")
brassica_uuid <- get_uuid("Brassica rapa")
abies_uuid <- get_uuid("Abies magnifica")
sequoia_uuid <- "51270d2b-6d63-43fc-b899-a9ff7a338d6c"
escontria_uuid <- "b95ba2c8-8b80-4346-acd9-b27ea03622f9"
calocedrus_uuid <- "db046be5-7aad-4889-b55c-6c0052dea91a"
quercus_uuid <- "b34746dd-ec23-4fa0-83ed-013f1d9dd6be"
araucaria_uuid <- "83feef9b-4096-4710-afa7-d4a978cfaa93"

# Create silhouette data frame
sil_df <- data.frame(
  PC1 = c(subset(phylo_pca_data, rownames(phylo_pca_data) == "Viscum album")$PC1,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Pedicularis furbishiae")$PC1,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Helianthus divaricatus")$PC1*1.2,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Brassica napus")$PC1,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Abies magnifica")$PC1,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Calocedrus macrolepis")$PC1*1.1,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Sequoia sempervirens")$PC1,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Escontria chiotilla")$PC1*1.2,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Quercus mongolica subsp. crispula")$PC1*2,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Araucaria cunninghamii")$PC1),
  PC2 = c(subset(phylo_pca_data, rownames(phylo_pca_data) == "Viscum album")$PC2*5,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Pedicularis furbishiae")$PC2*6.5,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Helianthus divaricatus")$PC2*4.5,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Brassica napus")$PC2*(-3),
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Abies magnifica")$PC2,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Calocedrus macrolepis")$PC2*1.1,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Sequoia sempervirens")$PC2*1.2,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Escontria chiotilla")$PC2*0.8,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Quercus mongolica subsp. crispula")$PC2*1.1,
          subset(phylo_pca_data, rownames(phylo_pca_data) == "Araucaria cunninghamii")$PC2*1.2),
  uuid = c(viscum_uuid,
           pedicularis_uuid,
           helianthus_uuid,
           brassica_uuid,
           abies_uuid,
           calocedrus_uuid,
           sequoia_uuid,
           escontria_uuid,
           quercus_uuid,
           araucaria_uuid),
  col = c(colours(6)[5],
          colours(6)[5],
          "lightgrey",
          "lightgrey",
          "lightgrey",
          "lightgrey",
          "lightgrey",
          "lightgrey",
          "lightgrey",
          "lightgrey"))

# Read in Thesium custom silhouette
img <- readPNG("Thesium_silhouette.png")
grob_img <- rasterGrob(img, interpolate = TRUE)

# Plot pPCA
Fig_5 <- ggplot(data = phylo_pca_data, aes(x = PC1, y = PC2)) +
  # Points for the PCA plot
  geom_point(color = alpha("black", 0.3), shape = 21, size = 2.3) +
  
  # Overlay images
  add_phylopic(uuid = sil_df$uuid,
               x = sil_df$PC1*1.2,
               y = sil_df$PC2*1.2,
               alpha = 0.8,
               height = 1.5,
               fill = sil_df$col) +
  
  # Add custom Thesium silhouette
  annotation_custom(
    grob = grob_img,
    xmin = -1.5, 
    xmax = -1,  
    ymin = 1, 
    ymax = 2
  ) +
  
  # Arrows for PCA loadings (rotation)
  geom_segment(data = phylo_pca_rotations, 
               aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")), 
               size = 1.5, color = "black") +
  
  # Arrows for PCA loadings (white outline)
  geom_segment(data = phylo_pca_rotations, 
               aes(x = 0, y = 0, xend = PC1 * 5, yend = PC2 * 5), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")), 
               size = 0.8, color = "white") +
  
  # Points for specific species and their labels
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Sequoia sempervirens"), 
             aes(x = PC1, y = PC2), color = "grey", size = 3) +
  
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Calocedrus macrolepis"),
             aes(x = PC1, y = PC2), color = "grey", size = 3) +
  
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Helianthus divaricatus"),
             aes(x = PC1, y = PC2), color = "grey", size = 3) +
  
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Abies magnifica"),
             aes(x = PC1, y = PC2), color = "grey", size = 3) +
  
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Escontria chiotilla"),
             aes(x = PC1, y = PC2), color = "grey", size = 3) +
  
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Quercus mongolica subsp. crispula"),
             aes(x = PC1, y = PC2), color = "grey", size = 3) +
  
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Brassica napus"),
             aes(x = PC1, y = PC2), color = "grey", size = 3) +
  
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Araucaria cunninghamii"),
             aes(x = PC1, y = PC2), color = "grey", size = 3) +
  
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Viscum album"), 
             aes(x = PC1, y = PC2), color = colours(6)[5], size = 3) +
  geom_text(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Viscum album"), 
            aes(x = PC1 * 1.4, y = PC2 * 1.9, label = "italic('V. album')"), parse = TRUE, 
            color = colours(6)[5], size = 5) +
  
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Thesium subsucculentum"), 
             aes(x = PC1, y = PC2), color = colours(6)[5], size = 3) +
  geom_text(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Thesium subsucculentum"), 
            aes(x = PC1 * 2.5, y = PC2 * -1.4, label = "italic('T. subsucculentum')"), parse = TRUE, 
            color = colours(6)[5], size = 5) +
  
  geom_point(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Pedicularis furbishiae"), 
             aes(x = PC1, y = PC2), color = colours(6)[5], size = 3) +
  geom_text(data = subset(phylo_pca_data, rownames(phylo_pca_data) == "Pedicularis furbishiae"), 
            aes(x = PC1 * 0.3, y = PC2 * 3.2, label = "italic('P. furbishiae')"), parse = TRUE, 
            color = colours(6)[5], size = 5) +
  
  
  # Axis labels and plot title with variance explained
  labs(x = paste0("PC1 - Fast-slow continuum (", round(variancePhyloPCA[1] * 100, 2), "%)"),
       y = paste0("PC2 - Reproductive strategy (", round(variancePhyloPCA[2] * 100, 2), "%)")) +
  
  # Customize theme for cleaner visuals
  theme_bw() +
  theme(
    legend.position = "none",  # Remove legend (not needed in this case)
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  )

# Manually add each expression using annotate()
for (i in seq_along(LHTSymbols)) {
  Fig_5 <- Fig_5 + annotate(
    "text",
    x = phyloPCA$L[,"PC1"][i]*5.5,
    y = phyloPCA$L[,"PC2"][i]*5.5,
    label = LHTSymbols[[i]],
    color = "black",
    size = 6,
    parse = TRUE
  )
}

Fig_5
ggsave("Figure 5.png", Fig_5)

