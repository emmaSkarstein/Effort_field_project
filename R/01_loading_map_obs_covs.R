
# This script loads the observations, map and covariates into the environment.


library(ggplot2)
library(sf) # for sf objects
library(sp)
library(plyr)
library(dplyr) # for smoother dataframe-manipulation
library(ggmap) # also for nice maps
library(maps)
library(PointedSDMs)
library(spatstat)
library(maptools)
library(INLA)
library(rgeos)
library(fields)
library(viridis)


# MAP ---------------------------------------------------------------
# This is always the same
norway <- ggplot2::map_data("world", region = "Norway(?!:Svalbard)")
norway <- dplyr::setdiff(norway, filter(norway, subregion == "Jan Mayen"))
Projection <- CRS("+proj=longlat +ellps=WGS84")
norwayfill <- map("world", "norway", fill=TRUE, plot=FALSE, 
                  ylim=c(58,72), xlim=c(4,32))
IDs <- sapply(strsplit(norwayfill$names, ":"), function(x) x[1])
norway.poly <- map2SpatialPolygons(norwayfill, IDs = IDs, 
                                   proj4string = Projection)


# LOADING DATA AND COVARIATES ---------------------------------------------------------------------------------

# Covariates
# This section varies based on which covariates are used
covariateData <- readRDS("data/environmental_covariates.RDS")
covariateData <- covariateData[complete.cases(covariateData$decimalLatitude,covariateData$decimalLongitude,covariateData$area_km2,covariateData$HFP),]
covariateData <- covariateData %>% 
  mutate(log_area = log(area_km2), 
         log_perimeter = log(perimeter_m), 
         log_catchment = log(catchment_area_km2)) %>% 
  dplyr::select(-c(ebint, no_vatn_lnr, eb_waterregionID))
env_covariateData <- covariateData %>% dplyr::select(-c(distance_to_road, HFP))


Covariates <- SpatialPointsDataFrame(coords = covariateData[,c("decimalLongitude","decimalLatitude")],
                                     data = covariateData, 
                                     proj4string = Projection)
env_covariates <- SpatialPointsDataFrame(coords = env_covariateData[,c("decimalLongitude","decimalLatitude")],
                                         data = env_covariateData, 
                                         proj4string = Projection)
Covariates@data <- data.frame(apply(Covariates@data, 2, scale))  # scale the covariates 
env_covariates <- data.frame(apply(env_covariates@data, 2, scale))  # scale the covariates 

# Observations
# The data is always the same, but will possibly be split into different sections for validation
Data_survey_df <- readRDS("data/survey_clean.rds")
Data_survey <- SpatialPointsDataFrame(coords = Data_survey_df[,c("decimalLongitude","decimalLatitude")], 
                                      data = Data_survey_df[,c("occurrenceStatus","species")],
                                      proj4string = Projection)

Data_artsobs_df <- readRDS("data/artsobs_clean.rds")
Data_artsobs <- SpatialPointsDataFrame(coords = Data_artsobs_df[,c("decimalLongitude","decimalLatitude")], 
                                       data = Data_artsobs_df[,c("occurrenceStatus","species")],
                                       proj4string = Projection)



# Separating by species:
perch_survey_df <- filter(Data_survey_df, grepl('Perca fluviatilis', species))
trout_survey_df <- filter(Data_survey_df, grepl('Salmo trutta', species))
char_survey_df <- filter(Data_survey_df, grepl('Salvelinus alpinus', species))
pike_survey_df <- filter(Data_survey_df, grepl('Esox lucius', species))

perch_artsobs_df <- filter(Data_artsobs_df, grepl('Perca fluviatilis', species))
trout_artsobs_df <- filter(Data_artsobs_df, grepl('Salmo trutta', species))
char_artsobs_df <- filter(Data_artsobs_df, grepl('Salvelinus alpinus', species))
pike_artsobs_df <- filter(Data_artsobs_df, grepl('Esox lucius', species))


MakeSpDF <- function(df){
  Projection <- CRS("+proj=longlat +ellps=WGS84")
  sp_df <- SpatialPointsDataFrame(coords = df[,c("decimalLongitude","decimalLatitude")], 
                                  data = df[,c("occurrenceStatus","species")],
                                  proj4string = Projection)
  sp_df
}

perch_survey <- MakeSpDF(perch_survey_df)
trout_survey <- MakeSpDF(trout_survey_df)
char_survey <- MakeSpDF(char_survey_df)
pike_survey <- MakeSpDF(pike_survey_df)

perch_artsobs <- MakeSpDF(perch_artsobs_df)
trout_artsobs <- MakeSpDF(trout_artsobs_df)
char_artsobs <- MakeSpDF(char_artsobs_df)
pike_artsobs <- MakeSpDF(pike_artsobs_df)

# Now we have the covariates in 'Covariates', 
# as well as two types of data sets in 'Data_survey' and 'Data_artsobs', 
# and eight individual species dataframes.


# MESH --------------------------------------------------------------
Meshpars <- list(cutoff=0.08, max.edge=c(0.6, 3), offset=c(1,1))
#Meshpars <- list(cutoff=0.08, max.edge=c(1, 3), offset=c(1,1))

Mesh <- MakeSpatialRegion(data=NULL, bdry=norway.poly, meshpars=Meshpars,
                          proj = Projection)




