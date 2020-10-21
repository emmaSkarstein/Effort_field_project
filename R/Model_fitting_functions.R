# Model fitting functions


# PointedSDMs plus second spatial field

# Based on procedure from Emily Simmonds

# To fit integrated SDMs based on multiple data sets, with a second spatial field that describes variation 
# unique to the citizen science data.

# Uses Simpson approach for PP data
# Binomial model for PA data
# Using cloglog

library(ggplot2)
library(sf) # for sf objects
library(plyr)
library(dplyr) # for smoother dataframe-manipulation
library(ggmap) # also for nice maps
library(maps)
library(RColorBrewer)
library(PointedSDMs)
library(sp)
library(spatstat)
library(maptools)
library(INLA)
library(reshape2)
library(rgeos)
library(fields)
library(viridis)
library(spatialEco)







MakeStacks <- function(data_structured, data_unstructured, env_covariates, all_covariates, Mesh){
  
  norway <- ggplot2::map_data("world", region = "Norway(?!:Svalbard)")
  norway <- setdiff(norway, filter(norway, subregion == "Jan Mayen"))
  Projection <- CRS("+proj=longlat +ellps=WGS84")
  norwayfill <- map("world", "norway", fill=TRUE, plot=FALSE, 
                    ylim=c(58,72), xlim=c(4,32))
  IDs <- sapply(strsplit(norwayfill$names, ":"), function(x) x[1])
  norway.poly <- map2SpatialPolygons(norwayfill, IDs = IDs, 
                                     proj4string = Projection)
  
  # INTEGRATION STACK -------------------------------------------------
  stk.ip <- MakeIntegrationStack(mesh = Mesh$mesh, data = all_covariates, 
                                 area = Mesh$w, tag ='ip', InclCoords=TRUE)
  
  
  # SURVEY, STRUCTURED STACK -------------------------------------------
  # This depends on observation input. It also takes quite some time to match covariates to observation locations.
  NearestCovs_str <- GetNearestCovariate(points = data_structured, covs = env_covariates)
  NearestCovs_str@data[ , "int.survey"] <- 1 # add intercept 
  
  # Projector matrix from mesh to unstructured data
  projmat.str <- inla.spde.make.A(mesh = Mesh$mesh, loc = as.matrix(data_structured@coords))
  
  # If presences are Boolean, reformat
  if(is.logical(data_structured@data[,"occurrenceStatus"])) {
    data_structured@data[,"occurrenceStatus"] <- as.integer(data_structured@data[,"occurrenceStatus"])
    data_structured@data[,"Ntrials"] <- rep(1, nrow(data_structured@data))
  }
  
  stk.survey <- inla.stack(data=list(resp = cbind(NA, data_structured@data[,"occurrenceStatus"] ), 
                                     Ntrials = data_structured@data[,"Ntrials"]), 
                           A = list(1, projmat.str), 
                           tag = "survey",
                           effects = list(NearestCovs_str@data, 
                                          list(i = 1:Mesh$mesh$n,
                                               id.iid = 1:Mesh$mesh$n)))
  
  
  
  # ARTSOBS, UNSTRUCTURED STACK -----------------------------------------
  # Finding the covariates that are closest to the observation points
  NearestCovs_unstr <- GetNearestCovariate(points = data_unstructured, covs = all_covariates)
  NearestCovs_unstr@data[ , "int.artsobs"] <- 1 # add intercept 
  
  # Projector matrix from mesh to unstructured data
  projmat.artsobs <- inla.spde.make.A(mesh = Mesh$mesh, loc = as.matrix(data_unstructured@coords))
  
  stk.artsobs <- inla.stack(data = list(resp = cbind(rep(1,nrow(NearestCovs_unstr)), NA),
                                        e = rep(0, nrow(NearestCovs_unstr))), # why is this zero?
                            #e = rep(1, nrow(NearestCovs_unstr))),
                            #e = exp(rep(NearestCovs_unstr@data$log_area))),
                            A = list(1, projmat.artsobs), 
                            tag = "artsobs",
                            effects = list(NearestCovs_unstr@data, 
                                           list(shared_field = 1:Mesh$mesh$n, 
                                                bias_field = 1:Mesh$mesh$n, # This is for the second spatial field!
                                                id.iid = 1:Mesh$mesh$n))) 
  
  # PREDICTIONS -----------------------------------------------------------
  Nxy.scale <- 0.1 # use this to change the resolution of the predictions
  Nxy <- round(c(diff(norway.poly@bbox[1,]), diff(norway.poly@bbox[2,]))/Nxy.scale)
  stk.pred <- MakeProjectionGrid(nxy=Nxy, mesh=Mesh$mesh, data=all_covariates, 
                                 tag='pred', boundary=norway.poly)
  
  # Return list of all stacks
  return(list(ip = stk.ip, survey = stk.survey, artsobs = stk.artsobs, pred = stk.pred))
}

MakeStructuredStack <- function(data_structured, covariates, Mesh){
  
  # SURVEY, STRUCTURED STACK -------------------------------------------
  # This depends on observation input. It also takes quite some time to match covariates to observation locations.
  NearestCovs_str <- PointedSDMs::GetNearestCovariate(points = data_structured, covs = covariates)
  NearestCovs_str@data[ , "int.survey"] <- 1 # add intercept 
  
  # Projector matrix from mesh to unstructured data
  projmat.str <- INLA::inla.spde.make.A(mesh = Mesh$mesh, loc = as.matrix(data_structured@coords))
  
  # If presences are Boolean, reformat
  if(is.logical(data_structured@data[,"occurrenceStatus"])) {
    data_structured@data[,"occurrenceStatus"] <- as.integer(data_structured@data[,"occurrenceStatus"])
    data_structured@data[,"Ntrials"] <- rep(1, nrow(data_structured@data))
  }
  
  stk.survey <- INLA::inla.stack(data=list(resp = cbind(NA, data_structured@data[,"occurrenceStatus"] ), 
                                           Ntrials = data_structured@data[,"Ntrials"]), 
                                 A = list(1, projmat.str), 
                                 tag = "survey",
                                 effects = list(NearestCovs_str@data, 
                                                list(i = 1:Mesh$mesh$n,
                                                     id.iid = 1:Mesh$mesh$n)))
  stk.survey
}

MakeTestStack <- function(data_test, covariates, Mesh){
  
  data_test@data <- dplyr::mutate(data_test@data, NACol = rep(NA)) 
  data_test@data <- dplyr::mutate(data_test@data, zeroCol = rep(0))
  
  # TEST, STRUCTURED STACK -------------------------------------------
  NearestCovs_str <- PointedSDMs::GetNearestCovariate(points = data_test, covs = covariates)
  NearestCovs_str@data[ , "int.test"] <- 1 # add intercept 
  
  # Projector matrix from mesh to unstructured data
  projmat.str <- INLA::inla.spde.make.A(mesh = Mesh$mesh, loc = as.matrix(data_test@coords))
  
  # If presences are Boolean, reformat
  if(is.logical(data_test@data[,"occurrenceStatus"])) {
    data_test@data[,"occurrenceStatus"] <- as.integer(data_test@data[,"occurrenceStatus"])
    data_test@data[,"Ntrials"] <- rep(1, nrow(data_test@data))
  }
  
  stk.survey <- INLA::inla.stack(data=list(resp = cbind(NA, data_test@data[,"NACol"] ), 
                                           Ntrials = data_test@data[,"zeroCol"]), 
                                 A = list(1, projmat.str), 
                                 tag = "test",
                                 effects = list(NearestCovs_str@data, 
                                                list(shared_field = 1:Mesh$mesh$n,
                                                     id.iid = 1:Mesh$mesh$n)))
  stk.survey
}

# FITTING MODEL --------------------------------------------------------------------

MakeFormula <- function(cov_names, second_sp_field = FALSE, overdispersion = FALSE){
  intercepts <- "int.survey + int.artsobs - 1"
  env_effects <- paste(cov_names, collapse = ' + ')
  random_effects <- "f(shared_field, model = mesh.shared) + f(i, copy = 'shared_field')"
  
  if(second_sp_field){
    random_effects <- paste(random_effects, "+ f(bias_field, model = mesh.bias)")
  }
  if(overdispersion){
    random_effects <- paste(random_effects, "+ f(id.iid, model = 'iid')")
  }
  
  formula1 <- as.formula(paste(c("resp ~ 0 ", intercepts, env_effects, random_effects), collapse = " + "))
  
  if(is.null(cov_names)){
    formula1 <- as.formula(paste(c("resp ~ 0 ", intercepts, random_effects), collapse = " + "))
  }
  
  formula1
}

MakeFormulaSingleDataset <- function(cov_names, second_sp_field = FALSE, overdispersion = FALSE){
  intercepts <- "int.survey - 1"
  env_effects <- paste(cov_names, collapse = ' + ')
  random_effects <- "f(shared_field, model = mesh.shared) + f(i, copy = 'shared_field')"
  
  if(second_sp_field){
    random_effects <- paste(random_effects, "+ f(bias_field, model = mesh.bias)")
  }
  if(overdispersion){
    random_effects <- paste(random_effects, "+ f(id.iid, model = 'iid')")
  }
  
  formula1 <- as.formula(paste(c("resp ~ 0 ", intercepts, env_effects, random_effects), collapse = " + "))
  formula1
}



MakePred <- function(stk.pred, model){
  Projection <- CRS("+proj=longlat +ellps=WGS84")
  
  Pred <- SpatialPixelsDataFrame(points = stk.pred$predcoords, 
                                 data = model$predictions, 
                                 proj4string = Projection)
  Pred@data$precision <- Pred@data$stddev^-2
  
  return(Pred)
}

######## FitModel variations #########
FitModelCustom <- function(..., Formula, mesh, prior.range = c(10, 0.1), 
                           prior.sigma = c(0.1, 0.1)){
  stck <- inla.stack(...)
  
  #mesh <- inla.spde2.matern(Mesh$mesh)
  mesh.shared <- INLA::inla.spde2.pcmatern(mesh,
                                           prior.range = prior.range,
                                           #prior.sigma = c(0.5, 0.1)
                                           prior.sigma = prior.sigma
                                           )
  mesh.bias <- INLA::inla.spde2.pcmatern(mesh,
                                         prior.range = prior.range,
                                         prior.sigma = prior.sigma)
  
  # Fit model including predictions
  mod <- inla(Formula, family=c('poisson','binomial'),
              control.family = list(list(link = "log"), list(link = "cloglog")),
              data=inla.stack.data(stck), verbose=FALSE,
              control.results=list(return.marginals.random=FALSE,
                                   return.marginals.predictor=FALSE),
              control.predictor=list(A=inla.stack.A(stck), link=NULL, compute=TRUE),
              control.fixed = list(mean=0),
              Ntrials=inla.stack.data(stck)$Ntrials, E=inla.stack.data(stck)$e,
              control.compute=list(waic=FALSE, dic=FALSE))
  
  # For predictions
  id <- inla.stack.index(stck, "pred")$data
  pred <- data.frame(mean=mod$summary.fitted.values$mean[id],
                     stddev=mod$summary.fitted.values$sd[id])
  
  return(list(model = mod, predictions = pred))
}


CalcLinPred <- function(Model, resp){
  glm(resp ~ 1 + offset(Model$model$summary.linear.predictor[inla.stack.index(stk.test,"test")$data,"mean"]))
}


TrainTest = function(sb, data, k){
  Projection <- CRS("+proj=longlat +ellps=WGS84")
  training_fold_id <- list(sapply(1:k, function(s) sb$blocks@data$layer[sb$blocks@data$folds==s]))
  training_id <- list(sapply(1:k, function(s) which(sapply(slot(sb$blocks, 'polygons'), function(i) slot(i, 'ID')) %in% training_fold_id[[1]][[s]])))
  training_blocks <- list(sapply(1:k, function(s) SpatialPolygons(sb$blocks@polygons[training_id[[1]][[s]]], proj4string=Projection)))
  training_polydf <- list(sapply(1:k, function(s)
    SpatialPolygonsDataFrame(training_blocks[[1]][[s]], data.frame(row.names=sapply(slot(training_blocks[[1]][[s]], 'polygons'), function(i) slot(i, 'ID'))))))
  training_points <- list(sapply(1:k, function(s) point.in.poly(data, training_polydf[[1]][[s]]) %>% sp.na.omit()))
  return("train_test_folds" = training_points)
}












