# Full model fitting and validation


# Model validation

library(blockCV)
library(foreach)
library(doParallel)

setwd("Citizen_Science_Skarstein_master") # This is for running on external servers, can be removed when running locally

# Loading observations and covariates
source("R/01_loading_map_obs_covs.R")

# Loading functions to fit models
source("R/Model_fitting_functions.R")

# Loading functions to visualize models
source("R/Model_visualization_functions.R")



# DIVIDE INTO TRAINING AND TEST FOLDS FOR MODEL VALIDATION
k <- 5
sb <- spatialBlock(speciesData = trout_survey, # sf or SpatialPoints
                   species = "occurrenceStatus", # the response column (binomial or multi-class)
                   #rasterLayer = validation_raster, # a raster for background (optional)
                   theRange = 150000, # size of the blocks in meters
                   k = k, # number of folds
                   selection = "random",
                   iteration = 100, # find evenly dispersed folds
                   biomod2Format = FALSE)

# SELECTING TRAINING AND TEST SETS
survey_traintest <- TrainTest(sb, trout_survey, sb$k)
artsobs_traintest <- TrainTest(sb, trout_artsobs, sb$k)


# Setting up parallel backend
cl <- parallel::makeForkCluster(4)
parallel <- TRUE
if(parallel){
  doParallel::registerDoParallel(cl)
}else{
  foreach::registerDoSEQ()
}


modelList <- foreach::foreach(i = 1:k) %dopar% {
  
  survey_train <- do.call("rbind", unlist(sapply((1:k)[-i], function(s) survey_traintest[[1]][[s]])))
  artsobs_train <- do.call("rbind", unlist(sapply((1:k)[-i], function(s) artsobs_traintest[[1]][[s]])))
  survey_test <- survey_traintest[[1]][[i]]
  # artsobs_test: we don't use the CS data for testing.
  
  # MAKE STACKS
  stks <- MakeStacks(data_structured = survey_train, data_unstructured = artsobs_train,
                     env_covariates = env_covariates, all_covariates = Covariates, Mesh = Mesh)

  stk.survey <- stks$survey
  stk.artsobs <- stks$artsobs
  stk.ip <- stks$ip
  stk.pred <- stks$pred
  
  stk.test <- MakeTestStack(survey_test, env_covariates, Mesh)
  
  # CONSTRUCT FORMULAS 
  Use <- c("decimalLongitude","decimalLatitude", "log_area", 
           "log_catchment", "eurolst_bio10", "SCI")
  Use_CS <- c(Use, "distance_to_road", "HFP")
  
  overdisp <- FALSE
  
  formula1 <- MakeFormula(cov_names = Use, second_sp_field = FALSE, overdispersion = overdisp)
  formula2 <- MakeFormula(cov_names = Use, second_sp_field = TRUE, overdispersion = overdisp)
  formula3 <- MakeFormula(cov_names = Use_CS, second_sp_field = FALSE, overdispersion = overdisp)
  formula4 <- MakeFormula(cov_names = Use_CS, second_sp_field = TRUE, overdispersion = overdisp)
  formula0 <- MakeFormulaSingleDataset(cov_names = Use, second_sp_field = FALSE, overdispersion = overdisp)
  formula5 <- MakeFormula(cov_names = NULL, second_sp_field = TRUE, overdispersion = overdisp)
  
  # FIT MODELS
  model1 <- FitModelCustom(stk.survey, stk.artsobs, stk.ip, stk.pred$stk, stk.test,
                         Formula = formula1, mesh = Mesh$mesh)

  model2 <- FitModelCustom(stk.survey, stk.artsobs, stk.ip, stk.pred$stk, stk.test,
                         Formula = formula2, mesh = Mesh$mesh)

  model3 <- FitModelCustom(stk.survey, stk.artsobs, stk.ip, stk.pred$stk, stk.test,
                         Formula = formula3, mesh = Mesh$mesh)

  model4 <- FitModelCustom(stk.survey, stk.artsobs, stk.ip, stk.pred$stk, stk.test,
                         Formula = formula4, mesh = Mesh$mesh)
  
  model0 <- FitModelCustom(stk.survey, stk.ip, stk.pred$stk, stk.test,
                         Formula = formula0, mesh = Mesh$mesh)
  
  model5 <- FitModelCustom(stk.survey, stk.artsobs, stk.ip, stk.pred$stk, stk.test,
                         Formula = formula5, mesh = Mesh$mesh)
  
  # CALCULATE DIC
  resp <- survey_test$occurrenceStatus
  mod_res1 <- CalcLinPred(model1, resp)
  mod_res2 <- CalcLinPred(model2, resp)
  mod_res3 <- CalcLinPred(model3, resp)
  mod_res4 <- CalcLinPred(model4, resp)
  mod_res0 <- CalcLinPred(model0, resp)
  mod_res5 <- CalcLinPred(model5, resp)
  
  
  list(model0 = mod_res0$deviance,    
       model1 = mod_res1$deviance, 
       model2 = mod_res2$deviance, 
       model3 = mod_res3$deviance, 
       model4 = mod_res4$deviance,
       model5 = mod_res5$deviance)
}

parallel::stopCluster(cl)

saveRDS(modelList, "R/output/cv_output_4mods.RDS")

res <- readRDS("R/output/cv_output_4mods.RDS")

res_mat <- data.frame(do.call("cbind", res))
res_mat[] <- lapply(res_mat, as.numeric)
rowMeans(res_mat)


