get_predictions_dfsd <- function(model_outs, 
                            model, # model that was run to create the models in model_outs
                            env_data) {
  
  # choose the type and index for predict function
  if(model == 'lr'|model == 'gam'){
    
    type <- "response"
    index <- NULL
    
  } else if(model == 'rf'){
    
    type <- "prob"
    index <- 2
    
  }
  
  # number of bootstraps that were run
  k = length(model_outs)
  sdm = model_outs
  
  ## bootstrapped models
  print(paste0('#####   predicting from bootstrapped models   #####')) 
  
  ## predict from each of the bootstrapped models
  ## different workflow for lrReg and other methods
  if(model != 'lrReg') {
    
    # predict from each of the bootstrapped models and stack them together
    boots_out <- lapply(sdm$Bootstrapped_models, FUN = function(x) predict(x, newdata = env_data, type=type, index=index))
    
    ## quantiles
    print(paste0('#####   getting quantiles   #####'))
    mean_preds <- Reduce("+", boots_out)/length(boots_out) # the mean
    sd_preds <- rowVars(simplify2array(boots_out), na.rm=TRUE, std = TRUE) 
   
    
    
  } else if(model == 'lrReg') { 
    
    print(paste0('#####   predicting for lrReg bootstrapped models   #####'))
    
    ## convert variables to matrix
    covsMat <- as.matrix(rasterToPoints(env_data)) # convert variables to matrix
    
    ## predict from lrReg model
    boots <- stack(lapply(sdm$Bootstrapped_models, FUN = function(x) {
      
      pred <- predict(x, covsMat[, 3:ncol(covsMat)], type = "response") # predict from matrix
      
      pred <- as.matrix(cbind(covsMat[, 1:2], pred)) # combine predictions with east - norths
      
      if (any(is.na(pred[, 3]))) pred <- pred[-which(is.na(pred[,3])), ] # get rid of NAs
      
      pred_rast <- rasterize(pred[, 1:2], env_dat[[1]], field = pred[, 3]) # turn into a raster of probabilities
      
      ## return predictions
      return(pred_rast)
      
    }))
    
    ## check for intercept-only models which are 0.5 probability in all cells
    uniqueVals <-lapply(1:k,
                        function(x) { length(cellStats(boots[[x]], unique)) })
    
    drop <- which(uniqueVals <= 2) ## i.e. the mean and NA, because with intercept-only models the only values are 0.5 and NA
    
    # assign intercept only models an AUC of NULL - important for weighted average later
    if(any(drop)){
      
      sdm$AUC[drop] <- NA
      
      print(paste("Dropping", length(drop), "intercept-only model(s). Intercept-only models are given an AUC value of NA so they can be identified.
                  Where 1:(k-1) models are intercept only, only the non-intercept models are included in the final average. Where all models are intercept-only,
                  their predictions are returned but should not be used."))
      
    }
    
    
    if (length(drop) == k | length(drop) == 0) {
      
      ## quantiles
      print(paste0('#####   getting quantiles lrReg   #####'))
      mean_preds <- mean(boots) # where all models are intercept-only, takes the mean to avoid errors later but AUC scores are NA which means they are dropped for final ensembles
      sd_preds <- rowVars(simplify2array(boots_out), na.rm=TRUE, std = TRUE)
     
      
      
    } else {
      
      print(paste0('#####   getting quantiles lrReg   #####'))
      mean_preds <- mean(boots[[-drop]])
      sd_preds <- rowVars(simplify2array(boots_out), na.rm=TRUE, std = TRUE)
      
    }
    
    ## if some models were intercept-only then recalculate k (number of models used)
    k <- k - length(drop)
    
  }
  
  return(list(mean_predictions = mean_preds,
              sd_predictions = sd_preds))
  
  }