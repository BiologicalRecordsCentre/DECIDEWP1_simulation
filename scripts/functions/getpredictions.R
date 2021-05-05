get_predictions <- function(model_outs, 
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
    boots_out <- raster::stack(lapply(sdm$Bootstrapped_models, FUN = function(x) predict(env_data, x, type=type, index=index)))
    
    ## quantiles
    print(paste0('#####   getting quantiles   #####'))
    mean_preds <- calc(boots_out, fun = mean, na.rm = T) # the mean
    quant_preds <- calc(boots_out, fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)}) # get the quantile variation
    rnge <- quant_preds[[2]]-quant_preds[[1]] # get the range of max - min
    
    
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
      quant_preds <- calc(boots, fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)})
      rnge <- quant_preds[[2]]-quant_preds[[1]]
      
      
    } else {
      
      print(paste0('#####   getting quantiles lrReg   #####'))
      mean_preds <- mean(boots[[-drop]])
      quant_preds <- calc(boots[[-drop]], fun = function(x) {quantile(x, probs = c(0.05, 0.95), na.rm = TRUE)})
      rnge <- quant_preds[[2]]-quant_preds[[1]]
      
    }
    
    ## if some models were intercept-only then recalculate k (number of models used)
    k <- k - length(drop)
    
  }
  
  return(list(mean_predictions = mean_preds,
              quant_minmax = quant_preds, 
              quant_range = rnge))
  
  }