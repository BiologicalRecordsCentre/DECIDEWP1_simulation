
###' Fit SMDs - main function for fitting SDMs to data from create presence absence
###'

fsdm <- function(species, model, climDat, spData, k, write, outPath, #inters = F, prediction = TRUE,
                 knots_gam = -1){ 
  
  # select species
  ind <- which(names(spData) == species)
  spDat <- spData[[ind]]
  
  if (is.null(spDat)) {
    out <- NULL
  }
  else {
    
    #### extract data ####
    pres <- data.frame(val = 1, raster::extract(x = climDat, y = spDat$Presence, na.rm = T), spDat$Presence)
    
    if (any(is.na(pres))) { # this might not be needed with the new screenraster argument for the create pseudoabs function
      
      warning("!!!   NAs in presences   !!!")
      
      # pres <- na.omit(pres)
    }
    
    nRec <- nrow(pres)
    
    print(paste("Occurrence records:", nRec))
    
    ab <- data.frame(val = 0, raster::extract(x = climDat, y = spDat$pseudoAbsence, na.rm = T), spDat$pseudoAbsence)
    
    if (any(is.na(ab))) {
      
      warning("!!!   NAs in pseudoabsences   !!!")
      
      # ab <- na.omit(ab)
    }
    
    # get a data frame with the lat-lon coordinates
    allDat <- rbind(pres[!names(pres) %in% c("lon", "lat")], ab[!names(ab) %in% c("lon", "lat")])
    allDat_loc <- rbind(pres, ab)
    
    
    # #### remove any columns with only 0 values - useless for the models
    # # I suppose I could change the threshold for 'importance', 
    # # i.e. remove variables that have very few unique values although 
    # # would have to be careful with not removing the presence absence column
    # # !!  Not sure this works with the predict function (only tested for lrReg)
    # if(any(colSums(abs(allDat)) == 0)){
    #   
    #   print(paste("!!  Removing variables that have only 0 values. variables removed = ", 
    #               names(allDat)[colSums(abs(allDat)) == 0], " !!"))
    #   
    #   allDat <- allDat[, colSums(abs(allDat)) != 0] # only keep columns that have some values in 
    #   allDat_loc <- allDat_loc[, colSums(abs(allDat_loc)) != 0] # only keep columns that have some values in 
    #   
    # }
    
    
    # ## raster layer needed as matrix for lrReg model
    # # don't need at the moment because I am predicting outside of the fitSDM function
    # if (model == "lrReg") {
    #   
    #   covsMat <- as.matrix(rasterToPoints(climDat))
    #   
    # }
    
    ## determine the weights argument for all models
    if ((model != "lrReg" & model != 'lr' & model != "gam" & model != "rf") & nRec != nrow(ab)) { 
      
      warning("Prevalence is not 0.5 and no weights are applied to account for this. Currently weights are only applied where model = lrReg, lr, gam or rf")
      
      if(model == 'me'){
        
        options(java.parameters = "-Xmx3g")
        
        warning("Maxent ('me') models work best with equal number of presences and absences, matching the number of presences and absences")
        
        if(nrow(pres) > nrow(ab)){
          
          warning("Number of presences > number of absences, reducing number of presences to match absences")
          
          pres <- pres[sample(x = 1:nrow(pres), size = nrow(ab)),]
          
          # match the number of absences to presences in the spDat 
          spDat$Presence <- spDat$Presence[sample(x = 1:nrow(spDat$Presence),
                                                  size = nrow(spDat$pseudoAbsence)),]
          
        } else if(nrow(pres) < nrow(ab)){
          
          warning("Number of abesences > number of presences, reducing number of presences to match absences")
          
          ab <- ab[sample(x = 1:nrow(ab), size = nrow(pres)),]
          
          # match the number of absences to presences in the spDat 
          spDat$pseudoAbsence <- spDat$pseudoAbsence[sample(x = 1:nrow(spDat$pseudoAbsence),
                                                            size = nrow(spDat$Presence)),]
          
        }
        
        # get a data frame with the lat-lon coordinates
        allDat <- rbind(pres[!names(pres) %in% c("lon", "lat")], ab[!names(ab) %in% c("lon", "lat")])
        allDat_loc <- rbind(pres, ab)
        
        
        # match the number of absences to presences in the spDat 
        spDat$pseudoAbsence <- spDat$pseudoAbsence[sample(x = 1:nrow(spDat$pseudoAbsence),
                                                          size = nrow(spDat$Presence)),]
        
        
        
      } else if(model != "lrReg" | model != 'lr' | model != "gam" | model != "rf"){
        
        stop("Model is not one of lr, lrReg, me, gam or rf. Stopping model run")
        
      }
      
      
    }
    
    if ((model == "lrReg"|model == 'lr'|model == "gam"|model == "rf") & nRec != nrow(ab)) {
      
      print("Prevalence is not 0.5. Weighting absences to simulate a prevalence of 0.5")
      
      nAb <- nrow(ab)
      
      prop <- nRec / nAb
      
      print(paste("Absence weighting:", prop))
      
    } else if(model == 'me') {
      
      prop <- NULL 
      
    } else { 
      
      stop('Model specified not accepted by fitSDM (problem with prop calculation)')
      
    }
    
    
    ###### Move straight to 'bootstrapping' ######
    
    print('######     Bootstrapping    ######')
    
    folds <- c(kfold(pres, k), kfold(ab, k))
    folds_me_pres <- kfold(spDat$Presence, k)
    folds_me_ab <- kfold(spDat$pseudoAbsence, k)
    
    e <- vector('list', length = k)
    mods_out <- vector('list', length = k)
    
    for (i in 1:k) {
      
      
      if (model == 'me') {
        
        train_me_pres <- spDat$Presence[folds_me_pres != i, ]
        train_me_abs <- spDat$pseudoAbsence[folds_me_ab != i, ]
        
        test_me_pres <- spDat$Presence[folds_me_ab == i, ]
        test_me_abs <- spDat$pseudoAbsence[folds_me_ab == i, ]
        
        mod <- maxent(x = climDat, p = data.frame(train_me_pres)[,1:2], 
                      a = data.frame(train_me_abs)[,1:2])
        
        e[[i]] <- dismo::evaluate(p = test_me_pres, 
                                  a = test_me_abs, 
                                  x = climDat,
                                  mod, tr = seq(0, 1, length.out = 200))
        
      }
      
      if(model == 'lr' | model == 'rf' | model == 'gam'){
        
        train <- allDat[folds != i, ]
        
        ## set the weights argument for models
        if ((model == "lrReg"|model == 'lr'|
             model == "gam"|model == "rf") & nRec != nrow(ab)){
          
          weights <- c(rep(1, length(train$val[train$val == 1])), rep(prop, length(train$val[train$val == 0])))
          
        } else if(nRec == nrow(ab)){ weights <- NULL } 
        
        test <- allDat[folds == i, ]
        
        if (model == "lr") {
          mod <- glm(val ~ ., data = train, 
                     family = binomial(link = "logit"),
                     weights = weights)
          
          ## weight the minority class by the ratio of majority/minority.
          ## minority is always the smaller class in my datasets, so give the 
          ## presence values a larger weight that the minorities
          ## round the weights argument to make sure those that are close to equilibrium
          ## aren't affected...?
        }
        else if (model == "rf") {
          
          # if records and absences are matched, run without weights (weights = NULL doesn't work)
          # if they aren't then implement weights argument
          if(nRec == nrow(ab)){
            
            mod <- randomForest(x = train[, 2:ncol(train)], 
                                y = as.factor(train[, 1]), 
                                importance = T, 
                                norm.votes = TRUE)
            
          } else if(nRec != nrow(ab)){
            
            mod <- randomForest(x = train[, 2:ncol(train)], 
                                y = as.factor(train[, 1]), 
                                importance = T, 
                                norm.votes = TRUE,
                                classwt = list(unique(weights)[1],
                                               unique(weights)[2])) ## must be list(presences, absences)
          }
          
        }
        else if (model == "gam"){
          
          ## create formula for gam
          l <- sapply(allDat, unique)
          ks <- rownames_to_column(data.frame(k = round(sapply(l, length))[-1]),
                                   var = "variable")
          
          
          # drop variables according to number of knots asked for
          # -1 is basically 9 knots
          if(knots_gam == -1) {
            v_keep <- ks[ks$k > 11,]
            
            print(paste("variable dropped =", ks$variable[!ks$variable %in% v_keep$variable]))
          }
          
          # any others just keep the variables with over the number of knots
          if(knots_gam > 0) {
            v_keep <- ks[ks$k > (knots_gam+3),]
            print(paste("variable dropped =", ks$variable[!ks$variable %in% v_keep$variable]))
          }
          # v_keep
          
          form <- as.formula(paste0("val ~ s(", paste(v_keep$variable,
                                                      ", k = ", knots_gam) %>% #,
                                      # v_keep$knts) %>%
                                      paste0(collapse = ") + s("), ")"))
          # form
          
          
          mod <- gam(formula = form, data = train, 
                     family = binomial(link = 'logit'), 
                     select = TRUE, method = 'REML', gamma = 1.4,
                     weights = weights)
        }
        
        e[[i]] <- dismo::evaluate(p = test[test$val == 1, ], 
                                  a = test[test$val == 0, ], 
                                  mod, tr = seq(0, 1, length.out = 200))
        
        
      }
      
      ## now implement lasso regression
      if(model == "lrReg"){
        
        train <- allDat[folds != i, ]
        
        ## set the weights argument for models
        if (model == "lrReg" & nRec != nrow(ab)) weights <- c(rep(1, length(train$val[train$val == 1])), rep(prop, length(train$val[train$val == 0])))
        
        test <- allDat[folds == i, ]
        
        ## test weights for the testing data 
        if (model == "lrReg" & nRec != nrow(ab)) testweights <- c(rep(1, length(test$val[test$val == 1])), rep(prop, length(test$val[test$val == 0])))
        
        
        mod <- glmnet::cv.glmnet(x = as.matrix(train[, 2:ncol(train)]),
                                 y = train[, 1],
                                 family = "binomial",
                                 nfolds = 3,
                                 weights = weights)
        
        # evaluate model on the testing data
        eval <- assess.glmnet(mod, newx = as.matrix(test[,2:ncol(test)]), newy = test[,1], weights = testweights)
        roc <- roc.glmnet(mod, newx = as.matrix(test[,2:ncol(test)]), newy = test[,1], weights = testweights)
        
        e[[i]] <- list(eval, roc)
        
        names(e[[i]]) <- c('evaluation', 'roc')
        
      }
      
      # store all the bootstrapped models
      mods_out[[i]] <- mod
      
    }
    
    ## get the auc from each run of the bootstrapping
    if(model != 'lrReg'){
      
      auc_val <- sapply(e, function(x) {
        slot(x, "auc")
      })
      
    } else if(model == 'lrReg'){
      
      auc_val <- sapply(c(1:k), function(x) e[[x]]$evaluation$auc)
      
    }
    
    
    ## get the mean AUC across all models
    meanAUC <- mean(auc_val)
    
    
    ## old storage when predicting from full model
    # out <- NULL
    # out <- list(species, nRec, fullMod, auc, k, pred, allDat_loc, e, mods_out)
    # names(out) <- c("Species", "Number of records", "Model", 
    #                 "AUC", "Number of folds for validation", "Predictions", 
    #                 "Data", "Model_evaluation", "Bootstrapped_models")
    
    out <- NULL
    out <- list(species, nRec,
                auc_val, meanAUC, k, 
                allDat_loc, ## taken out to see if it reduces file sizes
                e, mods_out)
    names(out) <- c("Species", "Number of records", 
                    "AUC", "meanAUC", "Number of folds for validation",
                    "Data",
                    "Model_evaluation", "Bootstrapped_models")
    
    if (write == TRUE) {
      print(species)
      save(out, file = paste0(outPath, species, "_", model, 
                              ".rdata"))
    }
  }
  return(out)
}



#### Presence absence function
cpa <- function (spdat, species, minYear, maxYear, nAbs, matchPres = FALSE, 
                 recThresh, replace = F, screenRaster = NULL) 
{
  dat <- spdat[spdat$year >= minYear & spdat$year <= maxYear,]
  pres <- dat[dat$species == species, c("lon", "lat")]
  
  if (nrow(pres) < recThresh) {
    warning("Number of records does not exceed recThresh")
    out <- NULL
  }
  
  else {
    
    ab <- dat[dat$species != species, c("lon", "lat")]
    
    ab <- ab[!ab %in% pres] # remove all the abesences from the same places as presences
    
    if (nrow(ab) < nrow(pres)) {
      warning(paste("More presences than possible locations for absences. Consider lowering the number of pseudo absences. Reducing the number of presences for:", species, 
                    "to:", nrow(ab)))
      pres <- pres[sample(x = 1:nrow(pres), size = nrow(ab)),]
    }
    
    if (matchPres == TRUE) {
      sampInd <- sample(1:nrow(ab), nrow(pres))
    } else {
      
      if (nAbs <= nrow(pres)){ # if presences > nAbs argument, then increase nAbs to match number of presences  
        sampInd <- sample(1:nrow(ab), nrow(pres))
      } else if (nAbs <= nrow(ab) & nAbs >= nrow(pres)) {
        sampInd <- sample(1:nrow(ab), nAbs)
      } else {
        warning(paste0("Fewer than", nAbs, "locations available for pseudo absences when using the target group approach. Setting nAbs to the maximum number possible (", nrow(ab), ")."))
        sampInd <- 1:nrow(ab)
      }
      
    }
    
    ab <- ab[sampInd, ]
    out <- list(pres, ab)
    names(out) <- c("Presence", "pseudoAbsence")
    
    ## if screenRaster is specified, check if any presence or absence points fall outside of the raster extent (i.e. they are NA).
    ## If some data fall outside of the extent of the covariates, drop them, and drop the equivalent number of absences or presences
    ## to ensure they are equal in number.
    ## if matchPres == FALSE and screenRaster = TRUE then some species might have number Presence > number pseudoabsence
    ## because absences might happen to fall in NA values
    ## This is addressed in the fsdm() function
    
    if (!is.null(screenRaster)) {
      
      for (i in 1:dim(screenRaster)[3]) {
        
        presDrop <- raster::extract(screenRaster[[i]], out$Presence)
        
        abDrop <- raster::extract(screenRaster[[i]], out$pseudoAbsence)
        
        if (any(is.na(presDrop))) out$Presence <- out$Presence[-which(is.na(presDrop)), ]
        
        if (any(is.na(abDrop))) out$pseudoAbsence <- out$pseudoAbsence[-which(is.na(abDrop)), ]
        
      }
      
      if (matchPres == TRUE) {
        
        if (nrow(out$Presence) > nrow(out$pseudoAbsence)) {
          
          out$Presence <- out$Presence[1:nrow(out$pseudoAbsence), ]
          
        } else if (nrow(out$Presence) < nrow(out$pseudoAbsence)) {
          
          out$pseudoAbsence <- out$pseudoAbsence[1:nrow(out$Presence), ]
          
        }
        
      }
      
    }
  }
  
  return(out)
  
}


##### get predictions from the bootstrapped model outputs
# takes the model object from fsdm()
# the model used to run it ('lr', 'rf' etc...)
# and the environmental data to predict onto
#' 
#' Takes the bootstrapped models from the fsdm() function above,
#' the model used to run the function above
#' and the environmental data you want to predict onto (with the 
#' same layers as used to run the models above)
#' and outputs a three item list:
#' the mean predictions across all of the bootstrapped modes
#' the min and max quantiles (0.05 and 0.95)
#' and the range between them
#' 
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
  k = length(model_outs$Bootstrapped_models)
  
  ## bootstrapped models
  print(paste0('#####   predicting from bootstrapped models   #####')) 
  
  ## predict from each of the bootstrapped models
  ## different workflow for lrReg and other methods
  if(model != 'lrReg') {
    
    # predict from each of the bootstrapped models and stack them together
    boots_out <- raster::stack(lapply(model_outs$Bootstrapped_models, FUN = function(x) predict(env_data, x, type=type, index=index)))
    
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
    boots <- stack(lapply(model_outs$Bootstrapped_models, FUN = function(x) {
      
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
      
      model_outs$AUC[drop] <- NA
      
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



#### Convert grid reference from BRC data to eastings and northings
c_en <- function(occ){
  
  coords <- gr_let2num(gridref = occ$TO_GRIDREF, centre = T)
  
  occ$cn <- gr_det_country(occ$TO_GRIDREF)
  
  occ <- cbind(occ, coords)
  
  if (length(unique(occ$cn) > 1)) {
    
    GBCRS <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")
    
    NICRS <- CRS("+proj=tmerc +lat_0=53.5 +lon_0=-8 +k=1 +x_0=200000 +y_0=250000 +ellps=airy +towgs84=482.5,-130.6,564.6,-1.042,-0.214,-0.631,8.15 +units=m +no_defs")
    
    NIocc <- occ[occ$cn == "OSNI",]
    
    GBocc <- occ[occ$cn == "OSGB",]
    
    NIcoords <- NIocc[, c("EASTING", "NORTHING")]
    
    coordinates(NIcoords) <- c("EASTING", "NORTHING")
    
    proj4string(NIcoords) <- NICRS
    
    NIcoords <- spTransform(NIcoords, GBCRS)
    
    NIocc[,c("EASTING", "NORTHING")] <- data.frame(NIcoords)
    
    occ <- rbind(GBocc, NIocc)
    
  }
  occ
}

### Thin records from Richard Hassall
thin_records = function(occs, templateGrid){
  require(raster)
  require(sp)
  
  # set where to store temporary raster files...
  oldTempDir = rasterOptions(setfileext=T)$tmpdir # include a default option to stop it printing information
  myTempDir = paste(oldTempDir, "thin_records", sep="/")
  assign("rasterOptions", rasterOptions(tmpdir=myTempDir), envir = .GlobalEnv)
  
  # convert records to presence on the template grid...
  rownames(occs) = 1:nrow(occs)
  R = rasterize(occs[,1:2], templateGrid, field=1, max)
  xy = coordinates(R)[!is.na(getValues(R)),] # presence grid cell coords
  
  # reset the temporary raster files and delete any created in the function
  assign("rasterOptions", rasterOptions(tmpdir=oldTempDir), envir = .GlobalEnv)
  unlink(myTempDir, recursive=TRUE)
  
  # return SpatialPoints of the presence grid cell coords...
  SpatialPoints(xy, proj4string=CRS(projection(templateGrid)))
  
}


### Get AUCs for different models
getAUC <- function(file) { # list of different SDMs
  
  f <- list()
  for(i in 1:length(file)){
    df <- file[[i]]
    f[[i]] <- data.frame(species = df$Species, 
                         auc = df$AUC, model = class(df$Model)[1])
  }
  
  f <- do.call("rbind", f)
  return(f)
}


### Average models with AUC weighting
sdm_av <- function (sdm_out, skillDat, species){
  
  rfSkill <- skillDat$auc[skillDat$species == species & skillDat$model == 
                            "randomForest"]
  lrSkill <- skillDat$auc[skillDat$species == species & skillDat$model == 
                            "glm"]
  meSkill <- skillDat$auc[skillDat$species == species & skillDat$model == 
                            "MaxEnt"]
  
  if (rfSkill >= 0.7 | lrSkill >= 0.7| meSkill >= 0.7) {
    
    lr_m <- sdm_out[[1]]
    rf_m <- sdm_out[[2]]
    me_m <- sdm_out[[3]]
    
    if(class(lr_m$Model)[1] != "glm" | 
       class(rf_m$Model)[1] != "randomForest" | 
       class(me_m$Model)[1] != "MaxEnt") {
      
      stop("list not in right order. Order must be -> glm, ranfor, maxent")
      
    }
    
    lrRast <- lr_m$Predictions
    rfRast <- rf_m$Predictions
    meRast <- me_m$Predictions
    
    stack <- stack(rfRast, lrRast, meRast)
    print(c(rfSkill, lrSkill, meSkill))
    if (rfSkill < 0.7) {
      rfSkill <- 0
    }
    if (lrSkill < 0.7) {
      lrSkill <- 0
    }
    if (meSkill < 0.7) {
      meSkill <- 0
    }
    print(c(rfSkill, lrSkill, meSkill))
    ensemble <- weighted.mean(x = stack, w = c(rfSkill, lrSkill, meSkill))
    ensemble
  }
}


# function to calculate partial dependence for a single predictor from eBird tutorial
calculate_pd <- function(predictor, model, data, 
                         x_res = 25, n = 1000) {
  # create prediction grid
  rng <- range(data[[predictor]], na.rm = TRUE)
  x_grid <- seq(rng[1], rng[2], length.out = x_res)
  grid <- data.frame(covariate = predictor, x = x_grid, 
                     stringsAsFactors = FALSE)
  names(grid) <- c("covariate", predictor)
  
  # subsample training data
  n <- min(n, nrow(data))
  s <- sample(seq.int(nrow(data)), size = n, replace = FALSE)
  data <- data[s, ]
  
  # drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)
  
  # predict
  p <- predict(model, data = grid)
  
  # summarize
  pd <- grid[, c("covariate", predictor)]
  names(pd) <- c("covariate", "x")
  pd$pred <- p$predictions[, 2]
  pd <- dplyr::group_by(pd, covariate, x) %>% 
    dplyr::summarise(pred = mean(pred, na.rm = TRUE)) %>% 
    dplyr::ungroup()
  
  return(pd)
}
