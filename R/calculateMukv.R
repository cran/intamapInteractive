calculateMukv = function(observations, predGrid, model, formulaString, ...) {

    prG = predGrid
    obs = observations
    if (missing(formulaString)) {
      numPredictors = dim(prG)[2] # JS: Not necessary with separate xy coordinates if predGrid is spdf - 2 # -2 for two xy coordinates 
      if (!is.null(numPredictors) && numPredictors > 0)  {
        predictors = names(prG[1:length(names(prG))])
      } else  {
        predictors = "1"
      }
      if ("data" %in% getSlots(class(obs))) {
        obs$dum = 1
      } else obs = SpatialPointsDataFrame(obs,data = data.frame(dum = rep(1,dim(coordinates(obs))[1])))
      predictors = paste(predictors, collapse =" + ")
      dum2 = "dum ~"
      eq = as.formula(paste(noquote(paste(dum2, predictors))))
    } else eq = formulaString
    
#    prG$dum = 1
#    test = noquote(paste(dum2, predictors))
    red_ann_gam_krig = krige(eq, obs, prG, model) 
    return(mean(red_ann_gam_krig$var1.var))
#    obs = as.data.frame(obs[,1]) 
#    red_ann_gam_krig = as.data.frame(red_ann_gam_krig)
}
