calculateMukv = function(observations, predGrid, model, formulaString, ...) {

# Names changed for easier debugging 
    prG = predGrid
    obs = observations
    if (missing(formulaString) || is.null(formulaString))
        eq = dum~1 else eq = formulaString

    if (!"data" %in% getSlots(class(obs)) & 
          (terms(eq)[[3]] == 1 || 
            all(all.vars(eq)[-1] %in% dimnames(coordinates(obs))[[2]]))){ 
        obs = SpatialPointsDataFrame(obs,data = data.frame(dum = rep(1,dim(coordinates(obs))[1])))
        names(obs) = as.character(eq[[2]])
     }
#    prG$dum = 1
#    test = noquote(paste(dum2, predictors))
    red_ann_gam_krig = krige(eq, obs, prG, model) 
    return(mean(red_ann_gam_krig$var1.var))
#    obs = as.data.frame(obs[,1]) 
#    red_ann_gam_krig = as.data.frame(red_ann_gam_krig)
}
