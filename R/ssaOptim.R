######################################################
# SCRIPT TO RUN SSA USING 1 CRITERION
#######################################################

ssaOptim = function(observations, predGrid, candidates, action, nDiff, model, 
                    nr_iterations, plotOptim, formulaString = NULL, ...){

######################Set up initial network(thinned/augmented)#############################
  nn = dim(coordinates(observations))[1]
  nGrd = dim(coordinates(predGrid))[1]
  netPtsInit = observations
  lmcfit = NULL
  if (action == "del"){
    sn=sample(nn)
    delPts = observations[which((sn)<(min(sn)+nDiff)),]
    netPts = observations[which((sn)>=(min(sn)+nDiff)),]
    addPts = NULL
#    netPts = rbind(netPts,addPts)
    nn = dim(coordinates(netPts))[1]
    crit1 = calculateMukv(observations = netPts, predGrid = predGrid, model = model, 
                          formulaString = formulaString, ...) # initial design

    res=ssaMap(candidates, predGrid, model,   
                max_points_shift = 1, # maximum number of points to move in SSA (should be 1)
                maxShiftFactorX = 0.2, # maximum shift for points to move in SSA as percentage of total study area in X direction
                minShiftFactorX = 0, # minimum shift for points to move in SSA as percentage of total study area in X direction
                maxShiftFactorY = 0.2, # maximum shift for points to move in SSA as percentage of total study area in Y direction
                minShiftFactorY = 0, # maximum shift for points to move in SSA as percentage of total study area in Y direction
                start_p = 0.2,
                countMax = 200, netPts, addPts, delPts, crit1, nn, action, nDiff, 
                netPtsInit, nr_iterations, plotOptim, formulaString = formulaString, ...)
    if ("data" %in% getSlots(class(netPts)) && !"data" %in% getSlots(class(res))) 
      res = SpatialPointsDataFrame(res,data = netPts@data[as.numeric(row.names(res@data)),])
    return(res)  
    }
  if (action == "add"){
# Checking if the formulaString contains auxiliary variables, except from "coordinates"
    cnames = dimnames(coordinates(observations))[[2]]
    if (!is.null(formulaString) && terms(formulaString)[[3]] != 1 &&
         !all(all.vars(formulaString)[-1] %in% cnames)) {
      observations = SpatialPointsDataFrame(coordinates(observations),
                     data = model.frame(terms(formulaString),observations),coords.nrs = c(1,2))
      if (length(names(predGrid)) ==0) predGrid = SpatialPointsDataFrame(coordinates(predGrid),
                            data = data.frame(dum = 1:nGrd))
      predGrid = SpatialPointsDataFrame(coordinates(predGrid), 
                     data = model.frame(delete.response(terms(formulaString)),predGrid),coords.nrs = c(1,2))
# Adding "observations" to predGrid, for later estimate of linear model in calculateMukv
      predGrid$dum = krige(formulaString,observations,predGrid,model = model)$var1.pred
      lp = length(names(predGrid))
      names(predGrid)[lp] = names(observations)[1]
# Deleting coordinates from 
      dno = which(names(observations) %in% dimnames(coordinates(predGrid))[[2]])
      if (length(dno) > 0) observations = observations[,-dno] 
      dnp = which(names(predGrid) %in% dimnames(coordinates(predGrid))[[2]])
      if (length(dnp) > 0) predGrid = predGrid[,-dnp] 

      models = list()
      for (i in 1:length(names(observations))) {
        if (i == 1) {
          models[[i]] = (autofitVariogram(formulaString,observations))$var_model
        } else models[[i]] = (autofitVariogram(as.formula(paste(names(observations)[i],"~1")),predGrid))$var_model
      }

    } else {
      observations = SpatialPoints(observations)
      predGrid = SpatialPoints(predGrid)
    }
    sn=sample(nGrd)
    delPts = NULL
    addPts = predGrid[which(sn<=nDiff),]
    netPts = rbind(observations,addPts)
    nn = dim(coordinates(netPts))[1]
    crit1 = calculateMukv(observations = netPts, predGrid = predGrid, model = model, 
                           formulaString = formulaString, ...) # initial design

    res=ssaMap(candidates, predGrid, model,  
                max_points_shift = 1, # maximum number of points to move in SSA (should be 1)
                maxShiftFactorX = 0.2, # maximum shift for points to move in SSA as percentage of total study area in X direction
                minShiftFactorX = 0, # minimum shift for points to move in SSA as percentage of total study area in X direction
                maxShiftFactorY = 0.2, # maximum shift for points to move in SSA as percentage of total study area in Y direction
                minShiftFactorY = 0, # maximum shift for points to move in SSA as percentage of total study area in Y direction
                start_p = 0.2,
                countMax = 200, netPts, addPts, delPts, crit1, nn, 
                action, nDiff, netPtsInit, nr_iterations, plotOptim, formulaString = formulaString, 
                models = models, ...)  
    if ("data" %in% getSlots(class(netPts)) && !"data" %in% getSlots(class(res))) 
      res = SpatialPointsDataFrame(res,data = netPts@data[as.numeric(row.names(res@data)),])
    return(res)
  }
}



