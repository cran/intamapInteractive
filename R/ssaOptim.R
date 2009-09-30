######################################################
# SCRIPT TO RUN SSA USING 1 CRITERION
#######################################################

ssaOptim = function(observations, predGrid, candidates, action, nDiff, model, nr_iterations, plotOptim, ...){

######################Set up initial network(thinned/augmented)#############################
  nn = dim(coordinates(observations))[1]
  nGrd = dim(coordinates(predGrid))[1]
  netPtsInit = observations
  if (action == "del"){
    sn=sample(nn)
    delPts = observations[which((sn)<(min(sn)+nDiff)),]
    netPts = observations[which((sn)>=(min(sn)+nDiff)),]
    addPts = NULL
#    netPts = rbind(netPts,addPts)
    nn = dim(coordinates(netPts))[1]
    crit1 = calculateMukv(observations = netPts, predGrid = predGrid, model = model, ...) # initial design

    res=ssaMap(candidates, predGrid, model,   
                max_points_shift = 1, # maximum number of points to move in SSA (should be 1)
                maxShiftFactorX = 0.2, # maximum shift for points to move in SSA as percentage of total study area in X direction
                minShiftFactorX = 0, # minimum shift for points to move in SSA as percentage of total study area in X direction
                maxShiftFactorY = 0.2, # maximum shift for points to move in SSA as percentage of total study area in Y direction
                minShiftFactorY = 0, # maximum shift for points to move in SSA as percentage of total study area in Y direction
                start_p = 0.2,
                countMax = 200, netPts, addPts, delPts, crit1, nn, action, nDiff, netPtsInit, nr_iterations, plotOptim, ...)
    return(res)  
    }
  if (action == "add"){
    observations = SpatialPoints(observations)
    sn=sample(nGrd)
    delPts = NULL
    addPts = predGrid[which(sn<=nDiff),]
    netPts = rbind(observations,addPts)
    nn = dim(coordinates(netPts))[1]
    crit1 = calculateMukv(observations = netPts, predGrid = predGrid, model = model, ...) # initial design

    res=ssaMap(candidates, predGrid, model,  
                max_points_shift = 1, # maximum number of points to move in SSA (should be 1)
                maxShiftFactorX = 0.2, # maximum shift for points to move in SSA as percentage of total study area in X direction
                minShiftFactorX = 0, # minimum shift for points to move in SSA as percentage of total study area in X direction
                maxShiftFactorY = 0.2, # maximum shift for points to move in SSA as percentage of total study area in Y direction
                minShiftFactorY = 0, # maximum shift for points to move in SSA as percentage of total study area in Y direction
                start_p = 0.2,
                countMax = 200, netPts, addPts, delPts, crit1, nn, action, nDiff, netPtsInit, nr_iterations, plotOptim, ...)  
    return(res)
    }
  }


