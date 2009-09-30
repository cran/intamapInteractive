############################################################################
# Spatial simulated annealling with one criterion. 
# This function calls on 1 related functions (calculateMukv) 
# and returns criterionIterf (criterion for all iterations)
# Stopping criterion is given by a number of times 
# without improvement in search of a better design (countMax)
############################################################################

ssaMap = function(candidates, predGrid, model, max_points_shift, maxShiftFactorX, minShiftFactorX, 
                  maxShiftFactorY, minShiftFactorY, start_p, countMax, 
                  netPts, addPts, delPts, crit1, nn, action, nDiff, netPtsInit, nr_iterations, plotOptim,...) {

  largest_poly <- sapply(slot(candidates, "polygons"), function(x) slot(x, "plotOrder"))[1] # largest polygon
  studyareacoords <- coordinates(candidates@polygons[[1]]@Polygons[[largest_poly]]) 

  count = 0
# settings for simulated annealing
  x_bounds <- bbox(candidates)[1, ]
  y_bounds <- bbox(candidates)[2, ]
  x_extent <- x_bounds[2] - x_bounds[1]
  y_extent <- y_bounds[2] - y_bounds[1]
  max_shift_x <- maxShiftFactorX * x_extent # maximum shift in x-direction <Jan-Willem hoger>
  min_shift_x <- minShiftFactorX * x_extent  # minimum shift in x-direction
  max_shift_y <- maxShiftFactorY * y_extent  # maximum shift in y-direction <Jan-Willem hoger, to 0.50>
  min_shift_y <- minShiftFactorY * y_extent  # minimum shift in y-direction

  nr_designs <- 1  # counter for number of accepted designs

  oldpoints <- netPts  # save initial design because new design may be rejected
  criterionInitial <- crit1
  oldcriterion <- criterionInitial  # also save current criterion
  oldDelPoints <- delPts

  criterionIterf <- NULL

  for (k in 1:100000){
    
# scenario of deletion
    if (action == "del") {    

      selected_shifts_del <- sample(nDiff)  
      selected_shifts_net <- sample(nn)
      oldDelPt = oldDelPoints[which(selected_shifts_del<=max_points_shift),] 
      newDelPt = oldpoints[which(selected_shifts_net<=max_points_shift),]

      delPts=rbind(oldDelPoints[which(selected_shifts_del>max_points_shift),],newDelPt)
      netPts=rbind(oldpoints[which(selected_shifts_net>max_points_shift),],oldDelPt)
      criterion = calculateMukv(observations = netPts, predGrid = predGrid, model = model, ...)

      p = runif(1) # to allow accepting an inferior design
      criterionIterf <- c(criterionIterf,oldcriterion)            
      if (criterion <= oldcriterion){
        oldpoints = netPts
        netPts = netPts
        oldDelPoints = delPts
        delPts = delPts
        oldcriterion = criterion
        nr_designs = nr_designs+1
        count = 0
        cat("No improvement for",count,"iterations ")
        cat("[Will stop at",countMax,"iterations with no improvement]", "\n")
      } else if (criterion > oldcriterion & p <= (start_p*exp(-10*k/nr_iterations))){
        oldpoints = netPts
        netPts = netPts
        oldDelPoints = delPts
        delPts = delPts
        oldcriterion = criterion
        nr_designs = nr_designs+1
        count = count + 1
        cat("No improvement for",count,"iterations ")
        cat("[Will stop at",countMax,"iterations with no improvement]", "\n")
        } else{
        criterion = oldcriterion
        olddelPt = oldDelPoints
        netPts = oldpoints
        nr_designs = nr_designs
        count = count + 1
        cat("No improvement for",count,"iterations ")
        cat("[Will stop at",countMax,"iterations with no improvement]", "\n")
      }
      netPts <<- netPts
      if (count < countMax) 
        {
        if (plotOptim){
          if (!is.na(candidates))
            {                     
            plot(candidates)
            points(oldpoints, col = 1, pch = 19, cex = 0.7)
            points(oldDelPoints, col = 2, pch = "X", cex = 1.2)
            title('Simulated annealing', xlab=(paste("Criterion = ", round(criterion, digits=5))), ylab=(paste("Iterations = ", k)));
            } else {
            plot(oldpoints, col = 1, pch = 19, cex = 0.7)
            points(oldDelPoints, col = 2, pch = "X", cex = 1.2)
            title('Simulated annealing', xlab=(paste("Criterion = ", round(criterion, digits=5))), ylab=(paste("Iterations = ", k)));
            }
          }
        } else {break}
    }

# scenario of addition
    if (action=="add"){
      oldpoints = coordinates(oldpoints)
      netPts = as.data.frame(coordinates(netPts))
      selected_shifts <- c(rep(max_points_shift+1,nn-nDiff),sample(seq(1,nDiff)))
      for (i in 1:nn){
        if (selected_shifts[i] > max_points_shift){
          selected_shifts[i] = 0} else {
          selected_shifts[i] = 1
        } # no shift for this point SJM... these points
      }
      arraypos<-which.max(selected_shifts)
      infield <- 0
      while (infield<1){
        x_shift <- max_shift_x-k/nr_iterations*(max_shift_x-min_shift_x)  # possible shift in x-direction decreases linearly
        y_shift <- max_shift_y-k/nr_iterations*(max_shift_y-min_shift_y)  # possible shift in y-direction decreases linearly
        netPts[,1] <- oldpoints[,1] + x_shift*((as.matrix(rep(-1,nn))+2*as.matrix(runif(nn)))*as.matrix(selected_shifts))  # true shift in x-direction
        netPts[,2] <- oldpoints[,2] + y_shift*((as.matrix(rep(-1,nn))+2*as.matrix(runif(nn)))*as.matrix(selected_shifts))  # true shift in y_direction   
        net=netPts
        coordinates(net)=~x+y
        while(length(zerodist(net))[1]>0){
          netPts[,1] <- oldpoints[,1] + x_shift*((as.matrix(rep(-1,nn))+2*as.matrix(runif(nn)))*as.matrix(selected_shifts))  # true shift in x-direction
          netPts[,2] <- oldpoints[,2] + y_shift*((as.matrix(rep(-1,nn))+2*as.matrix(runif(nn)))*as.matrix(selected_shifts))  # true shift in y_direction
          net=netPts
          coordinates(net)=~x+y
          }
        pointx <- netPts[arraypos, 1]
        pointy <- netPts[arraypos, 2]
        pointxy <- netPts[arraypos,]
        tmpi <- coordinates(candidates@polygons[[1]]@Polygons[[1]])
        infield <- (point.in.polygon(pointx, pointy, tmpi[,1], tmpi[,2]))
        }
      criterion = calculateMukv(observations = netPts, predGrid = predGrid, model = model, ...)
      netPts <- as.data.frame(netPts) # need as dataframe for oldpoints
      p = runif(1) # to allow accepting an inferior design
      criterionIterf <- c(criterionIterf,oldcriterion)            
      if (criterion <= oldcriterion){
        oldpoints = netPts
        netPts = netPts
        oldcriterion = criterion
        nr_designs = nr_designs+1
        count = 0
        cat("No improvement for",count,"iterations ")
        cat("[Will stop at",countMax,"iterations with no improvement]", "\n")
        } else if (criterion > oldcriterion & p <= (start_p*exp(-10*k/nr_iterations))){
        oldpoints = netPts
        netPts = netPts
        oldcriterion = criterion
        nr_designs = nr_designs+1
        count = count + 1
        cat("No improvement for",count,"iterations ")
        cat("[Will stop at",countMax,"iterations with no improvement]", "\n")
        } else {
        criterion = oldcriterion
        netPts = oldpoints
        count = count + 1
        cat("No improvement for",count,"iterations ")
        cat("[Will stop at",countMax,"iterations with no improvement]", "\n")
        }
      netPts <<- netPts
      if (count < countMax) 
        {              
        if (plotOptim == TRUE){
          plot(candidates)
          points(netPtsInit, col=1, pch = 19, cex = 0.7)  
          points(netPts[(nn-nDiff+1):nn,],  col = "green", pch = 19)
          points(pointxy,  col = 2, pch = 19)
          title('Simulated annealing', xlab=(paste("Criterion = ", round(criterion, digits=5))), ylab=(paste("Iterations = ", k)))
          } 
        } else {break}
      }
    }
#  return(list(netPts = netPts, criterionIterf = criterionIterf, itNumber=length(criterionIterf)))
  
  if (!inherits(netPts,"Spatial")) {
    netPts = data.frame(x = netPts[,1], y= netPts[,2])
    coordinates(netPts) = ~x+y
  }
  return(netPts)
  }


  
  