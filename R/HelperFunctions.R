DFselection = function(Y, x, C) {
  
  CVgl = rep(NA, 3)
  G = dim(x)[2]
  n = dim(x)[1]
  start = proc.time()
  pc = 1
  if (dim(C)[2] > 1) {
    pc = dim(C)[2]
  }
  
  ## Run group lasso for 2 degrees of freedom
  try({
    mg = 1
    
    X = matrix(NA, nrow=n, ncol=(G*mg) + choose(G,2)*(mg^2))
    
    for (g in 1 : G) {
      X[,mg*(g-1) + 1] = x[,g]
    }
    
    counter = 1
    for (g1 in 2 : G) {
      for (g2 in 1 : (g1 - 1)) {
        splineTemp1 = x[,g1]
        splineTemp2 = x[,g2]
        
        designX = matrix(NA, n, mg^2)
        
        tempY = splineTemp1*splineTemp2
        tempMod = stats::lm(tempY ~ splineTemp1 + splineTemp2)
        designX[,1] = tempMod$residuals
        
        X[,((G*mg) + (counter-1)*(mg^2) + 1) : ((G*mg) + counter*(mg^2))] = designX
        counter = counter + 1
      }
    }
    
    penaltyFactor = c(rep(1, dim(X)[2]),  rep(0,pc))
    
    ## Group lasso
    cvGGL = glmnet::cv.glmnet(x=cbind(X, C), y=Y, penalty.factor=penaltyFactor)
    
    CVgl[1] = cvGGL$cvm[which.min(cvGGL$cvm)]
    
  })
  
  ## Run group lasso for 2 degrees of freedom
  try({
    mg = 2
    
    X = matrix(NA, nrow=n, ncol=(G*mg) + choose(G,2)*(mg^2))
    
    for (g in 1 : G) {
      splineTemp = x[,g]
      X[,mg*(g-1) + 1] = x[,g]
      for (m in 2 : mg) {
        splineTemp = cbind(splineTemp, x[,g]^m)
        tempY = splineTemp[,m]
        tempX = X[,(mg*(g-1) + 1):(mg*(g-1) + m - 1)]
        modX = stats::lm(tempY ~ tempX)
        X[,mg*(g-1) + m] = modX$residuals
      }
    }
    
    counter = 1
    for (g1 in 2 : G) {
      for (g2 in 1 : (g1 - 1)) {
        splineTemp1 = x[,g1]
        splineTemp2 = x[,g2]
        
        for (m in 2 : mg) {
          splineTemp1 = cbind(splineTemp1, x[,g1]^m)
          splineTemp2 = cbind(splineTemp2, x[,g2]^m)
        }
        
        designX = matrix(NA, n, mg^2)
        counter2 = 1
        for (m1 in 1 : mg) {
          for (m2 in 1 : mg) {
            tempY = splineTemp1[,m1]*splineTemp2[,m2]
            tempMod = stats::lm(tempY ~ splineTemp1 + splineTemp2)
            designX[,counter2] = tempMod$residuals
            counter2 = counter2 + 1
          }
        }
        
        X[,((G*mg) + (counter-1)*(mg^2) + 1) : ((G*mg) + counter*(mg^2))] = designX
        counter = counter + 1
      }
    }
    
    penaltyFactor = c(rep(1, dim(X)[2]),  rep(0,pc))
    
    ## Group lasso
    cvGGL = glmnet::cv.glmnet(x=cbind(X, C), y=Y, penalty.factor=penaltyFactor)
    
    CVgl[2] = cvGGL$cvm[which.min(cvGGL$cvm)]
    
  })
  
  ## Run group lasso for 2 degrees of freedom
  try({
    mg = 3
    
    X = matrix(NA, nrow=n, ncol=(G*mg) + choose(G,2)*(mg^2))
    
    for (g in 1 : G) {
      splineTemp = x[,g]
      X[,mg*(g-1) + 1] = x[,g]
      for (m in 2 : mg) {
        splineTemp = cbind(splineTemp, x[,g]^m)
        tempY = splineTemp[,m]
        tempX = X[,(mg*(g-1) + 1):(mg*(g-1) + m - 1)]
        modX = stats::lm(tempY ~ tempX)
        X[,mg*(g-1) + m] = modX$residuals
      }
    }
    
    counter = 1
    for (g1 in 2 : G) {
      for (g2 in 1 : (g1 - 1)) {
        splineTemp1 = x[,g1]
        splineTemp2 = x[,g2]
        
        for (m in 2 : mg) {
          splineTemp1 = cbind(splineTemp1, x[,g1]^m)
          splineTemp2 = cbind(splineTemp2, x[,g2]^m)
        }
        
        designX = matrix(NA, n, mg^2)
        counter2 = 1
        for (m1 in 1 : mg) {
          for (m2 in 1 : mg) {
            tempY = splineTemp1[,m1]*splineTemp2[,m2]
            tempMod = stats::lm(tempY ~ splineTemp1 + splineTemp2)
            designX[,counter2] = tempMod$residuals
            counter2 = counter2 + 1
          }
        }
        
        X[,((G*mg) + (counter-1)*(mg^2) + 1) : ((G*mg) + counter*(mg^2))] = designX
        counter = counter + 1
      }
    }
    
    penaltyFactor = c(rep(1, dim(X)[2]),  rep(0,pc))
    
    ## Group lasso
    cvGGL = glmnet::cv.glmnet(x=cbind(X, C), y=Y, penalty.factor=penaltyFactor)
    
    CVgl[3] = cvGGL$cvm[which.min(cvGGL$cvm)]
    
  })
  dof = which.min(CVgl)
  end = proc.time()
  timet = end[3] - start[3]
  l = list("df" = dof, "time" = timet)
  return(l)
  
}