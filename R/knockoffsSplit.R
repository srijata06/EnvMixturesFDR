#' Identify important exposures and estimate their effects using knockoffs and data splitting.
#' 
#' This function identifies important exposures and their interactions and 
#' estimates their effects adjusting for confounders while also controlling FDR
#' at a given threshold. The method is based on the knockoffs approach and data splitting.
#'
#' @param Y                   The n by 1 outcome to be analyzed.
#' @param x                   An n by p matrix of exposures.
#' @param C                   An n by q matrix of additional covariates to adjust for.
#' @param xnew                A design matrix with p exposure values for estimation purpose.
#' @param q                   The false discovery rate threshold.
#'
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom stats lm rgamma runif rbeta cov var rnorm coef vcov
#' @importFrom knockoff create.second_order
#' @return A list of values with important main effects and interactions, estimates, 95\% CI and computation time.
#'
#' @export

knockoffsSplit = function(Y,x,C,xnew,q=0.2) {
  binary_alpha = q
  n = length(Y)
  G = dim(x)[2]
  start = proc.time()
  pc = 1
  if (!is.null(dim(C)[2])) {
    pc = dim(C)[2]
  }
  
  ## number of training samples
  n1 = floor(n/2)
  x1 = x[1:n1,]
  Y1 = Y[1:n1]
  
  if (dim(C)[2] == 1) {
    C1 = C[1:n1] 
  } else {
    C1 = C[1:n1,]
  }
  
  DFselect = DFselection(Y=Y, x=x, C=C)
  mg = DFselect$df
  time_df = DFselect$time
  
  X1 = matrix(NA, nrow=n1, ncol=(G*mg) + choose(G,2)*(mg^2))
  X1knock = matrix(NA, nrow=n1, ncol=(G*mg) + choose(G,2)*(mg^2))
  
  xk1 = create.second_order(X = x1, shrink=TRUE)
  
  for (g in 1 : G) {
    for (m in 1 : mg) {
      X1[,mg*(g-1) + m] = x1[,g]^m
      X1knock[,mg*(g-1) + m] = xk1[,g]^m
    }
  }
  
  counter = 1
  for (g1 in 2 : G) {
    for (g2 in 1 : (g1 - 1)) {
      tempX1 = matrix(NA, n1, mg)
      tempX2 = matrix(NA, n1, mg)
      
      for (m in 1 : mg) {
        tempX1[,m] = x1[,g1]^m
        tempX2[,m] = x1[,g2]^m
      }
      
      designX = matrix(NA, n1, mg^2)
      counter2 = 1
      for (m1 in 1 : mg) {
        for (m2 in 1 : mg) {
          tempY = (x1[,g1]^m1) * (x1[,g2]^m2)
          designX[,counter2] = lm(tempY ~ tempX1 + tempX2)$residuals
          counter2 = counter2 + 1
        }
      }
      
      X1[,((G*mg) + (counter-1)*(mg^2) + 1) : ((G*mg) + counter*(mg^2))] = designX
      
      counter = counter + 1
    }
  }
  
  counter = 1
  for (g1 in 2 : G) {
    for (g2 in 1 : (g1 - 1)) {
      tempX1 = matrix(NA, n1, mg)
      tempX2 = matrix(NA, n1, mg)
      
      for (m in 1 : mg) {
        tempX1[,m] = xk1[,g1]^m
        tempX2[,m] = xk1[,g2]^m
      }
      
      designX = matrix(NA, n1, mg^2)
      counter2 = 1
      for (m1 in 1 : mg) {
        for (m2 in 1 : mg) {
          tempY = (xk1[,g1]^m1) * (xk1[,g2]^m2)
          designX[,counter2] = lm(tempY ~ tempX1 + tempX2)$residuals
          counter2 = counter2 + 1
        }
      }
      
      X1knock[,((G*mg) + (counter-1)*(mg^2) + 1) : ((G*mg) + counter*(mg^2))] = designX
      
      counter = counter + 1
    }
  }
  
  ## Now scale design matrices
  mX = apply(X1, 2, mean)
  sX = apply(X1, 2, sd)
  
  for (jj in 1 : dim(X1)[2]) {
    X1[,jj] = (X1[,jj] - mX[jj]) / (sX[jj])
  }
  
  mXk = apply(X1knock, 2, mean)
  sXk = apply(X1knock, 2, sd)
  
  for (jj in 1 : dim(X1knock)[2]) {
    X1knock[,jj] = (X1knock[,jj] - mXk[jj]) / (sXk[jj])
  }
  
  ## create knockoff design matrices
  nGroups = G + choose(G, 2)
  
  XC1 = cbind(X1, X1knock, scale(C1))
  groups = c(rep(1:G, each=mg), rep(G + 1:(choose(G,2)), each=mg^2), 
             (G + choose(G, 2) + rep(1:G, each=mg)), 
             rep(G + choose(G, 2) + G + 1:(choose(G,2)), each=mg^2), 
             rep(G + choose(G, 2) + G + choose(G,2) + 1, pc))
  
  penaltyFactor = c(rep(1, 2*dim(X1)[2]),rep(0, pc))
  
  ## Cross validation
  ccGG = cv.glmnet(x=XC1, y=Y1, penalty.factor = penaltyFactor)
  GG = glmnet(x=XC1, y=Y1, 
              penalty.factor = penaltyFactor, 
              lambda=ccGG$lambda.min)
  
  ## For a given lambda, estimate the knockoff statistics
  tempBeta = as.numeric(GG$beta)
  
  Wj = rep(NA, nGroups)
  for (j in 1 : nGroups) {
    Wj[j] = sqrt(sum(tempBeta[groups == j]^2)) - sqrt(sum(tempBeta[groups == (j + nGroups)]^2))
  }
  
  WjFinal = Wj
  
  WjMain = WjFinal[1:G]
  WjInt = WjFinal[-c(1:G)]
  
  ## Do the rest of the process with the conservative model selection
  
  tVecMain = sort(abs(WjMain), decreasing=TRUE)
  FDPestMain = rep(NA, length(tVecMain))
  
  for (t in 1 : length(tVecMain)) {
    tempThresh = tVecMain[t]
    FDPestMain[t] = (1 + sum(WjMain <= -tempThresh)) / max(sum(WjMain >= tempThresh), 1)
  }
  
  finalThresholdsMain = tVecMain[which(FDPestMain < binary_alpha)]
  finalThresholdMain = finalThresholdsMain[length(finalThresholdsMain)]
  
  wGroupsMain = which(WjMain >= finalThresholdMain)
  
  
  tVecInt = sort(abs(WjInt), decreasing=TRUE)
  FDPestInt = rep(NA, length(tVecInt))
  
  for (t in 1 : length(tVecInt)) {
    tempThresh = tVecInt[t]
    FDPestInt[t] = (1 + sum(WjInt <= -tempThresh)) / max(sum(WjInt >= tempThresh), 1)
  }
  
  finalThresholdsInt = tVecInt[which(FDPestInt < binary_alpha)]
  finalThresholdInt = finalThresholdsInt[length(finalThresholdsInt)]
  
  wGroupsInt = G + which(WjInt >= finalThresholdInt)
  
  wGroups = c(wGroupsMain, wGroupsInt)
  wSelected = c(which(groups %in% wGroups), dim(X1)[2]+(1:pc))
  
  ############################################################################
  ##############################  Inference split #############################
  ############################################################################
  
  nnew = dim(xnew)[1]
  
  if (length(wGroups) == 0) {
    n2 = n - n1
    x2 = x[-c(1:n1),]
    Y2 = Y[-c(1:n1)]
    
    if (dim(C)[2] == 1) {
      C2 = C[-c(1:n1)]
    } else {
      C2 = C[-c(1:n1),]
    }
    
    mod = lm(Y2 ~ 1 + C2)
    
    est = as.numeric(rep(mod$coefficients[1], nnew))
    CIlower = as.numeric(rep(mod$coefficients[1] - 1.96*summary(mod)$coefficients[1,2], nnew))
    CIupper = as.numeric(rep(mod$coefficients[1] + 1.96*summary(mod)$coefficients[1,2], nnew))
  } else {
    n2 = n - n1
    x2 = x[-c(1:n1),]
    Y2 = Y[-c(1:n1)]
    
    if (dim(C)[2] == 1) {
      C2 = C[-c(1:n1)]
    } else {
      C2 = C[-c(1:n1),]
    }
    
    X2 = matrix(NA, nrow=n2, ncol=(G*mg) + choose(G,2)*(mg^2))
    Xnew = matrix(NA, nrow=nnew, ncol=(G*mg) + choose(G,2)*(mg^2))
    
    for (g in 1 : G) {
      for (m in 1 : mg) {
        X2[,mg*(g-1) + m] = x2[,g]^m
        Xnew[,mg*(g-1) + m] = xnew[,g]^m
      }
    }
    
    counter = 1
    for (g1 in 2 : G) {
      for (g2 in 1 : (g1 - 1)) {
        tempX1 = matrix(NA, n2, mg)
        tempX2 = matrix(NA, n2, mg)
        
        tempX1new = matrix(NA, nnew, mg)
        tempX2new = matrix(NA, nnew, mg)
        
        for (m in 1 : mg) {
          tempX1[,m] = x2[,g1]^m
          tempX2[,m] = x2[,g2]^m
          
          tempX1new[,m] = xnew[,g1]^m
          tempX2new[,m] = xnew[,g2]^m
        }
        
        designX = matrix(NA, n2, mg^2)
        designXnew = matrix(NA, nnew, mg^2)
        
        counter2 = 1
        for (m1 in 1 : mg) {
          for (m2 in 1 : mg) {
            tempY = (x2[,g1]^m1) * (x2[,g2]^m2)
            tempMod = lm(tempY ~ tempX1 + tempX2)
            designX[,counter2] = tempMod$residuals
            
            tempYnew = (xnew[,g1]^m1) * (xnew[,g2]^m2)
            designXnew[,counter2] = tempYnew - cbind(rep(1,nnew), tempX1new, tempX2new) %*%
              tempMod$coefficients
            counter2 = counter2 + 1
          }
        }
        
        X2[,((G*mg) + (counter-1)*(mg^2) + 1) : ((G*mg) + counter*(mg^2))] = designX
        Xnew[,((G*mg) + (counter-1)*(mg^2) + 1) : ((G*mg) + counter*(mg^2))] = designXnew
        
        counter = counter + 1
      }
    }
    
    groups = c(rep(1:G, each=mg), rep(G + 1:(choose(G,2)), each=mg^2), rep(G + choose(G, 2) + 1, pc))
    
    ## Now scale design matrices
    mX = apply(X2, 2, mean)
    sX = apply(X2, 2, sd)
    
    for (jj in 1 : dim(X2)[2]) {
      X2[,jj] = (X2[,jj] - mX[jj]) / (sX[jj])
      Xnew[,jj] = (Xnew[,jj] - mX[jj]) / (sX[jj])
    }
    
    XC2 = cbind(X2, C2)
    
    Xinf = XC2[,wSelected]
    
    mod = lm(Y2 ~ Xinf)
    
    coefEst = mod$coefficients[-c((length(mod$coefficients) - pc + 1):
                                    length(mod$coefficients))]
    
    varBeta = vcov(mod)[-c((length(mod$coefficients) - pc + 1):
                             length(mod$coefficients)),
                        -c((length(mod$coefficients) - pc + 1):
                             length(mod$coefficients))]
    
    ## keep just those coefficients associated with X only and were chosen
    XnewTemp = cbind(rep(1, dim(Xnew)[1]), 
                     Xnew[,wSelected[-which(wSelected > G*mg + choose(G,2)*mg^2)]])
    
    ## Look at inference for the global function
    est = XnewTemp %*% coefEst
    variance = XnewTemp %*% varBeta %*% t(XnewTemp)
    CIlower = est - 1.96*sqrt(diag(variance))
    CIupper = est + 1.96*sqrt(diag(variance))
  }
  
  importantMain = which(1:G %in% wGroups)
  importantInt = data.frame(X1 = NULL, X2 = NULL)
  
  counter = G + 1
  for (g1 in 2 : G) {
    for (g2 in 1 : (g1 - 1)) {
      if (counter %in% wGroups) importantInt = rbind(importantInt, c(g2,g1))
      counter = counter + 1
    }
  }
  if (nrow(importantInt) > 0) {
    names(importantInt) = c("X1", "X2")   
  }
  
  
  conservative = list(importantMain = importantMain,
                      importantInt = importantInt,
                      est = est,
                      lower = CIlower,
                      upper = CIupper)
  
  end = proc.time()
  timet = end[3] - start[3] + time_df 
  
  
  
  ## Do the rest of the process with the conservative model selection
  
  tVecMain = sort(abs(WjMain), decreasing=TRUE)
  FDPestMain = rep(NA, length(tVecMain))
  
  for (t in 1 : length(tVecMain)) {
    tempThresh = tVecMain[t]
    FDPestMain[t] = (sum(WjMain <= -tempThresh)) / max(sum(WjMain >= tempThresh), 1)
  }
  
  finalThresholdsMain = tVecMain[which(FDPestMain < binary_alpha)]
  finalThresholdMain = finalThresholdsMain[length(finalThresholdsMain)]
  
  wGroupsMain = which(WjMain >= finalThresholdMain)
  
  
  tVecInt = sort(abs(WjInt), decreasing=TRUE)
  FDPestInt = rep(NA, length(tVecInt))
  
  for (t in 1 : length(tVecInt)) {
    tempThresh = tVecInt[t]
    FDPestInt[t] = (sum(WjInt <= -tempThresh)) / max(sum(WjInt >= tempThresh), 1)
  }
  
  finalThresholdsInt = tVecInt[which(FDPestInt < binary_alpha)]
  finalThresholdInt = finalThresholdsInt[length(finalThresholdsInt)]
  
  wGroupsInt = G + which(WjInt >= finalThresholdInt)
  
  
  wGroups = c(wGroupsMain, wGroupsInt)
  wSelected = c(which(groups %in% wGroups), dim(X1)[2]+(1:pc))
  
  ############################################################################
  ##############################  Inference split #############################
  ############################################################################
  
  nnew = dim(xnew)[1]
  
  if (length(wGroups) == 0) {
    n2 = n - n1
    x2 = x[-c(1:n1),]
    Y2 = Y[-c(1:n1)]
    
    if (dim(C)[2] == 1) {
      C2 = C[-c(1:n1)]
    } else {
      C2 = C[-c(1:n1),]
    }
    
    mod = lm(Y2 ~ 1 + C2)
    
    est = as.numeric(rep(mod$coefficients[1], nnew))
    CIlower = as.numeric(rep(mod$coefficients[1] - 1.96*summary(mod)$coefficients[1,2], nnew))
    CIupper = as.numeric(rep(mod$coefficients[1] + 1.96*summary(mod)$coefficients[1,2], nnew))
  } else {
    n2 = n - n1
    x2 = x[-c(1:n1),]
    Y2 = Y[-c(1:n1)]
    
    if (dim(C)[2] == 1) {
      C2 = C[-c(1:n1)]
    } else {
      C2 = C[-c(1:n1),]
    }
    
    X2 = matrix(NA, nrow=n2, ncol=(G*mg) + choose(G,2)*(mg^2))
    Xnew = matrix(NA, nrow=nnew, ncol=(G*mg) + choose(G,2)*(mg^2))
    
    for (g in 1 : G) {
      for (m in 1 : mg) {
        X2[,mg*(g-1) + m] = x2[,g]^m
        Xnew[,mg*(g-1) + m] = xnew[,g]^m
      }
    }
    
    counter = 1
    for (g1 in 2 : G) {
      for (g2 in 1 : (g1 - 1)) {
        tempX1 = matrix(NA, n2, mg)
        tempX2 = matrix(NA, n2, mg)
        
        tempX1new = matrix(NA, nnew, mg)
        tempX2new = matrix(NA, nnew, mg)
        
        for (m in 1 : mg) {
          tempX1[,m] = x2[,g1]^m
          tempX2[,m] = x2[,g2]^m
          
          tempX1new[,m] = xnew[,g1]^m
          tempX2new[,m] = xnew[,g2]^m
        }
        
        designX = matrix(NA, n2, mg^2)
        designXnew = matrix(NA, nnew, mg^2)
        
        counter2 = 1
        for (m1 in 1 : mg) {
          for (m2 in 1 : mg) {
            tempY = (x2[,g1]^m1) * (x2[,g2]^m2)
            tempMod = lm(tempY ~ tempX1 + tempX2)
            designX[,counter2] = tempMod$residuals
            
            tempYnew = (xnew[,g1]^m1) * (xnew[,g2]^m2)
            designXnew[,counter2] = tempYnew - cbind(rep(1,nnew), tempX1new, tempX2new) %*%
              tempMod$coefficients
            counter2 = counter2 + 1
          }
        }
        
        X2[,((G*mg) + (counter-1)*(mg^2) + 1) : ((G*mg) + counter*(mg^2))] = designX
        Xnew[,((G*mg) + (counter-1)*(mg^2) + 1) : ((G*mg) + counter*(mg^2))] = designXnew
        
        counter = counter + 1
      }
    }
    
    groups = c(rep(1:G, each=mg), rep(G + 1:(choose(G,2)), each=mg^2), rep(G + choose(G, 2) + 1, pc))
    
    ## Now scale design matrices
    mX = apply(X2, 2, mean)
    sX = apply(X2, 2, sd)
    
    for (jj in 1 : dim(X2)[2]) {
      X2[,jj] = (X2[,jj] - mX[jj]) / (sX[jj])
      Xnew[,jj] = (Xnew[,jj] - mX[jj]) / (sX[jj])
    }
    
    XC2 = cbind(X2, C2)
    
    Xinf = XC2[,wSelected]
    
    mod = lm(Y2 ~ Xinf)
    
    coefEst = mod$coefficients[-c((length(mod$coefficients) - pc + 1):
                                    length(mod$coefficients))]
    
    varBeta = vcov(mod)[-c((length(mod$coefficients) - pc + 1):
                             length(mod$coefficients)),
                        -c((length(mod$coefficients) - pc + 1):
                             length(mod$coefficients))]
    
    ## keep just those coefficients associated with X only and were chosen
    XnewTemp = cbind(rep(1, dim(Xnew)[1]), 
                     Xnew[,wSelected[-which(wSelected > G*mg + choose(G,2)*mg^2)]])
    
    ## Look at inference for the global function
    est = XnewTemp %*% coefEst
    variance = XnewTemp %*% varBeta %*% t(XnewTemp)
    CIlower = est - 1.96*sqrt(diag(variance))
    CIupper = est + 1.96*sqrt(diag(variance))
  }
  
  importantMain = which(1:G %in% wGroups)
  importantInt = data.frame(X1 = NULL, X2 = NULL)
  
  counter = G + 1
  for (g1 in 2 : G) {
    for (g2 in 1 : (g1 - 1)) {
      if (counter %in% wGroups) importantInt = rbind(importantInt, c(g2,g1))
      counter = counter + 1
    }
  }
  if (nrow(importantInt) > 0) {
    names(importantInt) = c("X1", "X2")   
  }
  
  regular = list(importantMain = importantMain,
                 importantInt = importantInt,
                 est = est,
                 lower = CIlower,
                 upper = CIupper, time = timet)
  
  l = list(regular = regular,
           conservative = conservative)
  
  return(l)
  
}

