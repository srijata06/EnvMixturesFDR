#' Estimate effects of exposures and identify important ones using debiased lasso.
#'
#'This function estimates the effects of the exposures using debiased lasso 
#'and identifies the important exposures and their interactions adjusting for confounders while
#'controlling FDR at a given threshold. 
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
#' @return A list of values containing important main effects and interactions, estimate of the overall exposure effect and the 95\% CI of the estimated response.
#'
#' @export



DBL = function(Y, x, C, xnew, q=0.2){
  
  n = dim(x)[1]
  n2 = dim(xnew)[1]
  G = ncol(x)
  ng = G + choose(G,2)
  pc = 1
  if (!is.null(dim(C)[2])) {
    pc = dim(C)[2]
  }
  
  X1 = list()
  Xnew1 = list()
  conf_int = matrix(NA,nrow(xnew),2)
  for (mgg in 1:3) {
    X = matrix(NA, nrow = nrow(x), ncol=(G*mgg) + choose(G,2)*(mgg^2))
    Xnew = matrix(NA, nrow = nrow(xnew), ncol=(G*mgg) + choose(G,2)*(mgg^2))
    
    for (g in 1 : G) {
      X[,(mgg*(g-1) + 1)] = x[,g]
      Xnew[,mgg*(g-1) + 1] = xnew[,g]
      for (m in 2 : mgg){
        X[,(mgg*(g-1) + m) ] = x[,g]^m
        Xnew[,(mgg*(g-1) + m) ] = xnew[,g]^m
      }
    }
    
    counter = 1
    for (g1 in 1 : (G-1)) {
      for (g2 in (g1+1) : G) {
        designX = matrix(NA, n, mgg^2)
        designXnew = matrix(NA, n2, mgg^2)
        counter2 = 1
        for(m1 in 1:mgg){
          for(m2 in 1:mgg){
            
            designX[,counter2] = X[,(mgg*(g1-1) + m1)] * X[,(mgg*(g2-1) + m2)]
            designXnew[,counter2] = Xnew[,(mgg*(g1-1) + m1)] * Xnew[,(mgg*(g2-1) + m2)]
            counter2 = counter2 +1
            
          }
        }
        X[,((G*mgg) + (counter-1)*(mgg^2) + 1) : ((G*mgg) + counter*(mgg^2))] = designX
        Xnew[,((G*mgg) + (counter-1)*(mgg^2) + 1) : ((G*mgg) + counter*(mgg^2))] = designXnew
        
        counter = counter + 1
        
      }
    }
    X1[[mgg]] = X 
    Xnew1[[mgg]] = Xnew
  } 
  DFselect = DFselection(Y=Y, x=x, C=C)
  mg = DFselect$df
  
  X2 = X1[[mg]]
  Xnew2 = Xnew1[[mg]]
  np = mg*G + choose(G,2)*(mg^2)
  groups = c(rep(1:G, each=mg),rep((G+1):(choose(G,2)+G),each=mg**2))
  
  ## Now scale design matrices
  mX = apply(X2, 2, mean)
  sX = apply(X2, 2, sd)
  
  for (jj in 1 : dim(X2)[2]) {
    X2[,jj] = (X2[,jj] - mX[jj]) / (sX[jj])
    Xnew2[,jj] = (Xnew2[,jj] - mX[jj]) / (sX[jj])
  }
  
  
  ####### create design matrix with both X and C #######################
  X2 = cbind(X2, C)
  
  ###### Cross validation ##############################################
  ccGG = glmnet::cv.glmnet(x=X2, y=Y, penalty.factor = c(rep(1, dim(Xnew2)[2]), rep(0, pc)))
  GG = glmnet::glmnet(x=X2, y=Y, penalty.factor = c(rep(1, dim(Xnew2)[2]), rep(0, pc)), 
              lambda=ccGG$lambda.min)
  
  ######### debiasing procedure ##########################################
  p0 = dim(X2)[2]
  cf=matrix(NA,p0,p0)
  t=rep(0,p0)
  
  sigma = (t(X2) %*% X2)/n
  cvKeep = rep(NA, 10)
  
  for (k in 1:p0){
    if (k <= 10) {
      cvmod=glmnet::cv.glmnet(x=X2[,-k],y=X2[,k],intercept=FALSE, nfolds=5)
      cvKeep[k] = cvmod$lambda.1se
      mod=glmnet::glmnet(x=X2[,-k],y=X2[,k],lambda=cvmod$lambda.1se,intercept=FALSE) 
    } else {
      mod=glmnet::glmnet(x=X2[,-k],y=X2[,k],lambda=median(cvKeep),intercept=FALSE) 
    }
    cf[k,k] = 1
    cf[k,-k] = as.numeric(-mod$beta)
    
    t[k]= (1/n)*(sum((X2[,k]-(X2[,-k]%*%mod$beta))**2))+0.5*cvmod$lambda.1se*((sum(abs(cf[k,])))-1)
  }
  
  t1=diag(t,p0)
  theta=(solve(t1))%*%cf
  
  ########## debiased lasso estimate ##############################################
  
  debiased = as.numeric(GG$beta) + theta %*% t(X2) %*% (Y - (X2 %*% as.numeric(GG$beta)))/n
  
  ########## asymptotic variance of debiased lasso ################################
  
  varBeta = (theta %*% sigma %*% t(theta)) / n
  debiasednew = debiased[-((p0-pc+1):p0)]
  varBetanew = varBeta[-((p0-pc+1):p0),-((p0-pc+1):p0)]
  
  ######### FDR Controlling Procedure ############################################
  Tvec = rep(NA, ng)
  
  for (j in 1 : ng) {
    wcols = which(groups == j)
    Tvec[j] = debiasednew[wcols] %*% solve(varBetanew[wcols,wcols]) %*% debiasednew[wcols]
  }
  
  upperLimit = stats::qchisq(0.999,df = mg^2)
  upperLimit1 = stats::qchisq(0.999,df = mg)
  upperLimit2 = stats::qchisq(0.999,df = mg^2)
  
  tSeq = seq(0, upperLimit, length=1000)
  fdrSeq = seq(0, upperLimit, length=1000)
  fdrSeq1 = seq(0, upperLimit1, length=1000)
  fdrSeq2 = seq(0, upperLimit2, length=1000)
  
  
  for (j in 1 : length(tSeq)) {
    tempT = tSeq[j]
    fdrSeq[j] = (G*(1 - stats::pchisq(tempT, df=mg))+ choose(G,2)*(1-stats::pchisq(tempT,df=(mg^2)))) / max(1, length(which(Tvec > tempT)))
    fdrSeq1[j] = (G*(1 - stats::pchisq(tempT, df=mg))) / max(1, length(which(Tvec[1:G] > tempT)))
    fdrSeq2[j] = (choose(G,2)*(1-stats::pchisq(tempT,df=mg^2))) / max(1, length(which(Tvec[(G+1):ng] > tempT)))
    
  }
  
  Tfinal = tSeq[which(fdrSeq <= q)[1]]
  Tfinal1 = tSeq[which(fdrSeq1 <= q)[1]]
  Tfinal2 = tSeq[which(fdrSeq2 <= q)[1]]
  
  intercept = mean(Y - (X2 %*% debiased))
  estimate = (Xnew2 %*% debiasednew) +intercept
  variance = Xnew2 %*% varBetanew %*% t(Xnew2)
  
  conf_int[,1] = estimate - 1.96*sqrt(diag(variance))
  conf_int[,2] = estimate + 1.96*sqrt(diag(variance))
  
  ic_1 = c(which(Tvec[1:G] > Tfinal1))
  ic_2 = c(G+(which(Tvec[(G+1):ng] > Tfinal2)))
  ic_global = c(which(Tvec > Tfinal))
  
  importantMain1 = which(Tvec[1:G] > Tfinal1)
  importantMain2 = which(Tvec[1:G] > Tfinal)
  
  importantInt1 = importantInt2 = data.frame(X1 = NULL, X2 = NULL)
  counter = G + 1
  for (g1 in 1 : (G-1)) {
    for (g2 in (g1+1) : G) {
      if (counter %in% ic_2) importantInt1 = rbind(importantInt1, c(g1,g2))
      if (counter %in% ic_global) importantInt2 = rbind(importantInt2, c(g1,g2))
      counter = counter + 1
    }
  }
  
  if (nrow(importantInt1) > 0) {
    names(importantInt1) = c("X1", "X2")   
  }
  
  if (nrow(importantInt2) > 0) {
    names(importantInt2) = c("X1", "X2")   
  }
  
  return(list("importantMain" = importantMain1, 
              "importantInt" = importantInt1,
              "est" = estimate,"lower" = conf_int[,1],
              "upper" = conf_int[,2]))
  
}
