# EnvMixturesFDR

An R package for implementing the approaches in [Samanta and Antonelli (2021)](https://arxiv.org/abs/2103.10563) to simultaneously estimate 
the health effects of environmental mixtures and identify important exposures and interactions while
controlling FDR.  

### Installing EnvMixturesFDR

The first step is to install the library which can be done as follows:
```R
library(devtools)
install_github(repo = "srijata06/EnvMixturesFDR")
library(EnvMixturesFDR)
```  

### Data Generation

```R
n <- 500
G <- 10
pc <- 1
nnew <- 10

C = matrix(rnorm(n*pc), nrow=n) 
covsigma=matrix(NA,G,G)
      rho=0.5
      for(u in 1:G) {
        for(v in 1:G){
          covsigma[u,v]=rho**(abs(u-v))
        }
      }
      
X = mvtnorm::rmvnorm(n, mean=rep(0,G), sigma = covsigma)
Xnew = mvtnorm::rmvnorm(nnew, mean=rep(0,G), sigma = covsigma)

TrueF = function(x) {
    return(0.06*x[,7] + 0.10*(x[,8]^2)  - 0.04*x[,5] + 0.08*(x[,5]^2) + 0.25*x[,1]*x[,4]^2 + 0.18*x[,2]*x[,3])
}

Y = TrueF(X) + C + rnorm(n)
```  

This gives us data on `p = 10` exposures and 1 additional confounder variable (i,e, additional 
variable we adjust for example: age) and an outcome depending on these exposures and the confounder. 
The matrix Xnew is created to be used for inference (estimation) purposes, specifically to calculte
the MSE and confidence intervals of the estimated response.

### Using the functions

Now we describe how to implement the proposed approaches using the package. We have 3 main functions in this package, 
one each for an approach. For the two methods based on the knockoffs procedure we have the functions
`knockoffsSplit()` and `knockoffsFull()`. The former performs exposure/interaction selection using one 
half of the data and conditional inference using the other half of the data while the later uses 
the entire dataset for both the steps. The function `DBL()` performs estimation and variable selection
using debiased lasso. 

```R
kSplit_Results = knockoffsSplit(Y=Y, x=X, C=C, xnew=Xnew, q=0.25)        ##### default value of q is 0.2
kFull_Results = knockoffsFull(Y=Y, x=X, C=C, xnew=Xnew, q=0.25)          ##### default value of q is 0.2
DBLresults = DBL(Y=Y, x=X, C=C, xnew=Xnew, q=0.25)                       ##### default value of q is 0.2
```  

The functions `knockoffsSplit()`, `knockoffsFull()` and `DBL()` return important main effects, important interaction 
effects, the estimated outcome, computation time and confidence intervals for the estimated response. 
