gumbelMC1000 <- function(N,parFrom,parTo,stepsize){
  
  gumbel_mc_ml <- function(samplesize,par){
    
    #Check validity of parameters
    if(!(typeof(samplesize)%in% c("num","double", "int")))
      stop("samplesize is no integer")
    
    if(!(typeof(par)%in% c("num","double", "int")))
      stop("Parameter is no integer")
    
    #Load required packages for Copula Simulation
    require(VineCopula)
    require(stats)
    
    #Initialize empty result vector
    result <- vector()
    
    #Define separate power function (R-base operator can be difficult to handle with longer expressions)
    pow <- function(base,exponent){
      base^exponent
    }
    
    
    #For every loop: [1] Generate copula observations
    #                [2] Plug values into log copula density, calculate negative sum for optimization w.r.t theta
    #                [3] Find argmax(llik) and write the estimate for theta into the result vector
    for (i in 1:1000) {
      
      u1u2 <- BiCopSim(samplesize,4,par)
      u1 <- u1u2[,1]
      u2 <- u1u2[,2]
      
      gumbel_llik <- function(theta){
        e1 <- pow(-log(u1),theta)+pow(-log(u2),theta)
        -sum(-pow(e1,1.0/(theta))+(2.0/(theta)-2.0)*log(e1)+(theta-1.0)*log(log(u1)*log(u2))-log(u1*u2)+log1p((theta-1.0)*pow(e1,-1.0/(theta))))
      }
      
      result[i] <- optim(5,gumbel_llik, method = "L-BFGS-B", lower = 1, upper = 20)$par
      
    }
    #Return result
    result
  }
  
  coll_result <- data.frame(matrix(ncol = length(seq(parFrom,parTo,stepsize)), nrow = 1000))
  names_result <- sprintf("Theta = %s",seq(parFrom,parTo,stepsize))
  names(coll_result) <- names_result
  colnum <- 1
  
  for (i in seq(parFrom,parTo,stepsize)) {
    
    
    coll_result[,colnum] <- gumbel_mc_ml(N,i)
    colnum <- colnum + 1
  }
  coll_result
}