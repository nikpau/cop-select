gaussianMC1000 <- function(N,parFrom,parTo,stepsize){
  
  gaussian_mc_ml <- function(samplesize,par){
    
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
    #                [2] Plug values into log copula density, calculate negative sum for optimization w.r.t rho
    #                [3] Find argmax(llik) and write the estimate for rho into the result vector
    for (i in 1:1000) {
      
      u1u2 <- BiCopSim(samplesize,1,par)
      u1 <- u1u2[,1]
      u2 <- u1u2[,2]
      
      gaussian_llik <- function(rho){
        e1 <- qnorm(u1,0,1)
        e2 <- qnorm(u2,0,1)
        -sum(log(1/sqrt(1-pow(rho,2))*exp((pow(e1,2)+pow(e2,2))/2+(2*rho*e1*e2-pow(e1,2)-pow(e2,2))/(2*(1-pow(rho,2))))))
      }
      
      result[i] <- optim(.5,gaussian_llik, method = "L-BFGS-B", lower = -.98, upper = .98)$par
      
    }
    #Return result
    result
  }
  
  coll_result <- data.frame(matrix(ncol = length(seq(parFrom,parTo,stepsize)), nrow = 1000))
  names_result <- sprintf("Theta = %s",seq(parFrom,parTo,stepsize))
  names(coll_result) <- names_result
  colnum <- 1
  
  for (i in seq(parFrom,parTo,stepsize)) {
    
    
    coll_result[,colnum] <- gaussian_mc_ml(N,i)
    colnum <- colnum + 1
  }
  coll_result
}
