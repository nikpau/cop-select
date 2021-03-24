# Monte-Carlo Simulation for parameter estimation of the frank-copula.
# For each parameter step 1000 times, N copula samples are drawn and thus 1000 parameters are estimated.
# Will return a data frame containing all 1000 estimates per row and all respective parameter steps in the columns

frankMC1000 <- function(N,parFrom,parTo,stepsize){
  
  frank_mc_ml <- function(samplesize,par){
    
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
    
    #For every loop: [1] Generate copula observations
    #                [2] Plug values into log copula density, calculate negative sum for optimization w.r.t theta
    #                [3] Find argmax(llik) and write the estimate for theta into the result vector
    for (i in 1:1000) {
      
      u1u2 <- BiCopSim(samplesize,5,par)
      u1 <- u1u2[,1]
      u2 <- u1u2[,2]
      
      frank_llik <- function(theta){
        -sum(log((theta*(exp(theta)-1.0)*exp(theta*u1+theta*u2+theta))/(exp(theta*u2+theta*u1)-exp(theta*u2+theta)-exp(theta*u1+theta)+exp(theta))^2))
      }
      
      result[i] <- optim(5,frank_llik, method = "L-BFGS-B", lower =.001, upper = 20)$par
      
    }
    #Return result
    result
  }
  
  coll_result <- data.frame(matrix(ncol = length(seq(parFrom,parTo,stepsize)), nrow = 1000))
  names_result <- sprintf("Theta = %s",seq(parFrom,parTo,stepsize))
  names(coll_result) <- names_result
  colnum <- 1
  
  for (i in seq(parFrom,parTo,stepsize)) {
    
    
    coll_result[,colnum] <- frank_mc_ml(N,i)
    colnum <- colnum + 1
  }
  coll_result
}