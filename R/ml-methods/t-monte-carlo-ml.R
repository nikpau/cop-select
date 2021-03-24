tMC1000 <- function(N,parFrom,parTo,stepsize,df){
  
  t_mc_ml <- function(samplesize,par,df){
    
    #Check validity of parameters
    if(!(typeof(samplesize)%in% c("num","double", "int")))
      stop("samplesize is no integer")
    
    if(!(typeof(par)%in% c("num","double", "int")))
      stop("Parameter is no integer")
    
    if(!(typeof(df)%in% c("num","double", "int")))
      stop("Degrees of freedom is no integer")
    
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
      
      u1u2 <- BiCopSim(samplesize,2,par,df)
      u1 <- u1u2[,1]
      u2 <- u1u2[,2]
      
      t_llik <- function(rho){
        e1 <- qt(u1,df)
        e2 <- qt(u2,df)
        -sum(log((gamma((df+2)/2)/gamma(df/2))/(df*pi*sqrt(1-pow(rho,2))*dt(e1,df)*dt(e2,df))*pow(1+(pow(e1,2)+pow(e2,2)-2*rho*e1*e2)/(df*(1-pow(rho,2))),-(df+2)/2)))
      }
      
      result[i] <- optim(.5,t_llik, method = "L-BFGS-B", lower = -.98, upper = .98)$par
      
    }
    #Return result
    result
  }
  
  coll_result <- data.frame(matrix(ncol = length(seq(parFrom,parTo,stepsize)), nrow = 1000))
  names_result <- sprintf("Theta = %s",seq(parFrom,parTo,stepsize))
  names(coll_result) <- names_result
  colnum <- 1
  
  for (i in seq(parFrom,parTo,stepsize)) {
    
    
    coll_result[,colnum] <- t_mc_ml(N,i,df)
    colnum <- colnum + 1
  }
  coll_result
}
