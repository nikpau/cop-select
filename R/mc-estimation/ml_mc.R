# Estimation of copula parameter using classic maximum likelihood
#
# [1] For each cycle out of K, a sample of size N for a chosen family is generated.
#      Families are (1) Gaussian
#                   (2) Students t
#                   (3) Clayton
#                   (4) Gumbel
#                   (5) Frank
# [3] Likelihood functions for every pdf are created and parameters are estimates using optim.
# [4] Parameters for specified sequence are written into a data frame.

ML_MC <- function(N,K,family,parStart,parEnd,stepsize,df = NA,mse = F){
  
  #Load required packages for Copula Simulation
  require(VineCopula)
  require(stats)
  
  #Initialize empty result vector
  result <- vector()
  
  #Define separate power function (R-base operator can be difficult to handle with longer expressions)
  pow <- function(base,exponent){
    base^exponent
  }
  
  #Parameter estimation based on selected copula family
  mlParEst <- function(samplesize, family, par){
    
    if(family == 1){ #Gaussian
      
      for (i in 1:K) {
        
        u1u2 <- BiCopSim(samplesize,family,par)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        guess <- runif(1,1,10)
        
        gaussian_llik <- function(rho){
          e1 <- qnorm(u1,0,1)
          e2 <- qnorm(u2,0,1)
          -sum(log(1/sqrt(1-pow(rho,2))*exp((pow(e1,2)+pow(e2,2))/2+(2*rho*e1*e2-pow(e1,2)-pow(e2,2))/(2*(1-pow(rho,2))))))
        }
        
        result[i] <- optim(guess,gaussian_llik, method = "L-BFGS-B", lower = -.95, upper = .95)$par
        
      }
      
    }
    
    else if(family == 2){ #Student t
      
      for (i in 1:K) {
        
        u1u2 <- BiCopSim(samplesize,family,par,df)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        guess <- runif(1,1,10)
        
        t_llik <- function(rho){
          e1 <- qt(u1,df)
          e2 <- qt(u2,df)
          -sum(log((gamma((df+2)/2)/gamma(df/2))/(df*pi*sqrt(1-pow(rho,2))*dt(e1,df)*dt(e2,df))*pow(1+(pow(e1,2)+pow(e2,2)-2*rho*e1*e2)/(df*(1-pow(rho,2))),-(df+2)/2)))
        }
        
        result[i] <- optim(guess,t_llik, method = "L-BFGS-B", lower = -.98, upper = .98)$par
        
      }
      
    }
    
    else if(family == 3){ #Clayton
      
      for (i in 1:K) {
        
        u1u2 <- BiCopSim(samplesize,family,par)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        guess <- runif(1,1,10)
        
        clayton_llik <- function(theta){
          -sum(log((1+theta)*((u1*u2)^(-1-theta))*((-1+(u1^(-theta))+(u2^(-theta)))^(-2-(1/theta)))))
        }
        
        result[i] <- optim(guess,clayton_llik, method= "L-BFGS-B", lower =.01, upper = 20)$par
        
      }
      
    }
    
    else if(family == 4){ #Gumbel
      
      for (i in 1:K) {
        
        u1u2 <- BiCopSim(samplesize,family,par)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        guess <- runif(1,1,10)
        
        gumbel_llik <- function(theta){
          e1 <- pow(-log(u1),theta)+pow(-log(u2),theta)
          -sum(-pow(e1,1.0/(theta))+(2.0/(theta)-2.0)*log(e1)+(theta-1.0)*log(log(u1)*log(u2))-log(u1*u2)+log1p((theta-1.0)*pow(e1,-1.0/(theta))))
        }
        
        result[i] <- optim(guess,gumbel_llik, method = "L-BFGS-B", lower = 1, upper = 20)$par
        
      }
      
    }
    
    else if(family == 5){ #Frank
      
      for (i in 1:K) {
        
        u1u2 <- BiCopSim(samplesize,family,par)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        guess <- runif(1,1,10)
        
        frank_llik <- function(theta){
          -sum(log((theta*(exp(theta)-1.0)*exp(theta*u1+theta*u2+theta))/(exp(theta*u2+theta*u1)-exp(theta*u2+theta)-exp(theta*u1+theta)+exp(theta))^2))
        }
        
        result[i] <- optim(guess,frank_llik, method = "L-BFGS-B", lower =.01, upper = 20)$par
        
      }
      
    }
    
    result
    
  }
  
  coll_result <- data.frame(matrix(ncol = length(seq(parStart,parEnd,stepsize)), nrow = K))
  names_result <- sprintf("Theta = %s",seq(parStart,parEnd,stepsize))
  names(coll_result) <- names_result
  colnum <- 1
  
  for (i in seq(parStart,parEnd,stepsize)) {
    
    
    coll_result[,colnum] <- mlParEst(N,family,i)
    colnum <- colnum + 1
  }
  
  # if TRUE: Directly calculate the mse for every 1000 trails per cycle
  if(mse == T){
    mse_result <- MSE(coll_result,parStart,parEnd,stepsize)
    mse_result
  }
  else
    coll_result
  
}