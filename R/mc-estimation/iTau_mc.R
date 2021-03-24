# Estimation of copula parameter using the inverse of Kendall's Tau.
#
# [1] For each cycle a sample of size N for each family is generated.
#      Families are (1) Gaussian
#                   (2) Students t
#                   (3) Clayton
#                   (4) Gumbel
#                   (5) Frank (numerical inversion).
# [3] tau is estimated for each sample and the corresponding copula parameter is calculated.
# [4] Parameters for specified sequence are written into a data frame.


inverseTau_MC <- function(N,K,family,parStart,parEnd,stepsize,df = NA,mse = F){ #df is only used if family == 2
  
  require(VGAM)
  require(VineCopula)
  require(stats)
  require(gsl)
  
  result <- vector()

  
  parEst <- function(samplesize, family,par){
    
  
  if(family == 1){ #Gaussian
    
    for (i in 1:K) {
      
      u1u2 <- BiCopSim(samplesize,family,par)
      u1 <- u1u2[,1]
      u2 <- u1u2[,2]
      
      tau <- kendall.tau(u1,u2,T)
      
      result[i] <- sin(tau*(pi/2))
      
    }
  }
    else if(family == 2){ #Students t
      
      for (i in 1:K) {
        
        u1u2 <- BiCopSim(samplesize,family,par,df) #df = degrees of freedom for students t
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        
        tau <- kendall.tau(u1,u2,T)
        
        result[i] <- sin(tau*(pi/2))
        
      }
    }
    else if(family == 3){ # Clayton
      
      for (i in 1:K) {
        
        u1u2 <- BiCopSim(samplesize,family,par)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        
        tau <- kendall.tau(u1,u2,T)
        
        result[i] <- 2*(tau/(1-tau))
        
      }
    }
      else if(family == 4){ # Gumbel
       
         for (i in 1:K) {
          
          u1u2 <- BiCopSim(samplesize,family,par)
          u1 <- u1u2[,1]
          u2 <- u1u2[,2]
          
          tau <- kendall.tau(u1,u2,T)
          
          result[i] <- 1/(1-tau)
        
        }
      
      }
    
    else if(family == 5){ # Frank
      
      for (i in 1:K) {
        
        u1u2 <- BiCopSim(samplesize,family,par)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        
        tau <- kendall.tau(u1,u2,T)
        
        result[i] <- FrankTauToPar(tau)
        
      }
      
    }
    result
  }
  
  coll_result <- data.frame(matrix(ncol = length(seq(parStart,parEnd,stepsize)), nrow = K))
  names_result <- sprintf("Parameter = %s",seq(parStart,parEnd,stepsize))
  names(coll_result) <- names_result
  colnum <- 1
  
  for (i in seq(parStart,parEnd,stepsize)) {
    
    
    coll_result[,colnum] <- parEst(N,family,i)
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
