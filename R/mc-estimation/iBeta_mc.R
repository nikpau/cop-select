# Estimation of copula parameter using the inverse of Blomqvist's Beta.
#
# [1] For each cycle a sample of size N for each family is generated.
#      Families are (1) Gaussian
#                   (2) Students t
#                   (3) Clayton
#                   (4) Gumbel
#                   (5) Frank (numerical inversion).
# [3] Beta is estimated for each sample and the corresponding copula parameter is calculated.
# [4] Parameters for specified sequence are written into a data frame.

inverseBeta_MC <- function(N,K,family,parStart,parEnd,stepsize,df = NA,mse = F){
  
  require(VGAM)
  require(VineCopula)
  require(stats)
  
  result <- vector()
  
  
  #Estimate Parameters per copula using the just defined Beta
  betaParEst <- function(samplesize,family,par){
    
    if(family == 1){ #Gaussian
      
      for(i in 1:K){
      u1u2 <- BiCopSim(samplesize,family,par)
      u1 <- u1u2[,1]
      u2 <- u1u2[,2]
      
      beta <- blomqvist.beta(u1,u2)
      result[i] <- sin(beta*(pi/2))
      }
    }
    
    else if(family == 2){
      
      for(i in 1:K){
        u1u2 <- BiCopSim(samplesize,family,par,df)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        
        beta <- blomqvist.beta(u1,u2)
        result[i] <- sin(beta*(pi/2))
      }
      
    }
    
    else if(family == 3){ #Clayton
      
      for(i in 1:K){
        u1u2 <- BiCopSim(samplesize,family,par)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        
        beta <- blomqvist.beta(u1,u2)
        result[i] <- ClaytonBetaToPar(beta)
      }
      
    }
    
    else if(family == 4){ #Gumbel
      
      
      
      for(i in 1:K){
        u1u2 <- BiCopSim(samplesize,family,par)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        
        beta <- blomqvist.beta(u1,u2)
        
        gumbel_beta <- function(beta){
          if(beta >= 0){
            -(log(2)/(log(-log(2)/log((beta+1)/4))))
          }
          else{
            1
          }
        }
        result[i] <- gumbel_beta(beta)
      }
      
    }
    
    else  if(family == 5) {
      
      for(i in 1:K){ #Frank
        u1u2 <- BiCopSim(samplesize,family,par)
        u1 <- u1u2[,1]
        u2 <- u1u2[,2]
        
        beta <- blomqvist.beta(u1,u2)
        result[i] <- FrankBetaToPar(beta)
      }
      
    }
    
    result
    
  }
  
  coll_result <- data.frame(matrix(ncol = length(seq(parStart,parEnd,stepsize)), nrow = K))
  names_result <- sprintf("Parameter = %s",seq(parStart,parEnd,stepsize))
  names(coll_result) <- names_result
  colnum <- 1
  
  for (i in seq(parStart,parEnd,stepsize)) {
    
    
    coll_result[,colnum] <- betaParEst(N,family,i)
    colnum <- colnum + 1
  }

  # if TRUE: Directly calculate the mse for every 1000 trails per cycle
  if(mse == T){
    mse_result <- MSE(coll_result,parStart,parEnd,stepsize) #function MSE found in /aux-methods.R
    mse_result
  }
  else
    coll_result
  
}
