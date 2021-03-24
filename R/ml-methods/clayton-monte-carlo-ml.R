# Monte-Carlo Simulation for parameter estimation of the clayton-copula.
# For each parameter step 1000 times, N copula samples are drawn and thus 1000 parameters are estimated.
# Will return a data frame containing all 1000 estimates per row and all respective parameter steps in the columns

claytonMC1000 <- function(N,parFrom,parTo,stepsize){

        clayton_mc_ml <- function(samplesize,par){

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

                        u1u2 <- BiCopSim(samplesize,3,par)
                        u1 <- u1u2[,1]
                        u2 <- u1u2[,2]

                        clayton_llik <- function(theta){
                                -sum(log((1+theta)*((u1*u2)^(-1-theta))*((-1+(u1^(-theta))+(u2^(-theta)))^(-2-(1/theta)))))
                        }

                        result[i] <- optim(5,clayton_llik, method= "L-BFGS-B", lower =.001, upper = 20)$par

                }
                #Return result
                result
        }

        coll_result <- data.frame(matrix(ncol = length(seq(parFrom,parTo,stepsize)), nrow = 1000)) #Init overall result matrix
        names_result <- sprintf("Theta = %s",seq(parFrom,parTo,stepsize)) #Define column names
        names(coll_result) <- names_result #Set column names
        colnum <- 1 #Column counter

        for (i in seq(parFrom,parTo,stepsize)) {

                #For every parameter step: Simulate N observations 1000 times and estimate one parameter each. 
                #Write every estimate per parameter step in one column. Do this for every step.
                coll_result[,colnum] <- clayton_mc_ml(N,i)
                colnum <- colnum + 1
        }
        coll_result
}
