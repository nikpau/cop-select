#Function that repeatedly simulates copula observations from a known copula to check how
#often the different estimators guess right

mcCopSelection <- function(N,K,family,parStart,parEnd,stepsize){

        #Load required packages
        require(VGAM)
        require(VineCopula)
        require(stats)
        require(gsl)

        #Function for selecting a copula per estimator for a given set of data
        multipleCopulaSelect <- function(samplesize,family,par){

                fam_names <- c("Gaussian","t","Clayton","Gumbel","Frank")

                #Counter for right estimations
                true_counter_ml <- 0
                true_counter_iTau <- 0
                true_counter_iBeta <- 0

                for (i in 1:K) {

                        if(family == 2){
                                u1u2 <- BiCopSim(samplesize,family,par,4)
                                u1 <- u1u2[,1]
                                u2 <- u1u2[,2]
                        }
                        else{
                                u1u2 <- BiCopSim(samplesize,family,par)
                                u1 <- u1u2[,1]
                                u2 <- u1u2[,2]
                        }

                        #Estimate a copula based on observations.
                        sel <- copulaSelection(u1,u2)#function can be found in /aux-methods.R

                        #increase counter per estimator if it estimates the correct copula
                        ifelse((sel[[2]][1,1] =="ml" & sel[[2]][1,2] == fam_names[family]) == T,
                               true_counter_ml <- true_counter_ml + 1,
                               NA)

                        ifelse((sel[[2]][2,1] =="iTau" & sel[[2]][2,2] == fam_names[family]) == T,
                               true_counter_iTau <- true_counter_iTau + 1,
                               NA)

                        ifelse((sel[[2]][3,1] =="iBeta" & sel[[2]][3,2] == fam_names[family]) == T,
                               true_counter_iBeta <- true_counter_iBeta + 1,
                               NA)


                }
                #Calculate share of right guesses
                intermediate <- c(true_counter_ml/K,true_counter_iTau/K,true_counter_iBeta/K)
                intermediate

        }

        results <- data.frame(matrix(ncol = 3, nrow = length(seq(parStart,parEnd,stepsize))))
        names(results) <- c("ML","iTau","iBeta")

        rownum <- 1

        par <- data.frame(Parameter = seq(parStart,parEnd,stepsize))

        results <- cbind(par,results)

        for (i in seq(parStart,parEnd,stepsize)) {


                results[rownum,2] <- multipleCopulaSelect(N,family,i)[1]
                results[rownum,3] <- multipleCopulaSelect(N,family,i)[2]
                results[rownum,4] <- multipleCopulaSelect(N,family,i)[3]
                rownum <- rownum + 1
        }

        results

}
