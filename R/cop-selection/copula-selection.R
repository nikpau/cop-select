# Copula selection process based on maximum likelihood, inversion of Kendall's tau and Blomqvist's Beta.
#
# [1]Parameters are estimated with the three different estimators
# [2]AIC is calculated for each estimate (for the inversion methods, the estimated parameter is
#    plugged back into the log-likelihood function (optim step is skipped))
# [3]Lowest AIC with corresponding copula per estimator is returned
#
#
copulaSelection <- function(u1,u2, benchmark = F){

        #Check whether input values are in the interval [0,1]. If not, stop code execution, else continue.
        ifelse(u1 > 1 | u2 > 1 | u1 < 0 | u2 < 0,
               stop("Copula observations must be in the interval of [0,1]"),
               NA)

        #Load required packages
        require(VGAM)
        require(VineCopula)
        require(stats)
        require(gsl)




        #Initialize result vector
        result <- data.frame(EstMethod = c(rep("ml",5),rep("iTau",5),rep("iBeta",5)),
                             Family = rep(c("Gaussian","t","Clayton","Gumbel","Frank"),3),
                             AIC = c(rep(NA,5*3)), 
                             estPar = c(rep(NA,5*3)))

        #Power function
        pow <- function(base,exponent){
                base^exponent
        }

        #Calculate Kendall's Tau for data
        tau <- kendall.tau(u1,u2,T)

        #Calculate Blomqvist's Beta for data
        beta <- blomqvist.beta(u1,u2)

        #Load in log-likelihood functions for copulae
        lliks <- list(gaussian_llik <- function(rho){ #Gaussian
                              e1 <- qnorm(u1,0,1)
                              e2 <- qnorm(u2,0,1)
                              -sum(log(1/sqrt(1-pow(rho,2))*exp((pow(e1,2)+pow(e2,2))/2+(2*rho*e1*e2-pow(e1,2)-pow(e2,2))/(2*(1-pow(rho,2))))))
                             },
                             t_llik <- function(rho){
                                     e1 <- qt(u1,4)
                                     e2 <- qt(u2,4)
                                     -sum(log((gamma((4+2)/2)/gamma(4/2))/(4*pi*sqrt(1-pow(rho,2))*dt(e1,4)*dt(e2,4))*pow(1+(pow(e1,2)+pow(e2,2)-2*rho*e1*e2)/(4*(1-pow(rho,2))),-(4+2)/2)))
                             }, #4 degrees of freedom are set fix
                             clayton_llik <- function(theta){
                                     -sum(log((1+theta)*((u1*u2)^(-1-theta))*((-1+(u1^(-theta))+(u2^(-theta)))^(-2-(1/theta)))))
                             },
                             gumbel_llik <- function(theta){
                                     e1 <- pow(-log(u1),theta)+pow(-log(u2),theta)
                                     -sum(-pow(e1,1.0/(theta))+(2.0/(theta)-2.0)*log(e1)+(theta-1.0)*log(log(u1)*log(u2))-log(u1*u2)+log1p((theta-1.0)*pow(e1,-1.0/(theta))))
                             },
                             frank_llik <- function(theta){
                                     -sum(log((theta*(exp(theta)-1.0)*exp(theta*u1+theta*u2+theta))/(exp(theta*u2+theta*u1)-exp(theta*u2+theta)-exp(theta*u1+theta)+exp(theta))^2))
                             })

        #List of Kendall's Tau inverses for respective copulae
        iTauFuncs <- list(gauss_iTau <- function(tau){
                                  sin(tau*(pi/2))
                             },
                             t_iTau <- function(tau){
                                     sin(tau*(pi/2))
                             },
                             clayton_iTau <- function(tau){
                                     2*(tau/(1-tau))
                             },
                             gumbel_iTau <- function(tau){
                                     1/(1-tau)
                             },
                             frank_iTau <- function(tau){FrankTauToPar(tau)})

        #List of Blomqvist's Beta inverses for respective copulae
        iBetaFuncs <- list(gauss_beta <- function(beta){
                                   sin(beta*(pi/2))
                             },
                             t_beta <- function(beta){
                                     sin(beta*(pi/2))
                             },
                             clayton_beta <- function(beta){ClaytonBetaToPar(beta)},
                             gumbel_beta <- function(beta){
                                     -(log(2)/(log(-log(2)/log((beta+1)/4))))
                             },
                             frank_beta <- function(beta){FrankBetaToPar(beta)})


        #Random initial parameter for optim function  
        guess <- function(copulatype){

                if(copulatype == "a"){#Archimedean Copulae
                        runif(1,1,10)
                }
                else if(copulatype == "e"){#Elliptical Copulae
                        runif(1,-1,1)
                }
        }

        #Define Akaike Information Criterion
        AkInCr <- function(llik,value){
                2*lliks[[llik]](value)+2*1
        }

        # ML-Estimation  
        ml_results <- vector()

        for (i in 1:2) {
                ml_results[i] <- optim(guess("e"),lliks[[i]],method = "L-BFGS-B", lower = -.95, upper = .95)$par
        }

        for (i in 3:5) {
                ml_results[i] <- optim(guess("a"),lliks[[i]],method = "L-BFGS-B", lower = 1, upper = 20)$par
        }

        for (i in 1:5) {
                result$AIC[i] <- AkInCr(i,ml_results[i])
        }


        #iTau Estimation
        iTau_results <- vector()


        for (i in 1:5) {
                iTau_results[i] <- iTauFuncs[[i]](tau)
        }

        for (i in 1:5) {
                result$AIC[i+5] <- AkInCr(i,iTau_results[i])
        }

        #iBeta Estimation
        iBeta_results <- vector()


        for (i in 1:5){
                iBeta_results[i] <- iBetaFuncs[[i]](beta)
        }

        for (i in 1:5){
                result$AIC[i+10] <- AkInCr(i,iBeta_results[i])
        }

        # Benchmarking of the above expression. Expression will
        # be executed 100 times. Average computation time will be shown.
        if(benchmark == T){
                mbm <- microbenchmark("ml" = {ml_results <- vector()

                                      for (i in 1:2) {
                                              ml_results[i] <- optim(guess("e"),lliks[[i]],method = "L-BFGS-B", lower = -.95, upper = .95)$par
                                      }

                                      for (i in 3:5) {
                                              ml_results[i] <- optim(guess("a"),lliks[[i]],method = "L-BFGS-B", lower = 1, upper = 20)$par
                                      }

                                      for (i in 1:5) {
                                              result$AIC[i] <- AkInCr(i,ml_results[i])
                             }},
                             "iTau" = {   

                                     iTau_results <- vector()


                                     for (i in 1:5) {
                                             iTau_results[i] <- iTauFuncs[[i]](tau)
                                     }

                                     for (i in 1:5) {
                                             result$AIC[i+5] <- AkInCr(i,iTau_results[i])
                             }},
                             "iBeta" = {

                                     iBeta_results <- vector()


                                     for (i in 1:5){
                                             iBeta_results[i] <- iBetaFuncs[[i]](beta)
                                     }

                                     for (i in 1:5){
                                             result$AIC[i+10] <- AkInCr(i,iBeta_results[i])
                             }})
        }

        #write all parameter estimates in the last column of result vector
        for (i in 1:5){
                result$estPar[i] <- ml_results[i]
                result$estPar[i+5] <- iTau_results[i]
                result$estPar[i+10] <- iBeta_results[i]
        }

        # Split result vector again by estimation method
        agg_ml <- result[result$EstMethod == "ml",]
        agg_iTau <- result[result$EstMethod == "iTau",]
        agg_iBeta <- result[result$EstMethod == "iBeta",]

        agg <- data.frame(EstMethod = rep(NA,3),
                          Family = rep(NA,3),
                          AIC = rep(NA,3), estPar = rep(NA,3))

        #Only write the row with lowest AIC per estimation method
        agg[1,] <- agg_ml[which.min(agg_ml$AIC),]
        agg[2,] <- agg_iTau[which.min(agg_iTau$AIC),]
        agg[3,] <- agg_iBeta[which.min(agg_iBeta$AIC),]

        #return results
        if(benchmark == T){
                list(Call = deparse(sys.call()), SelectionResults = agg, Benchmark = mbm)
        }
        else
                list(Call = deparse(sys.call()), SelectionResults = agg)
}
