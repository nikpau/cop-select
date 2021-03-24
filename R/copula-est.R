# Initialize dynamic working directory (automatically set if you use RStudio)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Due to the lack of my code being a whole package, I will load the Monte-Carlo estimation functions manually
maxlikelihoodEnv <- new.env()

source("aux-methods/aux-methods.R", local = maxlikelihoodEnv)
source("mc-estimation/iTAu_mc.R", local = maxlikelihoodEnv)
source("mc-estimation/iBeta_mc.R", local = maxlikelihoodEnv)
source("mc-estimation/ml_mc.R", local = maxlikelihoodEnv)
source("cop-selection/mcCopSelection.R", local = maxlikelihoodEnv)
attach(maxlikelihoodEnv, name = "Est-Functions")


#Load Packages (Most of them are already required inside the 
#               respective functions, but for a better overview I load them again

require(VineCopula)
require(stats4)
require(gsl)
require(microbenchmark)
require(ggplot2)

############################################################################################

#MONTE CARLO SIMULATION

############################################################################################

# Monte-Carlo parameter estimation for the Gaussian Copula with different estimation methods and different sample sizes.
#Maximum likelihood: Sample size: 30,50,100 | 1000 repetitions each | Rho from -0.9 to 0.9 in 0.1 steps
gaussianML30 <- ML_MC(30,1000,1,-.9,.9,.1,mse = T)
gaussianML50 <- ML_MC(50,1000,1,-.9,.9,.1,mse = T)
gaussianML100 <- ML_MC(100,1000,1,-.9,.9,.1,mse = T)


#Inversion of Kendall's Tau: Sample size: 30,50,100 | 1000 repetitions each | Rho from -0.9 to 0.9 in 0.1 steps
gaussianiTau30 <- inverseTau_MC(30,1000,1,-.9,.9,.1,mse = T)
gaussianiTau50 <- inverseTau_MC(50,1000,1,-.9,.9,.1,mse = T)
gaussianiTau100 <- inverseTau_MC(100,1000,1,-.9,.9,.1,mse = T)


#Inversion of Blomqvist's Beta: Sample size: 30,50,100 | 1000 repetitions each | Rho from -0.9 to 0.9 in 0.1 steps
gaussianiBeta30 <- inverseBeta_MC(30,1000,1,-.9,.9,.1,mse = T)
gaussianiBeta50 <- inverseBeta_MC(50,1000,1,-.9,.9,.1,mse = T)
gaussianiBeta100 <- inverseBeta_MC(100,1000,1,-.9,.9,.1,mse = T)


# Monte-Carlo parameter estimation for the t-Copula (df = 4,mse = T) with different estimation methods and different sample sizes.
#Maximum likelihood: Sample size: 30,50,100 | 1000 repetitions each | Rho from -0.9 to 0.9 in 0.1 steps
tML30 <- ML_MC(30,1000,2,-.9,.9,.1,df = 4,mse = T)
tML50 <- ML_MC(50,1000,2,-.9,.9,.1,df = 4,mse = T)
tML100 <- ML_MC(100,1000,2,-.9,.9,.1,df = 4,mse = T)

#Inversion of Kendall's Tau: Sample size: 30,50,100 | 1000 repetitions each | Rho from -0.9 to 0.9 in 0.1 steps
tiTau30 <- inverseTau_MC(30,1000,2,-.9,.9,.1,df = 4,mse = T)
tiTau50 <- inverseTau_MC(50,1000,2,-.9,.9,.1,df = 4,mse = T)
tiTau100 <- inverseTau_MC(100,1000,2,-.9,.9,.1,df = 4,mse = T)

#Inversion of Blomqvist's Beta: Sample size: 30,50,100 | 1000 repetitions each | Rho from -0.9 to 0.9 in 0.1 steps
tiBeta30 <- inverseBeta_MC(30,1000,2,-.9,.9,.1,df = 4,mse = T)
tiBeta50 <- inverseBeta_MC(50,1000,2,-.9,.9,.1,df = 4,mse = T)
tiBeta100 <- inverseBeta_MC(100,1000,2,-.9,.9,.1,df = 4,mse = T)

# Monte-Carlo parameter estimation for the Clayton Copula with different estimation methods and different sample sizes.
#Maximum likelihood: Sample size: 30,50,100 | 1000 repetitions each | Theta from 1 to 10 in 0.5 steps
claytonML30 <- ML_MC(30,1000,3,1,10,.5,mse = T)
claytonML50 <- ML_MC(50,1000,3,1,10,.5,mse = T)
claytonML100 <- ML_MC(100,1000,3,1,10,.5,mse = T)

#Inversion of Kendall's Tau: Sample size: 30,50,100 | 1000 repetitions each | Theta from 1 to 10 in 0.5 steps
claytoniTau30 <- inverseTau_MC(30,1000,3,1,10,.5,mse = T)
claytoniTau50 <- inverseTau_MC(50,1000,3,1,10,.5,mse = T)
claytoniTau100 <- inverseTau_MC(100,1000,3,1,10,.5,mse = T)

#Inversion of Blomqvist's Beta: Sample size: 30,50,100 | 1000 repetitions each | Theta from 1 to 10 in 0.5 steps
claytoniBeta30 <- inverseBeta_MC(30,1000,3,1,10,.5,mse = T)
claytoniBeta50 <- inverseBeta_MC(50,1000,3,1,10,.5,mse = T)
claytoniBeta100 <- inverseBeta_MC(100,1000,3,1,10,.5,mse = T)

# Monte-Carlo parameter estimation for the Gumbel Copula with different estimation methods and different sample sizes.
#Maximum likelihood: Sample size: 40,50,100 | 1000 repetitions each | Theta from 1 to 10 in 0.5 steps
gumbelML30 <- ML_MC(30,1000,4,1,10,.5,mse = T)
gumbelML50 <- ML_MC(50,1000,4,1,10,.5,mse = T)
gumbelML100 <- ML_MC(100,1000,4,1,10,.5,mse = T)

#Inversion of Kendall's Tau: Sample size: 40,50,100 | 1000 repetitions each | Theta from 1 to 10 in 0.5 steps
gumbeliTau30 <- inverseTau_MC(30,1000,4,1,10,.5,mse = T)
gumbeliTau50 <- inverseTau_MC(50,1000,4,1,10,.5,mse = T)
gumbeliTau100 <- inverseTau_MC(100,1000,4,1,10,.5,mse = T)

#Inversion of Blomqvist's Beta: Sample size: 40,50,100 | 1000 repetitions each | Theta from 1 to 10 in 0.5 steps
gumbeliBeta30 <- inverseBeta_MC(30,1000,4,1,10,.5,mse = T)
gumbeliBeta50 <- inverseBeta_MC(50,1000,4,1,10,.5,mse = T)
gumbeliBeta100 <- inverseBeta_MC(100,1000,4,1,10,.5,mse = T)

# Monte-Carlo parameter estimation for the Frank Copula with different estimation methods and different sample sizes.
#Maximum likelihood: Sample size: 30,50,100 | 1000 repetitions each | Theta from 1 to 10 in 0.5 steps
frankML30 <- ML_MC(30,1000,5,1,10,.5,mse = T)
frankML50 <- ML_MC(50,1000,5,1,10,.5,mse = T)
frankML100 <- ML_MC(100,1000,5,1,10,.5,mse = T)

#Inversion of Kendall's Tau: Sample size: 30,50,100 | 1000 repetitions each | Theta from 1 to 10 in 0.5 steps
frankiTau30 <- inverseTau_MC(30,1000,5,1,10,.5,mse = T)
frankiTau50 <- inverseTau_MC(50,1000,5,1,10,.5,mse = T)
frankiTau100 <- inverseTau_MC(100,1000,5,1,10,.5,mse = T)

#Inversion of Blomqvist's Beta: Sample size: 30,50,100 | 1000 repetitions each | Theta from 1 to 10 in 0.5 steps
frankiBeta30 <- inverseBeta_MC(30,1000,5,1,10,.5,mse = T)
frankiBeta50 <- inverseBeta_MC(50,1000,5,1,10,.5,mse = T)
frankiBeta100 <- inverseBeta_MC(100,1000,5,1,10,.5,mse = T)

######################################################################################
#Plotting
#######################################################################################

options(scipen = 20)

#Function for automated plotting
plot_my_stuff <- function(ml,itau,ibeta, title ="Title",ylim,ystep, ctype){

        xseq <- double()

        if(ctype == "e"){
                xseq <- seq(-.9,.9,.1)
        }
        else if(ctype == "a"){
                xseq <- seq(1,10,.5)
        }

        plot(ml$Param,ml$mse, type = "b", pch = 2, lwd = 3,lty = 1, col = "blue4",
             xaxt = "n", yaxt = "n", ylim = c(0,ylim),
             main = title, ylab="MSE",xlab="Parameter", cex.lab = 1.5)

        lines(itau$Param,itau$mse, type = "b", pch = 3, lwd = 3,lty = 2, col = "green4")

        lines(ibeta$Param,ibeta$mse, type = "b", pch = 4, lwd = 3,lty = 4, col = "red4")

        axis(side = 1, at = c(xseq), cex.axis = 1.5)
        axis(side = 2, las = 0, at = c(seq(0,ylim,ystep)), cex.axis = 1.5)

        grid(lwd = 1.5)

        legend("top", legend = c("ML","iTau","iBeta"), col = c("blue4","green4", "red4"), 
               lty = c(1,2,4),pch = c(2,3,4), cex = 1.5, lwd = 2)

}

plot_my_stuff(frankML30,frankiTau30,frankiBeta30,"Frank copula (n = 30)",2,.5, "a")

#######################################################################################
#Copula selection
#######################################################################################


# Simulate arbitrary Copula observations from 5 families (Gaussian(1),t(2),Clayton(3),Gumbel(4),Frank(5))
# See '?BiCopSim' for help
data <- BiCopSim(100,2,.1,4)
u1 <- data[,1]
u2 <- data[,2]


#Input arbitrary copula data of same length and get the best fitting copula back
copulaSelection(u1,u2,T)#function can be found in /aux-methods.R

#######################################################################################
#######################################################################################
#                         Monte-Carlo copula selection selection                      #

mcCop <- mcCopSelection(100,100,5,1,10,.5)


#Plotting for the selection process. Fig 3.4 / 3.5 in the paper
plot(mcCop[,1],mcCop[,2], type = "b", pch = 2, lwd = 3,lty = 1, col = "blue4",
     xaxt = "n", yaxt = "n", ylim = c(0,1),
     main = "Frank Copula (n = 100) ", ylab="Accuracy",xlab="Parameter", cex.lab = 1.5)

lines(mcCop[,1],mcCop[,3], type = "b", pch = 3, lwd = 3,lty = 2, col = "green4")

lines(mcCop[,1],mcCop[,4], type = "b", pch = 4, lwd = 3,lty = 4, col = "red4")

axis(side = 1, at = c(seq(1,10,.5)), cex.axis = 1.5)
axis(side = 2, las = 0, at = c(seq(0,1,.1)), cex.axis = 1.5)

grid(lwd = 1.5)

legend("bottomright", legend = c("ML","iTau","iBeta"), col = c("blue4","green4", "red4"), 
       lty = c(1,2,4),pch = c(2,3,4), cex = 1.5, lwd = 2)

#######################################################################################
#######################################################################################
#                                 Selection Bechmarking                               #

g <- copulaSelectionBenchmark(100,10,1,-.9,.9,.1)#function can be found in /aux-methods.R
t <- copulaSelectionBenchmark(100,10,2,-.9,.9,.1)
c <- copulaSelectionBenchmark(100,10,3,1,10,.5)
gu <- copulaSelectionBenchmark(100,10,4,1,10,.5)
f <- copulaSelectionBenchmark(100,10,5,1,10,.5)


#######################################################################################
#######################################################################################
#                                   BENCHMARKING                                      #

mbm_gauss <- microbenchmark( "ml" = {ML_MC(100,1000,1,-.9,.9,.1)},
                            "iTau" = {inverseTau_MC(100,1000,1,-.9,.9,.1)},
                            "iBeta" = {inverseBeta_MC(100,1000,1,-.9,.9,.1)}, 
                            times = 5, 
                            control = list(warmup = 1))


mbm_t <- microbenchmark( "ml" = {ML_MC(100,1000,2,-.9,.9,.1,4)},
                        "iTau" = {inverseTau_MC(100,1000,2,-.9,.9,.1,4)},
                        "iBeta" = {inverseBeta_MC(100,1000,2,-.9,.9,.1,4)}, 
                        times = 5,
                        control = list(warmup = 1))

mbm_clayton <- microbenchmark( "ml" = {ML_MC(100,1000,3,1,10,.5)},
                              "iTau" = {inverseTau_MC(100,1000,3,1,10,.5)},
                              "iBeta" = {inverseBeta_MC(100,1000,3,1,10,.5)}, 
                              times = 5, 
                              control = list(warmup = 1))

mbm_gumbel <- microbenchmark( "ml" = {ML_MC(100,1000,4,1,10,.5)},
                             "iTau" = {inverseTau_MC(100,1000,4,1,10,.5)},
                             "iBeta" = {inverseBeta_MC(100,1000,4,1,10,.5)}, 
                             times = 5, 
                             control = list(warmup = 1))

mbm_frank <- microbenchmark( "ml" = {ML_MC(100,1000,5,1,10,.5)},
                            "iTau" = {inverseTau_MC(100,1000,5,1,10,.5)},
                            "iBeta" = {inverseBeta_MC(100,1000,5,1,10,.5)}, 
                            times = 5, 
                            control = list(warmup = 1))
