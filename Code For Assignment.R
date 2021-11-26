library(nimble)
library(igraph)
library(coda)
library(R6)

set.seed(234)
#####Question 2#####

par(mfrow=c(1,2))
hist(rbeta(10000,21,19))
hist(rbeta(10000,82,2))


#####Question 4####
#creates list to store samples
PPVsamples <- c()


for (i in 1:10000){
  #takes samples from posteriors
  Se <- rbeta(1,21,19)
  Sp <- rbeta(1,82,2)
  logitpi <- rnorm(1,0,10^2) 
  pi <- exp(logitpi)/(1+exp(logitpi))
  #plugs samples into equation for PPV
  PPV <- (Se*pi)/(Se*pi+((1-pi)*(1-Sp)))
  PPVsamples <- c(PPVsamples,PPV) 
}

#finds mean of PPV samples
mean(PPVsamples)

par(mfrow=c(1,1))
#plots histogram of PPV samples
hist(PPVsamples)

#removes variables that could possibly interfere with subsequent questions
rm(i,pi,logitpi)

#



#####Question 5#####

load("~/MT4531 Bayesian Inference/Assignment/Data/dataQ5.Rdata")


####Model With Both Covariates specification#### 

#adds a column for prevalence to original data 
dataQ5$pi <- dataQ5$ncarriers/dataQ5$ni

# Specifies model containing both covariates
carriermodelcode <- nimbleCode({
  #Specifying the likelihood 
  for (i in 1:N){
    ncarriers[i] ~ dbinom(pi[i],ni[i])
  }
  for (i in 1:N){
    logit(pi[i]) <- beta0 + betaWh*((X_hp[i] - mean(X_hp[1:N]))/sd(X_hp[1:N])) + betaHi*((X_hi[i] - mean(X_hi[1:N]))/sd(X_hi[1:N]))
  }
  #Priors
  beta0 ~ dnorm(0, 0.04)
  betaWh ~ dnorm(0, 0.04)
  betaHi ~ dnorm(0, 0.04)
})

#model constants
carriermodel_constants <- list(N = 15, ni = dataQ5$ni, X_hp = dataQ5$X_hp, X_hi = dataQ5$X_hi)

#data to be used
carriermodel_data <- list(ncarriers = dataQ5$ncarriers,pi = dataQ5$pi)

#set initial values of model
initvalues <- list(beta0 = 0, betaWh = 0, betaHi = 0)

#builds the model
carriers <- nimbleModel(code = carriermodelcode, name = 'carriers', constants = carriermodel_constants, data = carriermodel_data, inits = initvalues) 

#compiles model
Compcarriers <- compileNimble(carriers)

# sets up the quantities to be monitored. 
carriersconfig <- configureMCMC(carriers, enableWAIC = FALSE, monitors = c('beta0','betaWh','betaHi'), print = TRUE)  

# builds the MCMC algorithm
carriersMCMC <- buildMCMC(carriersconfig)

# compile the MCMC chain 
CompcarriersMCMC <- compileNimble(carriersMCMC, project = carriers)

#Sets initial values for MCMC chains
carriersMCMC_Inits <- list(list(beta0=-5, betaWh=-5, betaHi=-5), 
                          list(beta0=5, betaWh=5, betaHi=5))



####Testing To find Burn-in Length etc####

#runs MCMC with 2 chains, small number of iterations and no burnin and calculates posterior densities 
posteriorcarriers_testinglargeiter <- runMCMC(CompcarriersMCMC, niter = 10000, thin=1, nburnin=0, 
                     summary = TRUE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = carriersMCMC_Inits) 

#creates list of values for both chains in order to plot easier
carriers_combinedchains_testinglargeiter <- mcmc.list(posteriorcarriers_testinglargeiter$samples$chain1, posteriorcarriers_testinglargeiter$samples$chain2)

#plots combined chains
plot(carriers_combinedchains_testinglargeiter[,c('beta0','betaWh','betaHi')])

#plots autocerrelation plots
autocorr.plot(posteriorcarriers_testinglargeiter$samples$chain1)
autocorr.plot(posteriorcarriers_testinglargeiter$samples$chain2)


#runs MCMC with 2 chains, small number of iterations and no burnin and calculates posterior densities 
posteriorcarriers_testingsmalliter <- runMCMC(CompcarriersMCMC, niter = 2000, thin=1, nburnin=0, 
                                              summary = TRUE, samples = TRUE, nchains=2, 
                                              samplesAsCodaMCMC=TRUE, inits = carriersMCMC_Inits) 

#creates list of values for both chains in order to plot easier
carriers_combinedchains_testingsmalliter <- mcmc.list(posteriorcarriers_testingsmalliter$samples$chain1, posteriorcarriers_testingsmalliter$samples$chain2)

#plots combined chains
plot(carriers_combinedchains_testingsmalliter[,c('beta0','betaWh','betaHi')])




####Finalised 2 Covariate Model Section ####

#runs MCMC with 2 chains and calculates posterior densities 
posteriorcarriers_bothcovariates <- runMCMC(CompcarriersMCMC, niter = 50000, thin=1, nburnin=500, 
                                              summary = TRUE, samples = TRUE, nchains=2, 
                                              samplesAsCodaMCMC=TRUE, inits = carriersMCMC_Inits) 

#creates list of values for both chains in order to plot easier
carriers_combinedchains_bothcovariates <- mcmc.list(posteriorcarriers_bothcovariates$samples$chain1, posteriorcarriers_bothcovariates$samples$chain2)

#plots combined chains
plot(carriers_combinedchains_bothcovariates[,c('beta0','betaWh','betaHi')])

#calculates the WAIC for the model
calculateWAIC(CompcarriersMCMC)

#calculates 95% credible intervals for each beta parameter
beta0quantities <- c(as.numeric(carriers_combinedchains_bothcovariates[[1]][,1]),as.numeric(carriers_combinedchains_bothcovariates[[2]][,1]))
betaWhquantities <- c(as.numeric(carriers_combinedchains_bothcovariates[[1]][,2]),as.numeric(carriers_combinedchains_bothcovariates[[2]][,2]))
betaHiquantities <- c(as.numeric(carriers_combinedchains_bothcovariates[[1]][,3]),as.numeric(carriers_combinedchains_bothcovariates[[2]][,3]))
beta0credinterval <- quantile(beta0quantities,probs=c(0.025,0.925))
betaWhcredinterval <- quantile(betaWhquantities,probs=c(0.025,0.925))
betaHicredinterval <- quantile(betaHiquantities,probs=c(0.025,0.925))

####sampling Pi and PPv section####

pisamples <- c()
for (chain in 1:2){
  for (i in 1:length(carriers_combinedchains_bothcovariates[[chain]][,1])){
    beta0 <- as.numeric(carriers_combinedchains_bothcovariates[[chain]][i,1])
    betaHi <- as.numeric(carriers_combinedchains_bothcovariates[[chain]][i,2])
    betaWh <- as.numeric(carriers_combinedchains_bothcovariates[[chain]][i,3])
    X_hp <- sum((dataQ5$X_hp*dataQ5$ni))/sum(dataQ5$ni)
    X_hi <- sum((dataQ5$X_hi*dataQ5$ni))/sum(dataQ5$ni)
    normalisedXhi <- ((X_hp - mean(dataQ5$X_hp[1:15]))/sd(dataQ5$X_hp[1:15]))
    normalisedXhp <- ((X_hi - mean(dataQ5$X_hi[1:15]))/sd(dataQ5$X_hi[1:15]))
    logitpi <- beta0 + betaHi*normalisedXhi+betaWh*normalisedXhp
    pi <- exp(logitpi)/(1+exp(logitpi))
    pisamples <- c(pisamples,pi)
  }
}


hist(pisamples)
mean(pisamples)

PPVsamplesQ5 <- c()
for (i in 1:length(pisamples)){
    #takes samples from posteriors
    Se <- rbeta(1,21,19)
    Sp <- rbeta(1,82,2)
    #plugs samples into equation for PPV
    PPVQ5 <- (Se*pisamples[i])/(Se*pisamples[i]+((1-pisamples[i])*(1-Sp)))
    PPVsamplesQ5 <- c(PPVsamplesQ5,PPVQ5) 
}

hist(PPVsamplesQ5)
mean(PPVsamplesQ5)

####Sensitivity of Beta Priors Models####

##Larger Standard Deviation
# Specifies model containing both covariates with different beta posterior
carriermodelcodelargesd <- nimbleCode({
  #Specifying the likelihood 
  for (i in 1:N){
    ncarriers[i] ~ dbinom(pi[i],ni[i])
  }
  for (i in 1:N){
    logit(pi[i]) <- beta0 + betaWh*((X_hp[i] - mean(X_hp[1:N]))/sd(X_hp[1:N])) + betaHi*((X_hi[i] - mean(X_hi[1:N]))/sd(X_hi[1:N]))
  }
  #Priors
  beta0 ~ dnorm(0, 0.000625)
  betaWh ~ dnorm(0, 0.000625)
  betaHi ~ dnorm(0, 0.000625)
})



#builds the model
carrierslargesd <- nimbleModel(code = carriermodelcodelargesd, name = 'carriers', constants = carriermodel_constants, data = carriermodel_data, inits = initvalues) 

#compiles model
Compcarrierslargesd <- compileNimble(carrierslargesd)

# sets up the quantities to be monitored. 
carrierslargesdconfig <- configureMCMC(carrierslargesd, enableWAIC = FALSE, monitors = c('beta0','betaWh','betaHi'), print = TRUE)  

# builds the MCMC algorithm
carrierslargesdMCMC <- buildMCMC(carrierslargesdconfig)

# compile the MCMC chain 
CompcarrierslargesdMCMC <- compileNimble(carrierslargesdMCMC, project = carrierslargesd)

#Sets initial values for MCMC chains
carrierslargesdMCMC_Inits <- list(list(beta0=-20, betaWh=-20, betaHi=-20), 
                           list(beta0=20, betaWh=20, betaHi=20))

#runs MCMC with 2 chains and calculates posterior densities 
posteriorcarrierslargesd_bothcovariates <- runMCMC(CompcarrierslargesdMCMC, niter = 50000, thin=1, nburnin=500, 
                                            summary = TRUE, samples = TRUE, nchains=2, 
                                            samplesAsCodaMCMC=TRUE, inits = carrierslargesdMCMC_Inits) 

#creates list of values for both chains in order to plot easier
carrierslargesd_combinedchains_bothcovariates <- mcmc.list(posteriorcarrierslargesd_bothcovariates$samples$chain1, posteriorcarrierslargesd_bothcovariates$samples$chain2)

#plots combined chains
plot(carrierslargesd_combinedchains_bothcovariates[,c('beta0','betaWh','betaHi')])


##Larger Mean
# Specifies model containing both covariates with different beta posterior
carriermodelcodelargemean <- nimbleCode({
  #Specifying the likelihood 
  for (i in 1:N){
    ncarriers[i] ~ dbinom(pi[i],ni[i])
  }
  for (i in 1:N){
    logit(pi[i]) <- beta0 + betaWh*((X_hp[i] - mean(X_hp[1:N]))/sd(X_hp[1:N])) + betaHi*((X_hi[i] - mean(X_hi[1:N]))/sd(X_hi[1:N]))
  }
  #Priors
  beta0 ~ dnorm(6, 0.04)
  betaWh ~ dnorm(6, 0.04)
  betaHi ~ dnorm(6, 0.04)
})



#builds the model
carrierslargemean <- nimbleModel(code = carriermodelcodelargemean, name = 'carriers', constants = carriermodel_constants, data = carriermodel_data, inits = initvalues) 

#compiles model
Compcarrierslargemean <- compileNimble(carrierslargemean)

# sets up the quantities to be monitored. 
carrierslargemeanconfig <- configureMCMC(carrierslargemean, enableWAIC = FALSE, monitors = c('beta0','betaWh','betaHi'), print = TRUE)  

# builds the MCMC algorithm
carrierslargemeanMCMC <- buildMCMC(carrierslargemeanconfig)

# compile the MCMC chain 
CompcarrierslargemeanMCMC <- compileNimble(carrierslargemeanMCMC, project = carrierslargemean)

#Sets initial values for MCMC chains
carrierslargemeanMCMC_Inits <- list(list(beta0=1, betaWh=1, betaHi=1), 
                                  list(beta0=9, betaWh=9, betaHi=9))

#runs MCMC with 2 chains and calculates posterior densities 
posteriorcarrierslargemean_bothcovariates <- runMCMC(CompcarrierslargemeanMCMC, niter = 50000, thin=1, nburnin=500, 
                                                   summary = TRUE, samples = TRUE, nchains=2, 
                                                   samplesAsCodaMCMC=TRUE, inits = carrierslargemeanMCMC_Inits) 

#creates list of values for both chains in order to plot easier
carrierslargemean_combinedchains_bothcovariates <- mcmc.list(posteriorcarrierslargemean_bothcovariates$samples$chain1, posteriorcarrierslargemean_bothcovariates$samples$chain2)

#plots combined chains
plot(carrierslargemean_combinedchains_bothcovariates[,c('beta0','betaWh','betaHi')])


####Model With One Covariate WAIC comparisons####


# Specifies model containing only betaWh covariates

carrierbetaWhmodelcode <- nimbleCode({
  #Specifying the likelihood 
  for (i in 1:N){
    ncarriers[i] ~ dbinom(pi[i],ni[i])
  }
  for (i in 1:N){
    logit(pi[i]) <- beta0 + betaWh*((X_hp[i] - mean(X_hp[1:N]))/sd(X_hp[1:N])) 
  }
  #Priors
  beta0 ~ dnorm(0, 0.04)
  betaWh ~ dnorm(0, 0.04)
})

#model constants
carrierbetaWhmodel_constants <- list(N = 15, ni = dataQ5$ni, X_hp = dataQ5$X_hp)

#data to be used
carrierbetaWhmodel_data <- list(ncarriers = dataQ5$ncarriers,pi = dataQ5$pi)

#set initial values of model
initvaluesbetaWh <- list(beta0 = 0, betaWh = 0)

#builds the model
carriersbetaWh <- nimbleModel(code = carrierbetaWhmodelcode, name = 'carriersbetaWh', constants = carrierbetaWhmodel_constants, data = carrierbetaWhmodel_data, inits = initvaluesbetaWh) 

#compiles model
CompcarriersbetaWh <- compileNimble(carriersbetaWh)

# sets up the quantities tp be monitored. 
carriersbetaWhconfig <- configureMCMC(carriersbetaWh, enableWAIC = TRUE, monitors = c('beta0','betaWh'), print = TRUE)  

# builds the MCMC algorithm
carriersbetaWhMCMC <- buildMCMC(carriersbetaWhconfig)

# compile the MCMC chain 
CompcarriersbetaWhMCMC <- compileNimble(carriersbetaWhMCMC, project = carriersbetaWh)

#sets initial values for chains
carriersbetaWhMCMC_Inits <- list(list(beta0=-5, betaWh=-5), 
                                 list(beta0=5, betaWh=5))
#run MCMC
posteriorbetaWh <- runMCMC(CompcarriersbetaWhMCMC, niter = 50000, thin=1, nburnin=500, 
                           summary = TRUE, samples = TRUE, nchains=2, 
                           samplesAsCodaMCMC=TRUE, inits = carriersbetaWhMCMC_Inits) 

#calculates WAIC
calculateWAIC(CompcarriersbetaWhMCMC)





# Specifies model containing only betaHi covariates

carrierbetaHimodelcode <- nimbleCode({
  #Specifying the likelihood 
  for (i in 1:N){
    ncarriers[i] ~ dbinom(pi[i],ni[i])
  }
  for (i in 1:N){
    logit(pi[i]) <- beta0 + betaHi*((X_hi[i] - mean(X_hi[1:N]))/sd(X_hi[1:N]))
  }
  #Priors
  beta0 ~ dnorm(0, 0.04)
  betaHi ~ dnorm(0, 0.04)
})

#model constants
carrierbetaHimodel_constants <- list(N = 15, ni = dataQ5$ni, X_hi = dataQ5$X_hi)

#data to be used
carrierbetaHimodel_data <- list(ncarriers = dataQ5$ncarriers,pi = dataQ5$pi)

#set initial values of model
initvaluesbetaHi <- list(beta0 = 0, betaHi = 0)

#builds the model
carriersbetaHi <- nimbleModel(code = carrierbetaHimodelcode, name = 'carriersbetaHi', constants = carrierbetaHimodel_constants, data = carrierbetaHimodel_data, inits = initvaluesbetaHi) 

#compiles model
CompcarriersbetaHi <- compileNimble(carriersbetaHi)

# sets up the quantities tp be monitored. 
carriersbetaHiconfig <- configureMCMC(carriersbetaHi, enableWAIC = TRUE, monitors = c('beta0','betaHi'), print = TRUE)  

# builds the MCMC algorithm
carriersbetaHiMCMC <- buildMCMC(carriersbetaHiconfig)

# compile the MCMC chain 
CompcarriersbetaHiMCMC <- compileNimble(carriersbetaHiMCMC, project = carriersbetaHi)

#sets initial values for chains
carriersbetaHiMCMC_Inits <- list(list(beta0=-5, betaHi=-5), 
                           list(beta0=5, betaHi=5))
#run MCMC
posteriorbetaHi <- runMCMC(CompcarriersbetaHiMCMC, niter = 50000, thin=1, nburnin=500, 
                     summary = TRUE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = carriersbetaHiMCMC_Inits) 


calculateWAIC(CompcarriersbetaHiMCMC)

#Question 6

posteriorpredictivesamples

meanpredictivesamples <- mean(posterior)
















####Question 6####

#gives random list of which chain and which iteration to use for each of the 10000 samples 
randomchain <- sample(1:2,10000,replace=TRUE)
randomiter <- sample(1:49500,10000,replace=FALSE)

nocarrierssamples <- list()
for (j in 1:10000){
  beta0Q6 <- as.numeric(carriers_combinedchains_bothcovariates[[randomchain[j]]][randomiter[j],1])
  betaHiQ6 <- as.numeric(carriers_combinedchains_bothcovariates[[randomchain[j]]][randomiter[j],2])
  betaWhQ6 <- as.numeric(carriers_combinedchains_bothcovariates[[randomchain[j]]][randomiter[j],3])
  numbercarrierssamplesiter <- c() 
  for (i in 1:15){
    normXhiQ6 <- (dataQ5$X_hi[i]-mean(dataQ5$X_hi[1:15]))/sd(dataQ5$X_hi[1:15])
    normXhpQ6 <- (dataQ5$X_hp[i]-mean(dataQ5$X_hp[1:15]))/sd(dataQ5$X_hp[1:15])
    logitpiQ6 <- beta0Q6+betaWhQ6*normXhpQ6+betaHiQ6*normXhiQ6
    piQ6 <- exp(logitpiQ6)/(1+exp(logitpiQ6))
    numbercarrierssamplesiter <- c(numbercarrierssamplesiter,rbinom(1,dataQ5$ni[i],piQ6))
  }
  nocarrierssamples <- append(nocarrierssamples,list(numbercarrierssamplesiter))
} 
nocarrierssamples

#finds and stores mean,median and max of samples
meansamples <- c()
mediansamples <- c()
maxsamples <- c()
for (i in 1:10000){
  meansamples <- c(meansamples,mean(nocarrierssamples[[i]]))
  mediansamples <- c(mediansamples,median(nocarrierssamples[[i]]))
  maxsamples <- c(maxsamples,max(nocarrierssamples[[i]]))
}

hist(meansamples)
hist(mediansamples)
hist(maxsamples)

originalmean <- mean(dataQ5$ncarriers)
originalmedian <- median(dataQ5$ncarriers)
originalmax <- max(dataQ5$ncarriers)

pvaluemeanleft <- length(meansamples[which(meansamples<originalmean)])/length(meansamples)
pvaluemeanright <- length(meansamples[which(meansamples>originalmean)])/length(meansamples)
pvaluemedianleft <- length(mediansamples[which(mediansamples<originalmedian)])/length(mediansamples)
pvaluemedianright <- length(mediansamples[which(mediansamples>originalmedian)])/length(mediansamples)
pvaluemaxleft <- length(maxsamples[which(maxsamples<originalmax)])/length(maxsamples)
pvaluemaxright <- length(maxsamples[which(maxsamples>originalmax)])/length(maxsamples)

pvaluemeanleft
pvaluemeanright
pvaluemedianleft
pvaluemedianright
pvaluemaxleft
pvaluemaxright


