#JAGS demo

library(rjags)
library(mcmcplots)

#model
fencemodel <- "
model{

  for(i in 1:N){
    y[i] ~ dpois(theta[i]*x[i])
    theta[i] ~ dgamma(alpha, beta)
  }
  
  #priors
  alpha ~ dgamma(2, 1)
  beta ~ dexp(1)
  
  #pop level params
  pmean_theta <- alpha/beta
  pstd_theta <- sqrt(alpha/beta^2)
}
"

#data
sf <- data.frame(y = c(138, 91, 132, 123, 173, 124, 109, 154, 138, 134),
                 x = c(72, 50, 55, 60, 78, 63, 54, 70, 80, 68))
#pass data as a list, each element must match names in model
datlist <- list(y = sf$y, 
                x = sf$x,
                N = nrow(sf))
#provide inits for params (not required but recommended)
initslist <- list(list(alpha = 0.5, beta = 1),
                  list(alpha = 4, beta = 2),
                  list(alpha = 30, beta = 1.5))

#fit model
jm <- jags.model(file = textConnection(fencemodel),
                 data = datlist,
                 inits = initslist,
                 n.chains = 3)

jm_samples <- coda.samples(model = jm,
                           variable.names = c("theta", "alpha", "beta",
                                              "pmean_theta", "pstd_theta"),
                           n.iter = 3000)
#trace plots
traplot(jm_samples, parms = c("theta","alpha","beta"))
#ACF and other plots
mcmcplot(jm_samples, parms = c("alpha","beta","pmean_theta"))
#Rhat
gelman.diag(jm_samples)
#caterpillar plots
caterplot(jm_samples, parms = "theta", reorder=F)

summary(jm_samples)[[1]] #chain 1
head(jm_samples)[[1]]
#convenience functions to get 95% CI
broom.mixed::tidyMCMC(jm_samples, conf.int = T)

#penalty is effective number of params
dic.samples(jm, n.iter=3000)


#moth example

#data
moths <- data.frame(site = rep(1:7, each = 2),
                    morph = rep(1:2, times = 7),
                    distance = rep(c(0, 7.2, 24.1, 30.2, 36.4, 41.5, 51.2), each = 2),
                    placed = rep(c(56, 80, 52, 60, 60, 84, 92), each = 2),
                    removed = c(13, 14, 28, 20, 18, 22, 9, 16, 16, 23, 20, 40, 24, 39))

mothmodel <- "
model{

  for(i in 1:N){
    y[i] ~ dbin(p[i],n[i]) #likelihood
    ypred[i] ~ dbin(p[i],n[i]) #posterior predictive
    p[i] <- ilogit(beta1+beta2*x[i]+beta3*m[i]+beta4*x[i]*m[i])
  }
  
  #priors
  beta1 ~ dnorm(0, 1/5)
  beta2 ~ dnorm(0, 1/5)
  beta3 ~ dnorm(0, 1/5)
  beta4 ~ dnorm(0, 1/5)
}
"

datlistmoth <- list(y = moths$removed, 
                    x = moths$distance,
                    n = moths$placed,
                    m = ifelse(moths$morph==1,0,1),
                    N = nrow(moths))


mothjm <- jags.model(file = textConnection(mothmodel),
                 data = datlistmoth,
                 n.chains = 3)

moth_samples <- coda.samples(model = mothjm,
                           variable.names = c("beta1","beta2","beta3","beta4"),
                           n.iter = 6000)


#trace plots
traplot(moth_samples)
#ACF and other plots
#mcmcplot(moth_samples)
#Rhat
gelman.diag(moth_samples)
#caterpillar plots
#slopes
caterplot(moth_samples, reorder=F, parms = c("beta2","beta4"))
#intercepts
caterplot(moth_samples, reorder=F, parms = c("beta1","beta3"))

broom.mixed::tidyMCMC(moth_samples, conf.int = T)

mothmodel2 <- "
model{

  for(i in 1:N){
    y[i] ~ dbin(p[i],n[i]) #likelihood
    p[i] <- ilogit(beta1[m[i]]+beta2[m[i]]*x[m[i]])
  }
  
  #priors
  for(i in 1:2) {
    beta1[i] ~ dnorm(0, 1/5)
    beta2[i] ~ dnorm(0, 1/5)
  }
}
"

datlistmoth <- list(y = moths$removed, 
                    x = moths$distance,
                    n = moths$placed,
                    m = moths$morph,
                    N = nrow(moths))


mothjm2 <- jags.model(file = textConnection(mothmodel2),
                     data = datlistmoth,
                     n.chains = 3)

moth_samples <- coda.samples(model = mothjm2,
                             variable.names = c("beta1","beta2"),
                             n.iter = 6000)
#slopes
caterplot(moth_samples, reorder=F, parms = c("beta2","beta4"))
#intercepts
caterplot(moth_samples, reorder=F, parms = c("beta1","beta3"))

