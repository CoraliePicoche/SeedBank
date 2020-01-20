################################################################################################################################################################
### Frederic Barraquand
### n-species Beverton-Holt type interaction model (mutualism possible) based on previous competition model
### also adds two environmental variables [do it now or later?]
### Competition models (potentially n-species), FB 13/01/2016. R and JAGS. 
### Simulates and fits nonlinear Beverton-Holt type of dynamics (simpler than Ricker at first) and then fits MAR models to data
### Stuff to do: large number of species with a gradient of r/K strategies (rank abundance distrib like plancton) + arcane observation models like plankton
################################################################################################################################################################


rm(list=ls())
graphics.off()

library("R2jags")      # Load R2jags package

set.seed(42) 

############## Simulation of data ###############################
n.time<-500  # Number of years or more generally temporal units // also tried 100 before, more difficult
index_time<-1:n.time
n.species<-10 # Number of species in the community
N<-matrix(0, nrow=n.time,n.species,byrow=TRUE) # Matrix of abundances or densities
epsilon<-matrix(0, nrow=n.time-1,n.species,byrow=TRUE) # Noise on growth rates
N1<-matrix(1,n.species,byrow=TRUE) # Vector of initial abundances - could be changed
N[1,]<-N1
r<-runif(n.species,0.5,1) # Maximal growth rates -- instead of (1,2) or (2,3)
#r<-runif(n.species,2,3) # Maximal growth rates
sigma<-matrix(0.2,n.species,byrow=TRUE)

### Environmental variables
#1. env variables themselves
seasonality<-2*sin(2*pi*n.time/24)	# must be enough to affect the growth rates
y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n.time,sd=sqrt(0.5) )
y2noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=n.time,sd=sqrt(1) ) # more noisy
y1<-seasonality+y1noise
y2<-seasonality+y2noise
#2. effect of environmental covariates
c1=runif(n.species,0.4,0.6)
c2=runif(n.species,-0.2,0.2)

# Interaction coefficient matrix
### alpha<-matrix(runif(n.species*n.species,-0.1,0.1), ncol=n.species) 
### I could also avoid almost zero interactions -- cf GCausality code

for(iter in 1:10){
print(iter)
alpha<-matrix(rnorm(n.species*n.species,0.0,0.05), ncol=n.species) 

for (i in 1:n.species){alpha[i,i]<-sum(abs(alpha[i,]))+0.25}
### Note: Not so different from Haydon's matrices in Haydon (2000) Ecology
### corrected to row sum for diagonal dominance
### Should I incorporate a trade off between r and alpha for coexistence? -> perhaps later


### Simulation of data

for (t in 1:(n.time-1))
	{
    epsilon[t,]<-rnorm(n.species,mean = 0, sd = sigma)# Noise
		N[t+1,]<-N[t,]*exp(r+ y1[t]*c1 + y2[t]*c2  + epsilon[t,])/pmax(0.0001,1+alpha %*% N[t,])
}

### Plotting interaction matrix
image(alpha)
image(alpha-diag(n.species)*diag(alpha))

### Equilibrium values - some algebra shows that we have 
N_star<-solve(alpha) %*% (exp(r)-1) ### why is there a negative equilibrium abundance? compute stuff for N=3?
print(N_star)
##eigen(alpha)$values ### But that's not the jacobian around the fixed point, compute that...
}

stop()

## Plotting time series
matplot(1:n.time,N,type="b")
matlines(1:n.time,N,type="b")

## Plotting time series
pdf("comp_and_mutualism/InformativePriors_andForcing/TimeSeriesLogScale.pdf")
matplot(1:n.time,log(N),type="b")
matlines(1:n.time,log(N),type="b")
dev.off()

## Plot growth rates as a function of epsilon (representing here environmental variation)
plot(epsilon,log(N[2:n.time,])-log(N[1:(n.time-1),]))

## SAD
barplot(sort(N_star,decreasing = TRUE))
pdf("comp_and_mutualism/InformativePriors_andForcing/SAD.pdf")
barplot(log(sort(N_star,decreasing = TRUE)))
dev.off()
hist(log(sort(N_star,decreasing = TRUE)))

############################################################################
### Fitting the nonlinear model 
############################################################################

# Bundle data
jags.data <- list(T=n.time,S=n.species,logN=log(N),y1=y1,y2=y2)

sink("comp_and_mutualism/InformativePriors_andForcing/BH.mat.interact.txt")
cat("
    model {
    
    # Priors and constraints
    for (i in 1:S)
      {
      logN[1,i] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    
      # Priors population dynamics
      r[i] ~ dnorm(1,0.001) # below the truth, rather flat prior
      sigma[i] ~ dunif(0.01,5) # rather vague 
      sigma2[i]<-pow(sigma[i], 2)
      tau[i]<-pow(sigma[i],-2)
      # Could be Wishart prior for variance-covariance matrix, I assumed no correlation there.
      # Priors effect of the environment
      c1[i] ~ dunif(0,1)
      c2[i] ~ dunif(-1,1)
    
      #Interaction matrix
      for (j in 1:S)
        {
          alpha[i,j]~dnorm(ifelse(i==j,0.5,0.0),10) ###Smaller mean intersp. interaction
          
        }
      }# end of first loop on i
    
      # Likelihood
      # state process
    
      for (t in 1:(T-1))
        {    
        for (i in 1:S)
          {
          logN[t+1,i] ~ dnorm(logNupdate[t,i],tau[i])
          logNupdate[t,i] <- logN[t,i] + r[i] - log(F[t,i]) + c1[i] * y1[t] +  c2[i] * y2[t]
          Nvec[t,i]<-exp(logN[t,i])          
          F[t,i]<- max(0.00001,1 + alpha[i,] %*% Nvec[t,])
          #We need the max correction because 1 + alpha[i,] %*% Nvec[t,] can get negative for certain values of the chain
        }#end of loop on i
      }# end of loop on t
  
    }
    ",fill=TRUE)
sink()

#doesnt work in JAGS
#F[t,i]<-1
#for (j in 1:S)
#  {
#      F[t,i]<-F[t,i] + alpha[i,j]*exp(logN[t,j])
#  }


# Initial values
#inits <- function () {  list(sigma=runif(n.species,0.1,2),r=runif(n.species,2,2.5), alpha=matrix(runif(n.species*n.species,0,1),n.species,n.species))}
## Doesn't work - perhaps due to runif() for alpha values

inits <- function () {
  list(sigma=runif(n.species,0.1,2),r=runif(n.species,0.6,0.9), c1=runif(n.species,0.5,1),c2=runif(n.species,-1,1), alpha=matrix(rnorm(n.species*n.species,0,0.05),n.species,n.species))}
## doesn't work -- perhaps due to rnorm() being too close to zero for diag values

inits <- function () {
  list(sigma=runif(n.species,0.1,2),r=runif(n.species,0.6,0.9), c1=runif(n.species,0.5,1),c2=runif(n.species,-1,1),alpha=matrix(rnorm(n.species*n.species,0,0.025)+0.5*diag(n.species),n.species,n.species))}
#does not change anything

# Parameters monitored
parameters<-c("r","alpha","sigma","c1","c2")

# MCMC settings
nc <- 3 #number of chains
nb <- 4000 #14000 # “burn in”
#ni <- 14000# “number of iterations” # that's for a symmetric distrib...
ni<-10000 #34000
nt <- 10 # “thinning”

# run model
out <- jags(jags.data, inits, parameters, "comp_and_mutualism/InformativePriors_andForcing/BH.mat.interact.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out, dig = 2)
x=as.mcmc(out$BUGSoutput$summary)
write.table(x,file="comp_and_mutualism/InformativePriors_andForcing/out_BH.txt",sep="\t")

# library(mcmcplots)
# traceplots
# traplot(out,parameters)
# posterior densities
# denplot(out,parameters)

### Output results

alpha_est=out$BUGSoutput$mean$alpha
pdf("comp_and_mutualism/InformativePriors_andForcing/EstimatedAlpha_asFunctionOf_simulatedAlpha.pdf")
plot(as.vector(alpha),as.vector(alpha_est),xlab="Simulated alpha",ylab="Estimated alpha")
abline(lm( as.vector(alpha_est) ~ as.vector(alpha) ) )
abline(a=0,b=1,col="blue")
dev.off()

#############################################################################################################################
### Fitting a MAR model - using the nonlinear model as simulation and comparing to the linear approximation around eq.
#############################################################################################################################

### First, computation of the Jacobian of the log-linearized model, that is, the expected coefficients of the B matrix 
### Diagonal terms 1-\frac{\alpha_{11}N_{1}}{1+\alpha_{11}N_{1}+\alpha_{12}N_{2}}
### Off-diagonal terms -\frac{\alpha_{12}N_{2}}{1+\alpha_{11}N_{1}+\alpha_{12}N_{2}}
### Generalizes to N species using the N_i^*
### Diagonal J_ij = 1 - alpha_ii N_i/ (1+ sum_k alpha_ik N_k)
### Off diagonal J_ij = - alpha_ij N_j/ (1+ sum_k alpha_ik N_k)
########################################################################################

J=matrix(0,n.species,n.species)
for (i in 1:n.species){
  SumCompet = alpha[i,] %*% N_star
  for (j in 1:n.species){
      if (i==j)
      {
        J[i,j] = 1 - alpha[i,j]*N_star[j] /(1+  SumCompet)
      }
    else {
      J[i,j] = - alpha[i,j] *N_star[j]  /(1+  SumCompet)
    }
  }
}
J #this is essentially what we should find
image(J)
image(J-diag(n.species))

###################################
### Fitting the MAR model now #####
###################################

# Bundle data
logNmean=log(N) - colMeans(log(N))
jags.data <- list(T=n.time,S=n.species,V=2,logN=logNmean,y=cbind(y1,y2))

### New version without the algebra
sink("comp_and_mutualism/InformativePriors_andForcing/MAR.compet.txt")
cat("
    
    var logNupdate[T,S]
    
    model {
    
    # Priors and constraints
    for (i in 1:S)
    {
    logN[1,i] ~ dnorm(0,0.01) # Prior on initial pop size on the log scale
    
    # Priors population dynamics
    a[i] ~ dnorm(1,0.001) #  rather flat prior for the initial log pop growth rate
    sigma[i] ~ dunif(0.01,5) # rather vague 
    sigma2[i]<-pow(sigma[i], 2)
    tau[i]<-pow(sigma[i],-2)
    # Could be Wishart prior for variance-covariance matrix, I assumed no correlation there.
    
    # B matrix priors and 
    for (j in 1:S)
    {
    B[i,j]~dnorm(0,1)
    }
    # C matrix priors
    for (j in 1:V)
    {
    C[i,j]~dnorm(0,1)
    }
    
    }# end of first loop on i
    
    
    
    # Likelihood
    # state process
    
    for (t in 1:(T-1))
    {
    for (i in 1:S)
    {
    logN[t+1,i] ~ dnorm(logNupdate[t,i],tau[i])
    logNupdate[t,i] <- a[i] + B[i,] %*% logN[t,] + C[i,] %*% y[t,]
    }#end of loop on i
    
    }# end of loop on t
    
    
    }
    ",fill=TRUE)
sink()


# Initial values
inits <- function () {
  list(sigma=runif(n.species,0.1,2),a=rnorm(n.species,0,0.1), B=matrix(rnorm(n.species*n.species,0,1),n.species,n.species), C=matrix(rnorm(2*n.species,0,1),n.species,2))}


# Parameters monitored
parameters<-c("a","B","C","sigma")

# MCMC settings
nc <- 3 #number of chains
nb <- 4000 # “burn in”
#ni <- 34000# “number of iterations”
ni<-14000
nt <- 10 # “thinning”

# run model
out2 <- jags(jags.data, inits, parameters, "comp_and_mutualism/InformativePriors_andForcing/MAR.compet.txt", n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb, working.directory = getwd())
print(out2, dig = 2)
x=as.mcmc(out2$BUGSoutput$summary)
write.table(x,file="comp_and_mutualism/InformativePriors_andForcing/out_MAR.txt",sep="\t")

### Output results

B_est=out2$BUGSoutput$mean$B

Id = diag(n.species)
pdf("comp_and_mutualism/InformativePriors_andForcing/EstimatedB.pdf")
par(mfrow=c(1,3),pty='s')
plot(as.vector(J),as.vector(B_est),xlab="True J",ylab="Estimated B")
points(as.vector(diag(J)),as.vector(diag(B_est)),col="red",pch=21,bg="red")
abline(lm( as.vector(B_est) ~ as.vector(J) ) )
plot(as.vector(alpha_est),as.vector(-B_est+Id),xlab="Estimated alpha",ylab="Estimated -B + Id")
points(as.vector(diag(alpha_est)),as.vector(-diag(B_est - Id )),xlab="Estimated alpha",ylab="Estimated -B + Id ",col="red",pch=21,bg="red")
plot(as.vector(alpha),as.vector(-B_est + Id ),xlab="Simulated alpha",ylab="Estimated -B + Id ")
points(as.vector(diag(alpha)),as.vector(-diag(B_est - Id )),xlab="Simulated alpha",ylab="Estimated -B + Id ",col="red",pch=21,bg="red")
#abline(lm( as.vector(B_est) ~ as.vector(alpha_est) ) )
dev.off()
