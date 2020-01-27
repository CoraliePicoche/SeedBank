#############
#23/01/2020 CP: Uses quadratic programming to get closer approximations of the interaction matrix and growth rates of species in presence, heavily relying on the code from Maynard et al. 2019 (https://github.com/dsmaynard/reconciling_coexistence/blob/master/code/functions_algae.R)
#We want solve AN*+r_2=0 where A is the BH interaction matrix, N* is the species abundance at equilibrium, and r_2=1-exp(r_mean) with r_mean the mean growth rates of all species, under the conditions N*>0, i.e. -(A)^-1*r_2>0 (feasible equilibrium) 
#Should we also have stg like and b_{ij} \approx -a_{ij}N_j/(1+\sum_l a_{il}N_l), knowing that it is already a starting point
#Maynard et al. 2019 used the conditions that growth_rate>0 and self-regulation>0 (a_ii>0).
#############

#rm(list=ls())
graphics.off()
set.seed(42)

#library(quadprog)
#library(lsei)
library(limSolve)
library(Matrix)
source("script/step_functions.r")
source("script/matrix_MAR_clean.r")

nspp=11
r <- runif(nspp)
A <- -matrix(runif(nspp^2), nspp,nspp)
diag(A) <- diag(A)*2
x_obs <- runif(nspp)

tol=0.1 #In the first example: 1000
sp=1:11

##Random initial values for growth rate #With this, we can have AN-(exp(r)-1)=0 after quadratic optim
r_mean=runif(sp,0.5,1.5)
#SV values for growth_rates #With this, we cannot have AN-(exp(r)-1)=0 after quadratic optim
#Niche area to compute growth rates + range of optimal temperatures

if(1==0){
A=10^(3.1)/365.25
T_min=288
T_max=298

T_opt=runif(length(sp),293,298) #Between 20 And 25 for most species

#b parameter to keep the same niche area
f_to_optimize_B=function(b,T_min,T_max,T_opt,A){
        f1=integrate(Vectorize(growth_rate),lower=T_min-5,upper=T_max+5,T_opt,b)$val
        tmp=abs(f1-A)
        return(tmp)
}

B=rep(NA,length(sp))
for(i in 1:length(sp)){
        B[i]=optimize(f_to_optimize_B,T_min,T_max,T_opt[i],A,interval=c(0,100))$minimum
}
r_mean=growth_rate(293,T_opt,B)
}

###Based on Bissinger expression
#With mean temperature
tab_hydro=read.table("param/Augerhydro.txt",header=T,sep=";")
mean_temp=mean(tab_hydro[,"TEMP"],na.rm=T)

r_mean=rep(growth_rate_Bissinger(mean_temp,0.5),nspp)

r_2=1-exp(r_mean)


##Interaction matrix
####RANDOM
inter_mat=matrix(rnorm(length(sp)^2,0,0.05),length(sp),length(sp)) 
diag(inter_mat)=-abs(diag(inter_mat))*2

####EXACT
load("param/Auger_pencen_null_regular_common_MO.RData")
name_spp=colnames(cis$call$model$B)
inter_mat=clean_matrix(cis,signif=F)
rownames(inter_mat)=colnames(inter_mat)=name_spp

f_exact_resolution=function(B,N){
        A=matrix(NA,nrow=nrow(B),ncol=ncol(B))
        for(i in 1:dim(B)[1]){
                sumBJ=sum(B[i,]*N)
                sumB=sum(B[i,])
                for(j in 1:dim(B)[2]){
                        A[i,j]=-1/N[j]*B[i,j]/(1+sumB)
                }
        }

        return(A)
}

#### RANDOM
aN=c(10^6,10^6,runif(length(sp)-2,10^2,10^4))
x_obs=aN

###EXACT
pop_table=read.table("param/abundances_Auger.txt")
aN=x_obs=pop_table[name_spp,]

A_coast=f_exact_resolution(inter_mat,x_obs)
Avec=as.numeric(t(A_coast))

parvec <- c(Avec,r_2)
# number of pars to fit
npar <- length(parvec)
# get the bounds

Alow <- Avec-abs(Avec*tol) 
Aupp <- Avec+abs(Avec*tol) 

#rlow <- rep(0,length(r_2))
rlow <- r_2-abs(r_2*tol) #In Maynard's code
#rupp <- r_2+abs(r_2*tol) #In Maynard's code
rupp=rep(0,length(r_2))

# constrain the diagonals to be negative, and some small value away from zero
#Aupp[as.logical(as.numeric(diag(nspp)))] <- min(-1e-4,max(diag(A_coast)))/4 #This was for LV
Alow[as.logical(as.numeric(diag(nspp)))] = rep(0,nspp)
# inequality constraints
h <- c(Alow,rlow,-Aupp,-rupp)
G <- rbind(diag(npar),-diag(npar))
# equality constraints, assuming at equilibirum
E <- cbind(as.matrix(bdiag(replicate(nspp,matrix(x_obs,nrow=1),simplify = F))),diag(nspp))
f <- rep(0,nrow(E))
# fit the qp model, returning a silent warning if it doesn't converge

#fit <- tryCatch(lsei(A=diag(npar),B=parvec,E=E,F=f,G=G,H=h),  error=function(e) e, warning=function(w) w)
#fit <-lsei(A=diag(npar),B=parvec,E=E,F=f,G=G,H=h)
fit=limSolve::lsei(A = diag(npar), B = parvec, E = E, F = f, G = G, H = h,type=1,verbose=TRUE,fulloutput=TRUE)

A_after=t(matrix(fit$X[1:nspp^2],nspp,nspp))
r_after=log(1-fit$X[(nspp^2+1):length(fit$X)])


plot(1:length(parvec),100*(abs(fit$X-parvec)/parvec),xlab="Initial Values",ylab="Absolue %Difference")


list_inter=list(A_coast,5*A_coast)
