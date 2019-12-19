#########
#17/12/2019 CP : Define parameters
#########
library(optimr)
source("code1/step_functions.r")

set.seed(42)

k_coast2ocean=1 #Was 5 before but crashed to negative abundances. Crashes even for 1

#k_inter2intra=1/0.1931 was the coefficient in impact~0.088+0.1931*self-regulation
k_inter2intra=1/(-0.1228) #is the coefficient in impact~-0.0296-0.1228*self-regulation
base_inter2intra=-0.0296

e=0.64
germination=0.01 #Could be 0.1, 0.001
A=10^(3.1)/365.25
T_min=288
T_max=298

####Parameter definition
#Morta
tab_morta=read.table("param_fake.csv",sep=";",header=T,dec=".")
M=tab_morta[,"mortality_rate"]
sp=tab_morta[,"Code"]

#Optimum temperature
T_opt=runif(length(sp),293,298) #Between 20 And 25 for most species
T_opt[sp=="SKE"|sp=="NAV"]=runif(2,288,293)

#b parameter
f_to_optimize=function(b,T_min,T_max,T_opt,A){
	f1=integrate(Vectorize(growth_rate),lower=T_min-5,upper=T_max+5,T_opt,b)$val
	tmp=abs(f1-A)
	return(tmp)
}

B=rep(NA,length(sp))
for(i in 1:length(sp)){
	B[i]=optimize(f_to_optimize,T_min,T_max,T_opt[i],A,interval=c(0,100))$minimum
}


Nstar_coast=matrix(-1,length(sp),length(sp))
Nstar_ocean=matrix(-1,length(sp),length(sp))
iter=0
while(sum(Nstar_coast<0)>0|sum(Nstar_ocean<0)>0){
#Interaction matrix, coastal
print(iter)
iter=iter+1

inter_mat=matrix(rnorm(length(sp)^2,0,0.01),length(sp),length(sp))
tmp_intra=rep(NA,length(sp),length(sp))
tmp_mat=inter_mat
diag(tmp_mat)=rep(NA,length(sp))
for(i in 1:length(sp)){
#        inter_mat[i,i]=mean(tmp_mat[,i],na.rm=T)*k_inter2intra-0.088 #based on impac
	inter_mat[i,i]=mean(tmp_mat[i,],na.rm=T)*k_inter2intra-base_inter2intra #based on vulnerability
#	inter_mat[i,i]=sum(abs(tmp_mat[i,]),na.rm=T)+0.25
}

inter_ocean=inter_mat*k_coast2ocean
if(any(inter_mat<(-1))|any(inter_ocean<(-1))){
        print("Overcompensation, you may want to stop")
}

#### Check if both interaction matrices can lead to a stable, positive equilibrium
#Compute mean growth rate
b_middle=optimize(f_to_optimize,T_min,T_max,293,A,interval=c(0,100))$minimum
r_mean=growth_rate(293,T_opt,B)
#r_mean=r<-runif(length(sp),0.5,1)

Nstar_coast=solve(inter_mat)%*%(exp(r_mean)-1)
Nstar_ocean=solve(inter_ocean)%*%(exp(r_mean)-1)

}


list_inter=list(inter_mat,inter_ocean)

#Sinking rate
S=rbeta(length(sp),0.55,1.25)*30

#Gamma
resuspension=0.5*S
Gamma=resuspension*germination


