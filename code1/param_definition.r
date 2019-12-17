#########
#17/12/2019 CP : Define parameters
#########
library(optimr)
source("code1/step_functions.r")

set.seed(42)

k_coast2ocean=5
k_inter2intra=1/0.1931

e=0.64
germination=0.01 #Could be 0.1, 0.001
A=10^(3.1)
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

#Interaction matrix, coastal
inter_mat=matrix(rnorm(length(sp)^2,0,0.01),length(sp),length(sp))
tmp_intra=rep(NA,length(sp),length(sp))
tmp_mat=inter_mat
diag(tmp_mat)=rep(NA,length(sp))
for(i in 1:length(sp)){
        inter_mat[i,i]=mean(tmp_mat[,i],na.rm=T)*k_inter2intra
}

inter_ocean=inter_mat*k_coast2ocean
if(any(inter_mat<(-1))|any(inter_ocean<(-1))){
        print("Overcompensation, you may want to stop")
}

list_inter=list(inter_mat,inter_ocean)

#Sinking rate
S=rbeta(length(sp),0.55,1.25)*30

#Gamma
resuspension=0.5*S
Gamma=resuspension*germination


