#########
#17/12/2019 CP : Define parameters
#########
library(optimr)
source("code1/step_functions.r")

set.seed(42)

k_coast2ocean=1.5 #Was 5 before but led to overcompensation

#base_inter2intra=-0.382872
#k_inter2intra=0.603406 was the coefficient in self-regulation~-0.382872+0.603406*impact
k_inter2intra=-0.493729 #is the coefficient in self-regulation~-0.365043-0.493729*vulnerability
base_inter2intra=-0.365043
add_modules=TRUE

k_sediment2resuspension=10^(-1)*0.5
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
iter=iter+1

inter_mat=matrix(rnorm(length(sp)^2,0,0.05),length(sp),length(sp)) #In REPHY littoral, normal law with mean=0.01 and sd=0.05 for interspecies coefficients
tmp_intra=rep(NA,length(sp),length(sp))
tmp_mat=inter_mat
diag(tmp_mat)=rep(NA,length(sp))
for(i in 1:length(sp)){
#        inter_mat[i,i]=mean(tmp_mat[,i],na.rm=T)*k_inter2intra-0.088 #based on impact relationship with self-regulation
	inter_mat[i,i]=mean(tmp_mat[i,],na.rm=T)*k_inter2intra+base_inter2intra #based on vulnerability relationship with self-regulation
#	inter_mat[i,i]=sum(abs(tmp_mat[i,]),na.rm=T)+0.25 #based on F. Barraquand's work on BH model with mixed interactions, might be similar to Haydon's matrices in Haydon (2000) Ecology

}

#Centric 4*4, Pennate 4*4, Dino 2*2
if(add_modules){
	inter_mat[5:10,1:4]=0
	inter_mat[1:4,5:10]=0
	inter_mat[5:8,9:10]=0
	inter_mat[9:10,5:8]=0
}

inter_ocean=inter_mat*k_coast2ocean

#### Check if both interaction matrices can lead to a stable, positive equilibrium
#Compute mean growth rate
b_middle=optimize(f_to_optimize,T_min,T_max,293,A,interval=c(0,100))$minimum
r_mean=growth_rate(293,T_opt,B)

#### Check if both interaction matrices can lead to a stable, positive equilibrium
Nstar_coast=solve(-1*inter_mat)%*%(exp(r_mean)-1)
Nstar_ocean=solve(-1*inter_ocean)%*%(exp(r_mean)-1)
}
print(iter)
if(any(inter_mat<(-1))|any(inter_ocean<(-1))){
        print("Overcompensation, you may want to stop")
}


list_inter=list(inter_mat,inter_ocean)

#Sinking rate
S=rbeta(length(sp),0.55,1.25)*30/100
names(S)=sp
#Manip so that Chaetoceros spp and THA have the highest sinking rate
tmp_S=S
or=order(S,decreasing=T)
tmp_S[c("CHD","CHS","THA")]=S[or[1:3]]
tmp_S[or[1:3]]=S[1:3]
S=tmp_S

#Gamma
#The ones that sink the most are the ones that get the least suspended
or_dec=order(S,decreasing=T)
or_inc=order(S,decreasing=F)
tmp_S=S
tmp_S[or_dec]=S[or_inc]

resuspension=k_sediment2resuspension*tmp_S
Gamma=resuspension*germination
