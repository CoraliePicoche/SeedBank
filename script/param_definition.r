#########
#17/12/2019 CP : Define parameters
#14/01/2020 CP: Corrected the formulation of the interaction matrix to go from MAR to BH
#########
rm(list=ls())
graphics.off()

library(optimr)
source("script/step_functions.r")

set.seed(42)

###################  Fixed parameters  ###################

#Conversion from coastal interactions to oceanic interactions, assuming competition is higher in the ocean.
k_coast2ocean=1.5# Was 5 before but led to overcompensation

#Conversion from interspecific interactions to intraspecific regulation, based on data in Picoche & Barraquand (2019)
#base_inter2intra=-0.382872
#k_inter2intra=0.603406 was the coefficient in self-regulation~-0.382872+0.603406*impact
k_inter2intra=-0.493729 #is the coefficient in self-regulation~-0.365043-0.493729*vulnerability
base_inter2intra=-0.365043

#Other information on the interaction matrix
add_modules=FALSE #Centric diatoms / Pennate diatoms / Dinoflagellates can only interact within their group
intra_only=FALSE #Only intraspecific interactions (mostly to debug) 

#Conversion from sinking/sedimentation rate to resuspension rate (might be based on the volume/surface/weight or shape of the cells
#k_sediment2resuspension=10^(-1)*0.5
k_sediment2resuspension=0.5

#Exchange rate between coastal and oceanic (corresponds to tide)
e=0.64

#Germination rate of seeds
germination=0.001 #Could be 0.1, 0.001

#Niche area to compute growth rates + range of optimal temperatures
A=10^(3.1)/365.25
T_min=288
T_max=298

##############  Species-specific parameters  ##############
#Seed mortality rate
tab_morta=read.table("param/param_fake.csv",sep=";",header=T,dec=".")
M=tab_morta[,"mortality_rate"]
sp=tab_morta[,"Code"]
names(M)=sp

####Definition of growth rates
#Optimum temperature
T_opt=runif(length(sp),293,298) #Between 20 And 25 for most species
T_opt[sp=="SKE"|sp=="NAV"]=runif(2,288,293)
names(T_opt)=sp

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
names(B)=sp


####Interaction matrices
f_to_optimize_A=function(A,B,N){
#compute the log-scale Jacobian of A (more precisely J-I)
        J=matrix(NA,nrow=nrow(B),ncol=ncol(B))
        A=matrix(A,nrow=nrow(B),ncol=ncol(B))
        for(i in 1:dim(B)[1]){
                for(j in i:dim(B)[2]){
                        J[i,j]=-A[i,j]*N[j]/(1+A[i,]%*%N)
                        J[j,i]=-A[j,i]*N[i]/(1+A[j,]%*%N)
                }
        }
        tmp=sum(abs(B-J))
        return(tmp)
}

#Observed equilibrium
#aN=c(10^6,10^6,runif(length(sp)-2,10^2,10^4))
aN=c(runif(length(sp),10^3,5*10^3))
#aN=c(runif(length(sp),1,10))

#Matrix at equilibrium
Nstar_coast=matrix(-1,length(sp),length(sp))
Nstar_ocean=matrix(-1,length(sp),length(sp))
#Looking for a feasible equilibrium
iter=0
conv_coast=1
conv_ocean=1
#while(sum(Nstar_coast<0)>0|sum(Nstar_ocean<0)>0|((conv_coast+conv_ocean)!=0)){ ##We will have to go back to that...
while((conv_coast+conv_ocean)!=0){
#Interaction matrix, coastal
iter=iter+1
print(iter)
inter_mat=matrix(rnorm(length(sp)^2,0,0.05),length(sp),length(sp)) #In REPHY littoral, normal law with mean=0.01 and sd=0.05 for interspecies coefficients
colnames(inter_mat)=rownames(inter_mat)=sp

#Centric 4*4, Pennate 4*4, Dino 2*2
if(add_modules){
	inter_mat[5:10,1:4]=0
	inter_mat[1:4,5:10]=0
	inter_mat[5:8,9:10]=0
	inter_mat[9:10,5:8]=0
}

tmp_mat=inter_mat
tmp_mat[tmp_mat==0]=NA
diag(tmp_mat)=NA
for(i in 1:length(sp)){
	inter_mat[i,i]=mean(tmp_mat[i,],na.rm=T)*k_inter2intra+base_inter2intra #based on vulnerability relationship with self-regulation
#	inter_mat[i,i]=sum(abs(tmp_mat[i,]),na.rm=T)+0.25 #based on F. Barraquand's work on BH model with mixed interactions, might be similar to Haydon's matrices in Haydon (2000) Ecology

}

inter_ocean=inter_mat*k_coast2ocean

if(intra_only){
	for(i in 2:length(sp)){
		for(j in 1:(i-1)){
			inter_mat[i,j]=0
			inter_mat[j,i]=0
			inter_ocean[i,j]=0
			inter_ocean[j,i]=0
		}
	}
}

print(paste("Start",Sys.time()))
tmpA=matrix(rnorm(length(sp)^2,0,0.1),length(sp),length(sp))
val_optim_coast=optim(tmpA,f_to_optimize_A,control=list(maxit=10000000),method="CG",B=inter_mat,N=aN) #c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")
#Nelder-Mead: about 30 seconds to convergence
#BFGS: can reach 8 minutes
#CG: between 10 and 30 seconds
#L-BFGS-B does not converge
#SANN: At least 15 minutes for only one
#Brent only for one dimension problem
conv_coast=val_optim_coast$convergence
print(paste("Intermediate",Sys.time()))
val_optim_ocean=optim(tmpA,f_to_optimize_A,control=list(maxit=10000000),method="CG",B=inter_ocean,N=aN)
conv_ocean=val_optim_ocean$convergence
print(paste("Stop",Sys.time()))
#Even after 100 iterations, we can find no case in which Nstar>0 (either ocean or coast)

A_coast=val_optim_coast$par
A_ocean=val_optim_ocean$par



#### Check if both interaction matrices can lead to a stable, positive equilibrium
#Compute mean growth rate
b_middle=optimize(f_to_optimize_B,T_min,T_max,293,A,interval=c(0,100))$minimum
r_mean=growth_rate(293,T_opt,B)

#### Check if both interaction matrices can lead to a stable, positive equilibrium
Nstar_coast=solve(A_coast)%*%(exp(r_mean)-1)
Nstar_ocean=solve(A_ocean)%*%(exp(r_mean)-1)

print(sum(Nstar_coast<0))
print(sum(Nstar_ocean<0))
print((conv_coast+conv_ocean))
print("#######################")
}
print(iter)
#if(any(inter_mat<(-1))|any(inter_ocean<(-1))){
#        print("Overcompensation, you may want to stop")
#}


list_inter=list(A_coast,A_ocean)

#Sinking rate
S=rbeta(length(sp),0.55,1.25)*30/100
names(S)=sp
#Manip so that Chaetoceros spp and THA have the highest sinking rate
tmp_S=S
or=order(S,decreasing=T)
tmp_S[c("CHD","CHS","THA")]=S[or[1:3]]
tmp_S[or[1:3]]=S[1:3]
S=tmp_S
names(S)=sp

#Gamma=germination*resuspension
#Resuspension
#The ones that sink the most are the ones that get the least suspended
or_dec=order(S,decreasing=T)
or_inc=order(S,decreasing=F)
tmp_S=S
tmp_S[or_dec]=S[or_inc]

resuspension=k_sediment2resuspension*tmp_S
Gamma=resuspension*germination
names(Gamma)=sp
