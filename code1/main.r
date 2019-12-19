###############
#17/12/2019 CP Main for the model
###############

graphics.off()
rm(list=ls())

source("code1/param_definition.r")

n_iter=100
N=array(NA,dim=c(n_iter,3,length(sp)))
N[1,,]=rep(10,length(sp)*3)



temp=rnorm(n_iter,273+20,5)

for(t in 1:(n_iter-1)){
	Ntmp=step1(N[t,,],list_inter,temp[t],T_opt,M,B)
### For now, we have negative abundances from the beginning!
### I think the problem is due to the interaction matrix. I need to check how it is taken into account in a classical BH model, and with positive interactions. Maybe check that it can be stable??
	N[t+1,,]=step2(Ntmp,S,Gamma,e)
}
