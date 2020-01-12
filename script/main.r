###############
#17/12/2019 CP Main for the model
###############

graphics.off()
rm(list=ls())

source("script/param_definition.r")

n_iter=10000
N=array(NA,dim=c(n_iter,3,length(sp)),dimnames=list(NULL,c("coast","ocean","seed"),sp))
N[1,,]=rep(10^6,length(sp)*3)

theta=1.3
mean_temp=273+20
sd_temp=2.5
temp=mean_temp+theta*sd_temp*sin(2*pi*1:n_iter/365.25)+rnorm(n_iter,0,sd_temp*sqrt(1-theta^2/2))

for(t in 1:(n_iter-1)){
	Ntmp=step1(N[t,,],list_inter,temp[t],T_opt,M,B)
	N[t+1,,]=step2(Ntmp,S,Gamma,e)
}

write.table(N[,1,],"output/out_coast.csv",sep=";",dec=".")
write.table(N[,2,],"output/out_ocean.csv",sep=";",dec=".")
write.table(N[,3,],"output/out_seed.csv",sep=";",dec=".")

source("script/diagnostics.r")
