###############
#17/12/2019 CP Main for the model
###############

graphics.off()
rm(list=ls())

source("script/param_definition.r")
source("script/quadprog_for_inter_growth.r")

n_iter=10000
N=array(NA,dim=c(n_iter,3,length(sp)),dimnames=list(NULL,c("coast","ocean","seed"),sp))
N[1,,]=rep(10^3,length(sp)*3)

###False temperature
theta=1.3
#theta=0 #No seasonality
mean_temp=273+20
sd_temp=2.5
temp_model=mean_temp+theta*sd_temp*sin(2*pi*1:n_iter/365.25)+rnorm(n_iter,0,sd_temp*sqrt(1-theta^2/2))

###Reconstructed temperature
tab_temp=read.table("param/Augerhydro.txt",sep=";",header=T)
temp=tab_temp[,"TEMP"]
min_temp=min(temp,na.rm=T)
max_temp=max(temp,na.rm=T)
sd_temp=sd(temp,na.rm=T)
mean_temp=mean(temp,na.rm=T)

temp_model=273+mean_temp+theta*sd_temp*sin(2*pi*1:n_iter/365.25)+rnorm(n_iter,0,sd_temp*sqrt(1-theta^2/2))


mean_irr=12/24
irr=mean_irr+sin(2*pi*1:n_iter/365.25)*4/24

for(t in 1:(n_iter-1)){
#	Ntmp=step1(N[t,,],list_inter,temp[t],T_opt,M,B,model="fixed",fixed_growth=exp(0.52))
	Ntmp=step1(N[t,,],list_inter,temp_model[t],T_opt,M,B,model="BH",gr="Bissinger",threshold=0.001,irradiance=irr[t])
#	Ntmp=step1(N[t,,],list_inter,temp[t],T_opt,M,B,model="Martorell")
	N[t+1,,]=step2(Ntmp,S,Gamma,e)
}

write.table(N[,1,],"output/out_coast.csv",sep=";",dec=".")
write.table(N[,2,],"output/out_ocean.csv",sep=";",dec=".")
write.table(N[,3,],"output/out_seed.csv",sep=";",dec=".")

source("script/diagnostics.r")
