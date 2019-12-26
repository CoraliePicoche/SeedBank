###############
#17/12/2019 CP Main for the model
###############

graphics.off()
rm(list=ls())

source("code1/param_definition.r")

n_iter=10000
N=array(NA,dim=c(n_iter,3,length(sp)))
N[1,,]=rep(10^6,length(sp)*3)

theta=1.3
mean_temp=273+20
sd_temp=2.5
temp=mean_temp+theta*sd_temp*sin(2*pi*1:n_iter/365.25)+rnorm(n_iter,0,sd_temp*sqrt(1-theta^2/2))

for(t in 1:(n_iter-1)){
	Ntmp=step1(N[t,,],list_inter,temp[t],T_opt,M,B)
	N[t+1,,]=step2(Ntmp,S,Gamma,e)
}

colo=rainbow(10)

par(mfrow=c(1,3))
plot(1:n_iter,log10(N[,1,1]),col=colo[1],t="p",pch=16,xlim=c(n_iter-1000,n_iter),ylim=range(log10(10^(-5)+N[,1,])))
for(i in 2:10){
points(1:n_iter,log10(N[,1,i]),col=colo[i],t="p",pch=16)
}

plot(1:n_iter,log10(N[,2,1]),col=colo[1],t="p",pch=16,xlim=c(n_iter,n_iter-1000),ylim=range(log10(10^(-5)+N[,2,])))
for(i in 2:10){
points(1:n_iter,log10(N[,2,i]),col=colo[i],pch=16)
}

plot(1:n_iter,log10(N[,3,1]),col=colo[1],t="p",pch=16,xlim=c(n_iter,n_iter-1000),ylim=range(log10(10^(-5)+N[,3,])))
for(i in 2:10){
points(1:n_iter,log10(N[,3,i]),col=colo[i],pch=16)
}


