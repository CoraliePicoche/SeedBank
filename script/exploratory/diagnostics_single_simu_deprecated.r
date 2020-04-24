############
## 22/04/20 CP: old diagnostics used to debug models without delays. Not used anymore
############

pdf(paste("timeseries_all_in_one.pdf",sep=""),width=16,height=16)
par(mfrow=c(3,1))

plot(id,transfo_N_coast[id,1],col=colo[1],t="o",pch=apch[1],ylim=range(transfo_N_coast[id,]),xaxt="n",ylab="Coast",xlab="",lty=alty[1])
for(i in 2:10){
points(id,transfo_N_coast[id,i],col=colo[i],t="o",pch=apch[i],lty=alty[i])
}

plot(id,transfo_N_ocean[id,1],col=colo[1],t="o",pch=apch[1],ylim=range(transfo_N_ocean[id,]),xaxt="n",ylab="Ocean",xlab="",lty=alty[1])
for(i in 2:10){
points(id,transfo_N_ocean[id,i],col=colo[i],pch=apch[i],t="o",lty=alty[i])
}

plot(id,transfo_N_seed[id,1],col=colo[1],t="p",pch=apch[1],ylim=range(transfo_N_seed[id,]),ylab="Seed",xlab="time",lty=alty[1])
for(i in 2:10){
points(id,transfo_N_seed[id,i],col=colo[i],pch=apch[i],t="o",lty=alty[i])
}

legend("bottomright",sp,col=colo,pch=apch,lty=alty)

dev.off()
###Growthes
#WARNING: for now, uses a file, corres_hernandez_Auger.txt, that has been removed in the new version of diagnostics it should be implemented back (as well as the interpolation of the time series, that can be copied from another file. Use seed(42), of course
if(1==0){
pdf(paste("growth_rate_daily.pdf",sep=""),width=10,height=10)
par(mfrow=c(4,3))
#id_hot=(temp>293)[id[1:(length(id)-1)]]
for(i in 1:length(sp)){
        growth_rate_coast=diff(log(tab_coast[id,i]))
        plot(log(tab_coast[id[1:(length(id)-1)],i]),growth_rate_coast,pch=16,xlab="log abundance",ylab="growth")
#        points(log(N_coast[id[1:(length(id)-1)],i][id_hot]),growth_rate_coast[id_hot],pch=16,col="red")
#        points(log(N_coast[id[1:(length(id)-1)],i][!id_hot]),growth_rate_coast[!id_hot],pch=16,col="blue")
}
dev.off()
pdf(paste("growth_rate_vs_observation_biweekly.pdf",sep=""),width=10,height=10)
par(mfrow=c(4,3))

######Take back here for growth rate every two weeks

#Simulation
id_every_2_weeks=min(id)+seq(1:26)*14
tab_coast_every_2_weeks=matrix(NA,nrow=length(id_every_2_weeks),ncol=length(name_spp))
abundances_every_2_weeks=matrix(NA,nrow=length(id_every_2_weeks),ncol=length(name_spp))
for(j in 1:(length(id_every_2_weeks)-1)){
	j=j+1
	tab_coast_every_2_weeks[j,]=as.numeric(log(tab_coast[id_every_2_weeks[j+1],])-log(tab_coast[id_every_2_weeks[j],]))
	abundances_every_2_weeks[j,]=as.numeric(log(tab_coast[id_every_2_weeks[j+1],]))
}


#Observation
obs=diff(log(tab_plankton))
abundance_obs=log(tab_plankton)
abundance_obs=abundance_obs[1:(nrow(abundance_obs)-1),]

for(i in 1:length(sp)){
	limix=range(c(abundances_every_2_weeks[,i]),c(abundance_obs[,i]),na.rm=T)
	#limiy=range(c(tab_coast_every_2_weeks[,i]),c(obs[,i]),na.rm=T)
#	limiy=range(c(tab_coast_every_2_weeks[,i]),na.rm=T)
	limiy=c(-2,2)

        plot(abundance_obs[,i],obs[,i],pch=1,col="red",xlab="log abundance",ylab="growth",ylim=limiy,xlim=limix)
        points(abundances_every_2_weeks,tab_coast_every_2_weeks,pch=16,col="blue")
}
dev.off()
} #End of 1==0 for growth rate. Might take that back later.


###SAD
#Completely arbitrary (once more): count species [0,1000],[1001,5000],[5001,10000],[10001,...]
pdf("SAD.pdf",width=10,height=15)
par(mfrow=c(4,3),mar=c(3,3,1,1))

day_in_year=0
end=0
deb=0
vec=rep(c(1,0),6)
for(m in 1:12){
#Observations
tab_categories=rep(NA,4)
tab_categories[1]=sum(tab_mean[m,]<=1000)
tab_categories[2]=sum(tab_mean[m,]>1000&tab_mean[1,]<=5000)
tab_categories[3]=sum(tab_mean[m,]>5000&tab_mean[1,]<=10000)
tab_categories[4]=sum(tab_mean[m,]>10000)
plot(1:4,tab_categories+0.1,t="p",pch=16,col="black",xaxt="n",xlab="Nb individuals",ylab="Nb species",lwd=4,main=paste("Month",m),ylim=c(0,7),cex=1.5)
axis(1,at=1:4,labels=c("[0,1000]","]1000,5000]","]5000,10000]","[10001,...]"))
legend("topright",c("Observations","Simulations"),col=c("black","blue"),pch=16,bty="n")

#Simulations

deb=end+1
if(m!=2){
end=deb+29+vec[m]
}else{
end=deb+28
}
day_in_year=end

tab_sim=apply(tab_coast[id[deb:end],],2,mean)
tab_categories[1]=sum(tab_sim<=1000)
tab_categories[2]=sum(tab_sim>1000&tab_sim<=5000)
tab_categories[3]=sum(tab_sim>5000&tab_sim<=10000)
tab_categories[4]=sum(tab_sim>10000)
points(1:4,tab_categories,pch=16,col="blue",cex=1.5)

}
dev.off()

###Output denominator
comp_coast=read.table("compet_coast.csv",sep=";",dec=".")
comp_coast_final=comp_coast[id,]
comp_ocean=read.table("compet_ocean.csv",sep=";",dec=".")
comp_ocean_final=comp_ocean[id,]

pdf("interaction_effect_coast.pdf",width=7.5,height=7.5)
par(mfrow=c(2,1))
##first 6 are centric diatoms
plot(1:366,rep(NA,366),t="n",ylim=c(0.25,1.25),xlab="Day of the year",ylab="BH denominator")
abline(h=1.0)
for(i in 1:3){
points(1:366,comp_coast_final[,i],pch=16,t="o",col=colo[i],cex=0.25)
}
legend("bottomleft",name_spp[1:3],col=colo[1:3],pch=16,lty=1,bty="n")

plot(1:366,rep(NA,366),t="n",ylim=c(0.25,1.25),xlab="Day of the year",ylab="BH denominator")
abline(h=1.0)
for(i in 4:6){
points(1:366,comp_coast_final[,i],pch=16,t="o",col=colo[i],cex=0.25)
}
legend("bottomleft",name_spp[4:6],col=colo[4:6],pch=16,lty=1,bty="n")

##Last 5 are pennate diatoms and dinoflagellates
plot(1:366,rep(NA,366),t="n",ylim=c(0.25,1.25),xlab="Day of the year",ylab="BH denominator")
abline(h=1.0)
for(i in 7:9){
points(1:366,comp_coast_final[,i],pch=16,t="o",col=colo[i],cex=0.25)
}
legend("bottomleft",name_spp[7:9],col=colo[7:9],pch=16,lty=1,bty="n")

plot(1:366,rep(NA,366),t="n",ylim=c(0.25,1.25),xlab="Day of the year",ylab="BH denominator")
abline(h=1.0)
for(i in 10:11){
points(1:366,comp_coast_final[,i],pch=16,t="o",col=colo[i],cex=0.25)
}
legend("bottomleft",name_spp[10:11],col=colo[10:11],pch=16,lty=1,bty="n")
dev.off()

pdf("interaction_effect_ocean.pdf",width=7.5,height=7.5)
par(mfrow=c(2,1))
##first 6 are centric diatoms
plot(1:366,rep(NA,366),t="n",ylim=c(0.25,1.25),xlab="Day of the year",ylab="BH denominator")
abline(h=1.0)
for(i in 1:3){
points(1:366,comp_ocean_final[,i],pch=16,t="o",col=colo[i],cex=0.25)
}
legend("bottomleft",name_spp[1:3],col=colo[1:3],pch=16,lty=1,bty="n")

plot(1:366,rep(NA,366),t="n",ylim=c(0.25,1.25),xlab="Day of the year",ylab="BH denominator")
abline(h=1.0)
for(i in 4:6){
points(1:366,comp_ocean_final[,i],pch=16,t="o",col=colo[i],cex=0.25)
}
legend("bottomleft",name_spp[4:6],col=colo[4:6],pch=16,lty=1,bty="n")

##Last 5 are pennate diatoms and dinoflagellates
plot(1:366,rep(NA,366),t="n",ylim=c(0.25,1.25),xlab="Day of the year",ylab="BH denominator")
abline(h=1.0)
for(i in 7:9){
points(1:366,comp_ocean_final[,i],pch=16,t="o",col=colo[i],cex=0.25)
}
legend("bottomleft",name_spp[7:9],col=colo[7:9],pch=16,lty=1,bty="n")

plot(1:366,rep(NA,366),t="n",ylim=c(0.25,1.25),xlab="Day of the year",ylab="BH denominator")
abline(h=1.0)
for(i in 10:11){
points(1:366,comp_ocean_final[,i],pch=16,t="o",col=colo[i],cex=0.25)
}
legend("bottomleft",name_spp[10:11],col=colo[10:11],pch=16,lty=1,bty="n")
dev.off()
