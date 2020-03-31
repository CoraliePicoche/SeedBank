###################
#19/02/2020 CP: diagnostics and plots for the model
###################

colo=c(rep(c("red","orange","green","blue"),2),"red","orange","orchid")
apch=c(rep(16,4),rep(17,4),rep(18,2))
alty=c(rep(1,4),rep(2,4),rep(3,2))

tab=read.table(paste("simu.csv",sep=""),sep=";",dec=".",header=T)
dataset=as.character(tab[tab[,1]=="dataset",2])

#Simulation
tab_coast=read.table(paste("out_coast.csv",sep=""),sep=";",dec=".")
tab_ocean=read.table(paste("out_ocean.csv",sep=""),sep=";",dec=".")
tab_seed=read.table(paste("out_seed.csv",sep=""),sep=";",dec=".")
name_spp=colnames(tab_coast)

#Observations
tab_mean=read.table("../../param/mean_value_for_reconstructed_abundances.txt",sep=";",dec=".",header=T)

#Variation due to quadratic programming
before=as.matrix(read.table(paste("matrix_A_before_quad.csv",sep=""),sep=";",dec=".",header=T))
before_A=before[,1:(ncol(before)-1)]
before_A_nodiag=before_A
diag(before_A_nodiag)=NA
after=as.matrix(read.table(paste("matrix_A_after_quad.csv",sep=""),sep=";",dec=".",header=T))
after_A=after[,1:(ncol(after)-1)]
after_A_nodiag=after_A
diag(after_A_nodiag)=NA

ratio_after=mean(abs(diag(after_A)))/mean(abs(after_A_nodiag))
ratio_before=mean(abs(diag(before_A)))/mean(abs(before_A_nodiag))

pdf(paste("calibration_A.pdf",sep=""),width=7.5,height=10)
par(mfrow=c(3,1),mar=c(4,4,1,1))
plot(c(before_A),c(after_A),pch=16,col="black",xlab="",ylab="After calibration")
abline(a=0,b=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
points(diag(before_A),diag(after_A),pch=16,col="red")
legend("topleft",c("Inter","Intra"),col=c("black","red"),pch=16)
text(1,0,paste("Before ",ratio_before,"\nAfter ",ratio_after,sep=""))

diff=(after_A/before_A)
id_diff=which(diff>10)
if(length(id_diff)>0){
id=0
val=c()
for(i in 1:ncol(after_A)){
	for(j in 1:ncol(after_A)){
		id=id+1
		tmp=which(id_diff==id)
		if(length(tmp)==1){
			val=c(val,paste("(",colnames(after_A)[i],",",colnames(after_A)[j],")",sep=""))
		}
	}
}
}

#plot(c(before_A),c(after_A),pch=16,col="black",xlab="Before calibration",ylab="After calibration",xlim=c(min(c(before_A)),10^-4),ylim=c(min(c(after_A)),10^-4),main="Zoom")
plot(c(before_A),c(after_A),pch=16,col="black",xlab="Before calibration",ylab="After calibration",xlim=c(-5*10^-5,5*10^-5),ylim=c(-7.5*10^(-5),7.5*10^-5),main="Zoom")
abline(a=0,b=1)
abline(h=0,lty=2)
abline(v=0,lty=2)
points(diag(before_A),diag(after_A),pch=16,col="red")

points(before_A[id_diff],after_A[id_diff],col="orange",pch="*",cex=2)

convert_c_to_tab=c()
for (i in 1:ncol(before_A)){
	deb=as.character(i)
	for(j in 1:ncol(before_A)){
		end=as.character(j)
		convert_c_to_tab=c(convert_c_to_tab,paste("(",deb,",",end,")",sep=""))
	}
}



plot(1:length(before_A),after_A/before_A,pch=16,col="black",xlab="",ylab="Ratio after/before",xaxt="n")
print(id_diff)
if(length(id_diff)>0){
text(id_diff,diff[id_diff],val,pos=1)
}
seq_axis=seq(1,length(before_A),by=4)
abline(h=10,lty=3)
axis(1,at=seq_axis,labels=convert_c_to_tab[seq_axis])

dev.off()

#Time series for the last year
sp=colnames(tab_coast)
n_iter=nrow(tab_coast)

transfo_N_coast=log10(tab_coast+10^(-5))
transfo_N_ocean=log10(tab_ocean+10^(-5))
transfo_N_seed=log10(tab_seed+10^(-5))

id=(n_iter-365):n_iter

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


pop_table=read.table("../../param/abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
x_obs=pop_table[name_spp,"Mean_abundances"]
names(x_obs)=name_spp

pdf(paste("timeseries_one_by_one.pdf",sep=""),width=10)
par(mfrow=c(1,1))

id_observed=seq(13,12*30.5,length.out=12)+min(id)


for(i in 1:length(sp)){
        aylim=range(c(transfo_N_coast[id,i]),c(transfo_N_ocean[id,i]),c(transfo_N_seed[id,i]),c(log10(tab_mean[,i]+10^(-5))))
        plot(id,transfo_N_coast[id,i],main=sp[i],col="lightblue",t="o",pch=16,lty=1,xlab="time",ylab="abundance",ylim=aylim,xaxt="n")
	axis(1,at=seq(id[1],id[length(id)],by=30),labels=seq(1,366,by=30))
        lines(id,transfo_N_ocean[id,i],col="darkblue",t="o",pch=16,lty=1)
        lines(id,transfo_N_seed[id,i],col="brown",t="o",pch=16,lty=1)
	abline(h=log10(x_obs[sp[i]]))
	points(id_observed,log10(tab_mean[,i]+10^(-5)),pch=17,col="red",cex=2)
	if(i==1){
		#legend("right",c("coast","ocean","seed","mean obs"),col=c("lightblue","darkblue","brown","red"),pch=c(16,16,16,17))
		legend("right",c("coast","ocean","seed"),col=c("lightblue","darkblue","brown"),pch=c(16,16,16))
	}
}
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
pdf("SAD.pdf",width=7.5,height=5)
par(mfrow=c(1,2))
#Observations
#In Winter
tab_categories=rep(NA,4)
tab_categories[1]=sum(tab_mean[1,]<=1000)
tab_categories[2]=sum(tab_mean[1,]>1000&tab_mean[1,]<=5000)
tab_categories[3]=sum(tab_mean[1,]>5000&tab_mean[1,]<=10000)
tab_categories[4]=sum(tab_mean[1,]>10000)
plot(1:4,tab_categories+0.1,t="p",pch=16,col="black",xaxt="n",xlab="Nb individuals",ylab="Nb species",lwd=4,main="January",ylim=c(0,7),cex=1.5)
axis(1,at=1:4,labels=c("[0,1000]","]1000,5000]","]5000,10000]","[10001,...]"))
legend("topright",c("Observations","Simulations"),col=c("black","blue"),pch=16,bty="n")

#Simulations
tab_sim=apply(tab_coast[id[1:31],],2,mean)
tab_categories[1]=sum(tab_sim<=1000)
tab_categories[2]=sum(tab_sim>1000&tab_sim<=5000)
tab_categories[3]=sum(tab_sim>5000&tab_sim<=10000)
tab_categories[4]=sum(tab_sim>10000)
points(1:4,tab_categories,pch=16,col="blue",cex=1.5)


#In Summer
tab_categories=rep(NA,4)
tab_categories[1]=sum(tab_mean[7,]<=1000)
tab_categories[2]=sum(tab_mean[7,]>1000&tab_mean[1,]<=5000)
tab_categories[3]=sum(tab_mean[7,]>5000&tab_mean[1,]<=10000)
tab_categories[4]=sum(tab_mean[7,]>10000)
plot(1:4,tab_categories+0.1,t="p",pch=16,col="black",xaxt="n",xlab="Nb individuals",ylab="Nb species",lwd=4,main="July",ylim=c(0,7),cex=1.5)
axis(1,at=1:4,labels=c("[0,1000]","]1000,5000]","]5000,10000]","[10001,...]"))

tab_sim=apply(tab_coast[id[180:210],],2,mean)
tab_categories[1]=sum(tab_sim<=1000)
tab_categories[2]=sum(tab_sim>1000&tab_sim<=5000)
tab_categories[3]=sum(tab_sim>5000&tab_sim<=10000)
tab_categories[4]=sum(tab_sim>10000)
points(1:4,tab_categories,pch=16,col="blue",cex=1.5)
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

