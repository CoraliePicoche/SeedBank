###################
#19/02/2020 CP: diagnostics and plots for the model
###################

rm(list=ls())
graphics.off()

library('lubridate')
library('zoo')

set.seed(42)
n_simulation=2

colo=c(rep(c("red","orange","green","blue"),2),"red","orange")
apch=c(rep(16,4),rep(17,4),rep(18,2))
alty=c(rep(1,4),rep(2,4),rep(3,2))

tab=read.table(paste("param/simu",n_simulation,".csv",sep=""),sep=";",dec=".",header=T)
dataset=as.character(tab[tab[,1]=="dataset",2])

#Simulation
tab_coast=read.table(paste("output/out_coast",n_simulation,".csv",sep=""),sep=";",dec=".")
tab_ocean=read.table(paste("output/out_ocean",n_simulation,".csv",sep=""),sep=";",dec=".")
tab_seed=read.table(paste("output/out_seed",n_simulation,".csv",sep=""),sep=";",dec=".")
name_spp=colnames(tab_coast)

#Observations
timestep=14
consecutif=2
abundances_tab=read.table(paste("param/","corres_hernandez_",dataset,".txt",sep=""),sep=";",header=T)
dates=as.Date(abundances_tab$Date)
abundances_tab=abundances_tab[year(dates)>=1996,name_spp]#Using data from 1996
dates=dates[year(dates)>=1996]
dates_bis=seq(dates[1],dates[length(dates)],timestep) #Regular time grid
tab_plankton=na.approx(abundances_tab,maxgap=consecutif,x=dates,xout=dates_bis,na.rm=FALSE) #Interpolation over regular time grid
#Replace missing values
for (s in name_spp){
	tab_plankton[is.na(tab_plankton[,s]),s]=runif(sum(is.na(tab_plankton[,s])),0,min(tab_plankton[,s],na.rm=TRUE))
	}
tab_mean=matrix(NA,ncol=length(name_spp),nrow=12,dimnames=list(c("01","02","03","04","05","06","07","08","09","10","11","12"),name_spp))
for(i in 1:nrow(tab_mean)){
	tab_mean[i,]=apply(tab_plankton[month(dates)==as.numeric(rownames(tab_mean)[i]),],2,mean)
}
rownames(tab_mean)=1:12
write.table(tab_mean,"param/mean_value_for_reconstructed_abundances.txt",sep=";",dec=".",row.names=T)

#Variation due to quadratic programming
before=as.matrix(read.table(paste("output/matrix_A_before_quad_",n_simulation,".csv",sep=""),sep=";",dec=".",header=T))
before_A=before[,1:(ncol(before)-1)]
before_A_nodiag=before_A
diag(before_A_nodiag)=NA
after=as.matrix(read.table(paste("output/matrix_A_after_quad_",n_simulation,".csv",sep=""),sep=";",dec=".",header=T))
after_A=after[,1:(ncol(after)-1)]
after_A_nodiag=after_A
diag(after_A_nodiag)=NA

ratio_after=mean(abs(diag(after_A)))/mean(abs(after_A_nodiag))
ratio_before=mean(abs(diag(before_A)))/mean(abs(before_A_nodiag))

pdf(paste("output/calibration_A_",n_simulation,".pdf",sep=""),width=7.5,height=10)
par(mfrow=c(3,1),mar=c(4,4,1,1))
plot(c(before_A),c(after_A),pch=16,col="black",xlab="",ylab="After calibration")
abline(a=0,b=1)
points(diag(before_A),diag(after_A),pch=16,col="red")
legend("topleft",c("Inter","Intra"),col=c("black","red"),pch=16)
text(1,0,paste("Before ",ratio_before,"\nAfter ",ratio_after,sep=""))

plot(c(before_A),c(after_A),pch=16,col="black",xlab="Before calibration",ylab="After calibration",xlim=c(min(c(before_A)),10^-4),ylim=c(min(c(after_A)),10^-4),main="Zoom")
abline(a=0,b=1)
points(diag(before_A),diag(after_A),pch=16,col="red")

convert_c_to_tab=c()
for (i in 1:ncol(before_A)){
	deb=as.character(i)
	for(j in 1:ncol(before_A)){
		end=as.character(j)
		convert_c_to_tab=c(convert_c_to_tab,paste("(",deb,",",end,")",sep=""))
	}
}

diff=(after_A/before_A)
id_diff=which(diff>10)

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

plot(1:length(before_A),after_A/before_A,pch=16,col="black",xlab="",ylab="Ratio after/before",xaxt="n")
print(id_diff)
#####text(id_diff,diff[id_diff],val,pos=1)
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

pdf(paste("output/timeseries_all_in_one_",n_simulation,".pdf",sep=""),width=16,height=16)
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


pop_table=read.table("param/phenology_Auger.csv",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
x_obs=pop_table[name_spp,"Mean_abundances"]
names(x_obs)=name_spp

pdf(paste("output/timeseries_one_by_one",n_simulation,".pdf",sep=""),width=10)
par(mfrow=c(1,1))

id_observed=seq(13,12*30.5,length.out=12)+min(id)


for(i in 1:length(sp)){
        aylim=range(c(transfo_N_coast[id,i]),c(transfo_N_ocean[id,i]),c(transfo_N_seed[id,i]),c(log10(tab_mean[,i]+10^(-5))))
        plot(id,transfo_N_coast[id,i],main=sp[i],col="lightblue",t="o",pch=16,lty=1,xlab="time",ylab="abundance",ylim=aylim,xaxt="n")
	axis(1,at=seq(id[1],id[length(id)],by=30),labels=seq(1,366,by=30))
        lines(id,transfo_N_ocean[id,i],col="darkblue",t="o",pch=16,lty=1)
        lines(id,transfo_N_seed[id,i],col="brown",t="o",pch=16,lty=1)
	abline(h=log10(x_obs[sp[i]]))
	#points(id_observed,log10(tab_mean[,i]+10^(-5)),pch=17,col="red",cex=2)
	if(i==1){
		#legend("right",c("coast","ocean","seed","mean obs"),col=c("lightblue","darkblue","brown","red"),pch=c(16,16,16,17))
		legend("right",c("coast","ocean","seed"),col=c("lightblue","darkblue","brown"),pch=c(16,16,16))
	}
}
dev.off()


#Comparison growth rates
#Observation
gr_observed=apply(log(tab_plankton),2,diff)

pdf(paste("output/growth_rate_",n_simulation,"_daily.pdf",sep=""),width=10,height=10)
par(mfrow=c(4,3))
#id_hot=(temp>293)[id[1:(length(id)-1)]]
for(i in 1:length(sp)){
        growth_rate_coast=diff(log(tab_coast[id,i]))
        plot(log(tab_coast[id[1:(length(id)-1)],i]),growth_rate_coast,pch=16,xlab="log abundance",ylab="growth")
#        points(log(N_coast[id[1:(length(id)-1)],i][id_hot]),growth_rate_coast[id_hot],pch=16,col="red")
#        points(log(N_coast[id[1:(length(id)-1)],i][!id_hot]),growth_rate_coast[!id_hot],pch=16,col="blue")
}
dev.off()
pdf(paste("output/growth_rate_",n_simulation,"vs_observation_biweekly.pdf",sep=""),width=10,height=10)
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

