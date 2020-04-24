###################
#19/02/2020 CP: diagnostics and plots for the model
#22/04/2020 CP: improved diagnostic figure for single species with summary statistics
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
tab_mean=read.table("../../param/mean_monthly_abundance.txt",sep=";",dec=".",header=T)

#Variation due to quadratic programming
before=as.matrix(read.table(paste("interaction_matrix_before_calibration.csv",sep=""),sep=";",dec=".",header=T))
before_A=before[,1:(ncol(before)-1)]
before_A_nodiag=before_A
diag(before_A_nodiag)=NA
after=as.matrix(read.table(paste("interaction_matrix_after_calibration.csv",sep=""),sep=";",dec=".",header=T))
after_A=after[,1:(ncol(after)-1)]
after_A_nodiag=after_A
diag(after_A_nodiag)=NA

ratio_after=mean(abs(diag(after_A)))/mean(abs(after_A_nodiag))
ratio_before=mean(abs(diag(before_A)))/mean(abs(before_A_nodiag))

pdf(paste("calibration_interaction_matrix.pdf",sep=""),width=7.5,height=10)
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


####Two options are possible for mean values
#If we use raw values of raw_abundances, we avoid the artefacts created by the interpolation and the random value used when gaps are over 2 points in the time series, but we increase the mean value artificially as cells are counted only when they are numerous. The inverse is true when using interpolated data. This is only a matter of choice.

#abundances_tab=read.table(paste("param/","raw_abundances_Auger.txt",sep=""),sep=";",header=T)
#dates=as.Date(abundances_tab$Date)
#abundances_tab=abundances_tab[year(dates)>=1996,name_spp]#Using data from 1996
#dates=dates[year(dates)>=1996]
#x_obs=apply(abundances_tab,mean,na.rm=T)

pop_table=read.table("../../param/mean_interpolated_abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
x_obs=pop_table[name_spp,"Mean_abundances"]
names(x_obs)=name_spp

pdf(paste("timeseries_one_by_one.pdf",sep=""),width=10)
par(mfrow=c(1,1))

id_observed=seq(13,12*30.5,length.out=12)+min(id)

tab_pheno=read.table("../../param/generalist_specialist_spp_added_amplitude_season.csv",sep=";",dec=".",header=T,row.names=1)

for(i in 1:length(sp)){
        aylim=range(c(transfo_N_coast[id,i]),c(transfo_N_ocean[id,i]),c(log10(tab_mean[,i]+10^(-5))))
        plot(id,transfo_N_coast[id,i],main=sp[i],col="lightblue",t="o",pch=16,lty=1,xlab="time",ylab="abundance",ylim=aylim,xaxt="n")
	axis(1,at=seq(id[1],id[length(id)],by=30),labels=seq(1,366,by=30))
        lines(id,transfo_N_ocean[id,i],col="darkblue",t="o",pch=16,lty=1)
#        lines(id,transfo_N_seed[id,i],col="brown",t="o",pch=16,lty=1)
	abline(h=log10(x_obs[sp[i]]))
	#points(id_observed,log10(tab_mean[,i]+10^(-5)),pch=17,col="red",cex=2)
	points(id_observed,log10(tab_mean[,i]+10^(-5)),pch="*",col="red")
	if(i==1){
		#legend("right",c("coast","ocean","seed","mean obs"),col=c("lightblue","darkblue","brown","red"),pch=c(16,16,16,17))
		legend("right",c("coast","ocean"),col=c("lightblue","darkblue"),pch=c(16,16,16),bty="n")
	}
	mtext(paste("Mean abundance sim",format(mean(transfo_N_coast[id,i]),digits=2),"obs",format(log10(x_obs[sp[i]]),digits=2),sep=" "),cex=0.75,side=3,line=3,adj=0)
	mtext(paste("Mean amplitude sim",format(diff(range(transfo_N_coast[id,i])),digits=2),"obs",format(tab_pheno[sp[i],"Mean_amplitude"],digits=2),sep=" "),cex=0.75,side=3,line=2,adj=0)
	tt=tab_pheno[sp[i],"Season"]
	if(tt==0){
		text_to_write="Bloom winter"
	}else if(tt==1){
		text_to_write="Early bloom"
	}else{
		text_to_write="Late bloom"
	}
	mtext(text_to_write,cex=0.75,side=3,line=1,adj=0)
}
dev.off()

