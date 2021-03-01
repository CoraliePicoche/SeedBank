###################
#19/02/2020 CP: diagnostics and plots for the model
#22/04/2020 CP: improved diagnostic figure for single species with summary statistics
#12/05/2020 CP: Replaced average observed values by examples of years for phytoplankton dynamics
#13/05/2020 CP: Improved the zoom
###################

library(lubridate)

colo=c(rep(c("red","orange","green","blue"),2),"red","orange","orchid")
apch=c(rep(16,4),rep(17,4),rep(18,2))
alty=c(rep(1,4),rep(2,4),rep(3,2))

tab=read.table(paste("../../param/simu.csv",sep=""),sep=";",dec=".",header=T)
dataset=as.character(tab[tab[,1]=="dataset",2])

#Simulation
tab_coast=read.table(paste("out_coast.csv",sep=""),sep=";",dec=".")
tab_ocean=read.table(paste("out_ocean.csv",sep=""),sep=";",dec=".")
tab_seed=read.table(paste("out_seed.csv",sep=""),sep=";",dec=".")
name_spp=colnames(tab_coast)

#Observations
tab_mean=read.table("../../param/mean_monthly_abundance.txt",sep=";",dec=".",header=T)

#####Comparing mean abundances
nb_year=2

pop_table=read.table("../../param/mean_interpolated_abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table[,1]
pop_table=pop_table[name_spp,]
x_obs=pop_table[name_spp,"Mean_abundances"]
names(x_obs)=name_spp
n_iter=nrow(tab_coast)
id=(n_iter-365*nb_year):n_iter
#Average abundance
mean_sim=apply(tab_coast[id,],2,mean)
vec_diff_abundances_minus5=log10(mean_sim+10^(-5))-log10(x_obs+10^(-5))
vec_diff_abundances_plus1=log10(mean_sim+1)-log10(x_obs+1)

vec_diff_raw_values=(mean_sim-x_obs)/x_obs*100

sum(vec_diff_raw_values^2)

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

#pdf(paste("calibration_interaction_matrix.pdf",sep=""),width=7.5,height=10)
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
before_A_no0=before_A
before_A_no0[before_A==0]=NA
after_A_no0=after_A
after_A_no0[after_A==0]=NA
axlim=c(min(c(before_A)),quantile(c(before_A_no0),na.rm=T)[4])
aylim=c(min(c(after_A)),quantile(c(after_A_no0),na.rm=T)[4])
plot(c(before_A),c(after_A),pch=16,col="black",xlab="Before calibration",ylab="After calibration",xlim=axlim,ylim=aylim,main="Zoom")
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


####Compare observed average monthly abundance to simulation
#pdf(paste("timeseries_one_by_one_average_obs.pdf",sep=""),width=10)
par(mfrow=c(1,1))

id_observed=seq(13,12*30.5,length.out=12)+min(id)

tab_pheno=read.table("../../param/generalist_specialist_spp_added_amplitude_season.csv",sep=";",dec=".",header=T,row.names=1)

for(i in 1:length(sp)){
        aylim=range(c(transfo_N_coast[id,sp[i]]),c(transfo_N_ocean[id,sp[i]]),c(log10(tab_mean[,sp[i]]+10^(-5))))
        plot(id,transfo_N_coast[id,sp[i]],main=sp[i],col="lightblue",t="o",pch=16,lty=1,xlab="time",ylab="abundance",ylim=aylim,xaxt="n")
	axis(1,at=seq(id[1],id[length(id)],by=30),labels=seq(1,366,by=30))
        lines(id,transfo_N_ocean[id,sp[i]],col="darkblue",t="o",pch=16,lty=1)
#        lines(id,transfo_N_seed[id,i],col="brown",t="o",pch=16,lty=1)
	abline(h=log10(x_obs[sp[i]]))
	#points(id_observed,log10(tab_mean[,i]+10^(-5)),pch=17,col="red",cex=2)
	points(id_observed,log10(tab_mean[,sp[i]]+10^(-5)),pch="*",col="red")
	if(i==1){
		#legend("right",c("coast","ocean","seed","mean obs"),col=c("lightblue","darkblue","brown","red"),pch=c(16,16,16,17))
		legend("right",c("coast","ocean"),col=c("lightblue","darkblue"),pch=c(16,16,16),bty="n")
	}
	mtext(paste("Mean abundance sim",format(mean(transfo_N_coast[id,sp[i]]),digits=2),"obs",format(log10(x_obs[sp[i]]),digits=2),sep=" "),cex=0.75,side=3,line=3,adj=0)
	mtext(paste("Mean amplitude sim",format(diff(range(transfo_N_coast[id,sp[i]])),digits=2),"obs",format(tab_pheno[sp[i],"Mean_amplitude"],digits=2),sep=" "),cex=0.75,side=3,line=2,adj=0)
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

#Compare real time series to simulation
abundances_raw=read.table(paste("../../param/","raw_abundances_Auger.txt",sep=""),sep=";",header=T)
dates=as.Date(abundances_raw$Date)
abundances_raw=log10(abundances_raw[year(dates)>=1996,name_spp]+10^(-5))#Using data from 1996
dates=dates[year(dates)>=1996]

#We use only a subset of possible years for the sake of clarity
years_to_show=c(2000,2005,2015)
col_yy=c("red","orange","darkred")

#pdf(paste("timeseries_one_by_one_real_obs.pdf",sep=""),width=10)
par(mfrow=c(1,1))

for(i in 1:length(sp)){
	subset_abundance=abundances_raw[,sp[i]]
        aylim=range(c(transfo_N_coast[id,sp[i]]),c(transfo_N_ocean[id,sp[i]]),subset_abundance,na.rm=T)
        plot(id,transfo_N_coast[id,sp[i]],main=sp[i],col="lightblue",t="o",pch=16,lty=1,xlab="time",ylab="abundance",ylim=aylim,xaxt="n",lwd=1.25,cex=1.25)
        axis(1,at=seq(id[1],id[length(id)],by=30),labels=seq(1,366,by=30))
        lines(id,transfo_N_ocean[id,sp[i]],col="darkblue",t="o",pch=16,lty=1)
        
	for(id_y in 1:length(years_to_show)){
		yy=years_to_show[id_y]
		subset_abundance_yy=subset_abundance[year(dates)==yy]
		dates_subset=dates[year(dates)==yy]
		dates_subset=dates_subset-as.Date(paste(yy,"-01-01",sep=""))+id[1]
		points(dates_subset,subset_abundance_yy,col=col_yy[id_y],pch=16,cex=0.75)
		lines(dates_subset,subset_abundance_yy,col=col_yy[id_y],lty=2,cex=0.75)
	}


	if(i==1){
                #legend("right",c("coast","ocean","seed","mean obs"),col=c("lightblue","darkblue","brown","red"),pch=c(16,16,16,17))
                legend("topright",c("coast","ocean"),col=c("lightblue","darkblue"),pch=16,lty=1,bty="n")
                legend("topleft",as.character(years_to_show),col=col_yy,pch=16,lty=2,,bty="n")
        }
        mtext(paste("Mean abundance sim",format(mean(transfo_N_coast[id,sp[i]]),digits=2),"obs",format(log10(x_obs[sp[i]]),digits=2),sep=" "),cex=0.75,side=3,line=3,adj=0)
        mtext(paste("Mean amplitude sim",format(diff(range(transfo_N_coast[id,sp[i]])),digits=2),"obs",format(tab_pheno[sp[i],"Mean_amplitude"],digits=2),sep=" "),cex=0.75,side=3,line=2,adj=0)
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

#We only show simulations
#Abundance classes: let's use 3 categories depending on the maximum: CHA, LEP, SKE, THP can reach nearly 5; GUI PSE PRO PRP can reach 4;  DIT, NIT PLE are the others (around 3)
high_ab=c("CHA","LEP","SKE","THP")
middle_ab=c("GUI","PSE","PRO","PRP")
low_ab=c("DIT","NIT","PLE")
#OR I can use the modules that interact with each other 
cen=c("CHA","DIT","GUI","LEP","SKE","THP")
pen=c("NIT","PLE","PSE")
din=c("PRO","PRP")

list_to_show=list(cen,pen,din)

colo=c("lightblue","blue","darkblue","orchid","lightslateblue","lightblue4")

pdf("time_series_simulation_per_cluster_with_temperature.pdf",width=7,height=15)
temp_model=as.matrix(read.table("../../param/reconstructed_temperature_Auger_with_corrected_phase_and_amplitude.txt", row.names=1,header=T,sep=";",dec="."))-273.15
name_groupe=c("Centric diatoms","Pennate diatoms","Dinoflagellates")
let=c("a","b","c")
par(mfrow=c(3,1))
for(tt in 1:length(list_to_show)){
ydelim=range(c(transfo_N_coast[id,list_to_show[[tt]]]))
	if(tt==3){
		#par(mar=c(4,4.5,3.5,1)) #Not including temperature
		par(mar=c(4,4.5,3.5,4)) #Including temperature, therefore yaxis on the right
		xl="time"
		pos_legend="topleft"
	}else{
		#par(mar=c(3,4.5,4.5,1)) #Not including temperature
		par(mar=c(3,4.5,4.5,4)) #Including temperature (and therefore, yaxis on the right)
		xl=""
		if(tt==1){
		#pos_legend=list(x=id[1]+30,y=(ydelim[2]-ydelim[1])*0.35+ydelim[1])
		pos_legend=list(x=id[1]+100,y=(ydelim[2]-ydelim[1])*0.35+ydelim[1])
		#pos_legend="bottomright"
		}else{
		pos_legend="topleft"
		}
	}
	plot(id,rep(NA,length(id)),t="n",ylim=ydelim,xlab=xl,ylab="Log10(abundance)",xaxt="n",main=name_groupe[tt],cex.lab=1.75,cex.main=1.75,cex.axis=1.75)
	text(id[1]-8,ydelim[1]+diff(ydelim)*1.07,let[tt],las=2,xpd=NA,cex=1.75,font=2)
	if(tt==3){
        axis(1,at=seq(id[1],id[length(id)],by=60),labels=seq(1,366,by=60),cex.axis=1.75)
	}else{
        axis(1,at=seq(id[1],id[length(id)],by=60),labels=rep(NA,length(seq(id[1],id[length(id)],by=60))),cex.axis=1.75)
	}
	for(s in 1:length(list_to_show[[tt]])){
        	lines(id,transfo_N_coast[id,list_to_show[[tt]][s]],col=colo[s],t="o",pch=16,lty=1)
	}
#	legend(x=id[1]+30,y=(ydelim[2]-ydelim[1])*0.35+ydelim[1],list_to_show[[tt]],col=colo[1:length(list_to_show[[tt]])],pch=16,bty="n",pt.cex=1.75,cex=1.5)
	legend(pos_legend,list_to_show[[tt]],col=colo[1:length(list_to_show[[tt]])],pch=16,bty="n",pt.cex=1.75,cex=1.5)
	if(tt==3){
	par(new = TRUE)
	plot(id, temp_model[id], type = "l", axes = FALSE, bty = "n", xlab = "",ylab="")
	axis(side=4,cex.axis=1.75)
	mtext("Temperature",side=4,line=3,cex=1.25)
	}
}
dev.off()

list_to_show=list(high_ab,middle_ab,low_ab)
#pdf("time_series_simulation_per_ab.pdf",width=7,height=15)
par(mfrow=c(3,1))
for(tt in 1:length(list_to_show)){
        if(tt==3){
                xl="time"
        }else{
                xl=""
        }
        plot(id,rep(NA,length(id)),t="n",ylim=range(c(transfo_N_coast[id,list_to_show[[tt]]])),xlab=xl,ylab="abundance",xaxt="n")
        axis(1,at=seq(id[1],id[length(id)],by=30),labels=seq(1,366,by=30))
        for(s in 1:length(list_to_show[[tt]])){
                lines(id,transfo_N_coast[id,list_to_show[[tt]][s]],col=colo[s],t="o",pch=16,lty=1)
        }
        legend("bottomright",list_to_show[[tt]],col=colo[1:length(list_to_show[[tt]])],pch=16,bty="n")
}
dev.off()
