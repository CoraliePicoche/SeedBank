## 11/03/2020 CP: phenology of species

rm(list=ls())
graphics.off()
library(lubridate)
library(zoo)

start_duration_bloom=function(time_series_inter,date_bis,hydro_interp,threshold){ #This function is supposed to detect the date of the beginning and duration of bloom, defined as the population going above a certain threshold/quantile, with a gradient positive (beginning) or negative.
#	dates_bis=seq(dates[1],dates[length(dates)],1) #Regular time grid
#	time_series_interp=na.approx(time_series,x=dates,xout=dates_bis,na.rm=FALSE)
	gradient=diff(time_series_interp)
	ab_quantile=quantile(time_series_interp,threshold,na.rm=T)
	dates_above=which(time_series_interp>ab_quantile)
	beg=c()
	end_tmp=c()
	for(i in 1:length(dates_above)){
		if(dates_above[i]==1){
			beg=c(beg,dates_bis[dates_above[i]])
		}else{
			if(!is.na(gradient[dates_above[i]-1])){
				if((gradient[dates_above[i]-1]>0)&(time_series_interp[dates_above[i]-1]<ab_quantile)){
					beg=c(beg,dates_bis[dates_above[i]])
				}
				else if ((gradient[dates_above[i]-1]<0)&(time_series_interp[dates_above[i]-1]>ab_quantile)){
					if(end_tmp[length(end_tmp)]-dates_bis[dates_above[i]]==1){

					}else{
						end_tmp=c(end_tmp,dates_bis[dates_above[i]])
					}
	
				}
			}
		}
	}
	
	dd_tmp=dates_bis[!is.na(time_series_inter)]
	dd1=dd_tmp[1] #First date for which we have no NA
	dd2=dd_tmp[length(dd_tmp)] #Final date for which we have no NA

	#Now compute duration
	only_beg=beg[1]
	only_end=end_tmp[length(end_tmp)]
	temp_beg=c()
	duration=c()
	for(i in 2:(length(beg)-1)){
		if((beg[i+1]-beg[i])>2){
			if((as.Date(beg[i])!=dd1)&&(as.Date(beg[i+1])!=dd2)){
				duration=c(duration,beg[i]-only_beg[length(only_beg)])
			}
			only_beg=c(only_beg,beg[i+1])
			only_end=
			temp_beg=c(temp_beg,hydro_interp[dates_bis==beg[i]])	
		}
	}
	temp_beg=c(temp_beg,hydro_interp[dates_bis==beg[i]])	

	return(list(only_beg,duration,temp_beg))
}

tab_coast=read.table(paste("output/out_coast1.csv",sep=""),sep=";",dec=".")
name_spp=colnames(tab_coast)

abundances_tab=read.table(paste("param/","corres_hernandez_Auger.txt",sep=""),sep=";",header=T)
dates=as.Date(abundances_tab$Date)
abundances_tab=abundances_tab[year(dates)>=1996,name_spp]#Using data from 1996
dates=dates[year(dates)>=1996]

tab_mean=read.table("param/mean_value_for_reconstructed_abundances.txt",sep=";",dec=".",header=T)

hydro_tab=read.table("param/Augerhydro.txt",sep=";",header=T)
dates_hydro=as.Date(hydro_tab$Date)
temp=hydro_tab[year(dates_hydro)>=1996,"TEMP"]#Using data from 1996
dates_hydro=dates_hydro[year(dates_hydro)>=1996]

yy=unique(year(dates))

#pdf("output/time_series_persp_peryear.pdf",width=17.5,height=15)
par(mfrow=c(5,5),mar=c(1,3,1,3))
for(sp in name_spp){
	plot(0,0,t="n",xlab="",ylab="")
	legend(-0.75,0.5,c("observed","interpolated","temperature"),pch=c(16,NA,NA),col=c("black","blue","red"),lty=c(1),bty="n")
	for(i in 1:1){
		plot(0,0,t="n",xlab="",ylab="")#This is only to have one species per page as we have 21 years and 25 panels
		text(0,0,sp)
	}
	for(y in yy){
		dates_y=dates[year(dates)==y]
		abundances_tab_y=log10(abundances_tab[year(dates)==y,sp])
		dates_bis=seq(dates_y[1],dates_y[length(dates_y)],1) #Regular time grid
        	time_series_interp=na.approx(abundances_tab_y,x=dates_y,xout=dates_bis,na.rm=FALSE)

		ab_quantile=quantile(time_series_interp,0.5,na.rm=T)

		plot(dates_y,abundances_tab_y,t="o",pch=16,ylim=range(log10(abundances_tab[,sp]),na.rm=T),xlim=c(as.Date(paste("01/01/",y,sep=""),format="%d/%m/%Y"),as.Date(paste("31/12/",y,sep=""),format="%d/%m/%Y")))
		lines(dates_bis,time_series_interp,col="blue")
#		points(as.Date(paste("15",1:12,y,sep="/"),format="%d/%m/%Y"),log10(tab_mean[,sp]),pch=16,col="red")
		abline(h=ab_quantile,lty=2)
		dates_y_hydro=dates_hydro[year(dates_hydro)==y]

		temp_y=temp[year(dates_hydro)==y]
        	hydro_interp=na.approx(temp_y,x=dates_y_hydro,xout=dates_bis,na.rm=FALSE)

		x=start_duration_bloom(time_series_interp,dates_bis,hydro_interp,0.5)
		legend("topleft",c(paste("Deb",month(as.Date((x[[1]][1])))),paste("Nb:",length(x[[1]])),paste("Dur:",format(mean(x[[2]]),digits=1)),paste("Temp:",format(min(x[[3]]),digits=1),"/",format(max(x[[3]]),digits=1))),pch=NA,bty="n",col="black",inset=c(0.,0))
#		abline(v=x[[1]],col="black",lty=2)
		
		par(new = TRUE)
		plot(dates_y_hydro,temp_y,col="red",t="l",axes=F,bty="n",xlab="",ylab="")
		axis(side=4, at = pretty(range(temp,na.rm=T)))
#		abline(h=x[[4]],lty=2,col="red")

		print(c(paste("Deb",month(as.Date((x[[1]][1])))),paste("Nb:",length(x[[1]])),paste("Dur:",mean(x[[2]])),paste("Temp:",format(min(x[[3]],digits=1)),"/",format(max(x[[3]]),digits=1))))

	}
	for(i in 1:2){
		plot(0,0,t="n",xlab="",ylab="")#This is only to have one species per page as we have 21 years and 25 panels
		text(0,0,sp)
	}
}

#dev.off()
