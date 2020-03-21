## 11/03/2020 CP: phenology of species
#21/03/2020 Corrected a bug on the beginning of the time series and added histogram per species, and overall


rm(list=ls())
graphics.off()
library(lubridate)
library(zoo)

threshold=0.5

start_duration_bloom_new=function(time_series_interp,date_bis,hydro_interp,threshold){ #This function takes interpolated time series (abundance and temperature here) with corresponding dates to detect blooming beginning and end (and thus duration) and temperature for bloom. Beginning of the bloom is defined as the date when abundances goes above a certain threshold
        gradient=diff(time_series_interp)
        ab_quantile=quantile(time_series_interp,threshold,na.rm=T)
        id_dates_above=which(time_series_interp>ab_quantile)
	ll=length(time_series_interp)
	if(ll!=length(dates_bis)){
		stop("Pb in time series length")
	}else if(max(id_dates_above)>ll){
		stop("Pb in identifying bloomin periods")
	}
	#Define the beginning and ends of the bloom
        beg_tmp=c()
        end_tmp=c()
	#We check the change in abundance from i-1 to i+1, so we need to check that both indices exist
	#If the abundance is already above the threshold when the time series begins, ie if we can't have i-1
	id_notna=which(!is.na(time_series_interp))
       	if(id_dates_above[1]==id_notna[1]){
                        beg_tmp=c(beg_tmp,dates_bis[id_dates_above[1]])
			deb=2
	}else{
		deb=1
	}
	for(i in deb:(length(id_dates_above))){
		if(id_dates_above[i]!=ll){ #If i+1 does not exist
	        	if((!(is.na(gradient[id_dates_above[i]-1])))&(!(is.na(gradient[id_dates_above[i]+1])))){
				#If we go from beyond the threshold to above the threshold, it's a beginning
        	        	if((gradient[id_dates_above[i]-1]>0)&(time_series_interp[id_dates_above[i]-1]<=ab_quantile)&((time_series_interp[id_dates_above[i]+1]>=ab_quantile))){
                	                        beg_tmp=c(beg_tmp,dates_bis[id_dates_above[i]])
				#If we go from above the threshold to beyond the threshold, it's an end
                        	}else if ((gradient[id_dates_above[i]-1]<0)&(time_series_interp[id_dates_above[i]-1]>=ab_quantile)&(time_series_interp[id_dates_above[i]+1]<=ab_quantile)){
                                        end_tmp=c(end_tmp,dates_bis[id_dates_above[i]])
                        	}
			}
		}
	}
	#Check the duration of each bloom	
	dur=c()
	#Sometimes, there were too much NA to have any date of end or beginning. In this case, there is no duration possible
	if((length(end_tmp)>0)&(length(beg_tmp)>0)){
	#If each beginning has an end
	if(length(beg_tmp)==length(end_tmp)){
		dur=end_tmp-beg_tmp
	}else{
		#If there is more beginning than ends: there was a bloom at the end of the year that did not decrease enough to cross the threshold
		if(length(beg_tmp)>length(end_tmp)){
			already_used_beg=c()
			for(i in 1:length(end_tmp)){
					ok=F
					g=length(beg_tmp)
				while(!ok&g>0){
					if((beg_tmp[g]<end_tmp[i])&!(beg_tmp[g]%in%already_used_beg)){
						already_used_beg=c(already_used_beg,beg_tmp[g])
						dur=c(dur,end_tmp[i]-beg_tmp[g])
						ok=T
					}
					g=g-1
				}
			}
			
		}else{
		#If there is more ends than beginnings: the time series begin with the end of a bloom
			already_used_end=c()
                        for(i in 1:length(beg_tmp)){
				ok=F
				g=1
                                while(!ok&g<=length(end_tmp)){
                                        if((beg_tmp[i]<end_tmp[g])&!(end_tmp[g]%in%already_used_end)){
                                                already_used_end=c(already_used_end,end_tmp[g])
                                                dur=c(dur,end_tmp[g]-beg_tmp[i])
						ok=T
                                        }
					g=g+1
                                }
                        }
		}
	}
	}
	#Find the temperature at the beginning of the bloom
	temp_beg=c()
	for(b in beg_tmp){
		temp_beg=c(temp_beg,hydro_interp[dates_bis==b])	
	}
	
	return(list(beg_tmp,dur,end_tmp,temp_beg))
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

#name_spp=c("CHA")
pdf("output/time_series_persp_peryear.pdf",width=17.5,height=15)
par(mfrow=c(5,5),mar=c(3,3,3,3))
#name_spp=c("LEP")
#yy=2001
dur_final=c()
total_bloom_final=c()
temp_min_final=c()
for(sp in name_spp){
	plot(0,0,t="n",xlab="",ylab="")
	legend(-0.75,0.5,c("observed","interpolated","temperature"),pch=c(16,NA,NA),col=c("black","blue","red"),lty=c(1),bty="n")
	text(0,0,sp)
	total_bloom=c()
	duration=c()
	temp_min=c()
	for(y in yy){
		print(y)
		dates_y=dates[year(dates)==y]
		abundances_tab_y=log10(abundances_tab[year(dates)==y,sp])
		dates_bis=seq(dates_y[1],dates_y[length(dates_y)],1) #Regular time grid
        	time_series_interp=na.approx(abundances_tab_y,x=dates_y,xout=dates_bis,na.rm=FALSE)

		ab_quantile=quantile(time_series_interp,threshold,na.rm=T)

		plot(dates_y,abundances_tab_y,t="o",pch=16,ylim=range(log10(abundances_tab[,sp]),na.rm=T),xlim=c(as.Date(paste("01/01/",y,sep=""),format="%d/%m/%Y"),as.Date(paste("31/12/",y,sep=""),format="%d/%m/%Y")),main=y)
		lines(dates_bis,time_series_interp,col="blue")
#		points(as.Date(paste("15",1:12,y,sep="/"),format="%d/%m/%Y"),log10(tab_mean[,sp]),pch=16,col="red")
		abline(h=ab_quantile,lty=2)
		dates_y_hydro=dates_hydro[year(dates_hydro)==y]

		temp_y=temp[year(dates_hydro)==y]
        	hydro_interp=na.approx(temp_y,x=dates_y_hydro,xout=dates_bis,na.rm=FALSE)

		x=start_duration_bloom_new(time_series_interp,dates_bis,hydro_interp,threshold)

		total_bloom=c(total_bloom, length(x[[1]]))
		duration=c(duration,x[[2]])
		temp_min=c(temp_min,min(x[[4]]))

		if(is.null(x[[1]])){
			x[[1]]=NA
		}
		if(is.null(x[[2]])){
			x[[2]]=NA
		}
		if(length(x[[3]])==0){
			x[[3]]=NA
		}
		if(is.null(x[[4]])){
			x[[4]]=NA
		}
		legend("topleft",c(paste("Deb",month(as.Date((x[[1]][1])),label=T)),paste("Nb:",length(x[[1]])),paste("Dur:",min(x[[2]]),max(x[[2]])),paste("Temp:",format(min(x[[4]]),digits=1),"/",format(max(x[[4]]),digits=1))),pch=NA,bty="n",col="black",inset=c(0.,0))
		abline(v=x[[1]],col="black",lty=2)
		abline(v=x[[3]],col="black",lty=3)
		
		par(new = TRUE)
		plot(dates_y_hydro,temp_y,col="red",t="l",axes=F,bty="n",xlab="",ylab="")
		axis(side=4, at = pretty(range(temp,na.rm=T)))
#		abline(h=x[[4]],lty=2,col="red")

		print(c(paste("Deb",month(as.Date((x[[1]][1])),label=T)),paste("Nb:",length(x[[1]])),paste("Dur:",min(x[[2]]),max(x[[2]])),paste("Temp:",format(min(x[[4]],digits=1)),"/",format(max(x[[4]]),digits=1))))

	}
	hist(total_bloom)
	hist(duration)
	hist(temp_min)
	
	total_bloom_final=c(total_bloom_final,total_bloom)
	dur_final=c(dur_final,duration)
	temp_min_final=c(temp_min_final,temp_min)
}

dev.off()

pdf("output/summary_all_indices.pdf",width=12)
par(mfrow=c(1,3))
	hist(total_bloom_final)
	abline(v=median(total_bloom_final),col="red")
	hist(dur_final)
	abline(v=median(dur_final),col="red")
	hist(temp_min_final)
	abline(v=median(temp_min_final),col="red")
dev.off()
