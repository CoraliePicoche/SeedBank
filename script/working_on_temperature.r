############
#CP 01/06/2020 Checking if temperature is ok
###########

rm(list=ls())
graphics.off()

set.seed(42)

library(lubridate)

#Data to use (Auger)
evt_tab=read.table(paste("param/Augerhydro.txt",sep=""),sep=";",header=T)
evt_tab$Date=as.Date(evt_tab$Date)
temp=evt_tab$TEMP[year(evt_tab$Date)>=1996]
dates=evt_tab$Date[year(evt_tab$Date)>=1996]
min_temp=min(temp,na.rm=T)
max_temp=max(temp,na.rm=T)
sd_temp=sd(temp,na.rm=T)
mean_temp=mean(temp,na.rm=T)

#We need to build a simulation for temperatures as we don't have enough data in the real dataset
n_iter=10000
theta=1.412
#temp_model=273.15+mean_temp+theta*sd_temp*sin(2*pi*1:n_iter/365+pi)+rnorm(n_iter,0,sd_temp*sqrt(1-theta^2/2)) ###WHAT IT SHOULD BE

pdf("output/comparison_temperature.pdf",width=16,height=16)
par(mfrow=c(4,4))
for(y in unique(year(dates))){
	temp_tmp=temp[year(dates)==y]
	dates_tmp=dates[year(dates)==y]
	temp_model=mean_temp+theta*sd_temp*sin(2*pi*(1:365)/365+pi+pi/3)+rnorm(365,0,sd_temp*sqrt(1-theta^2/2))
	date_in_day=dates_tmp-as.Date(paste(y,"01-01",sep="-"))
	plot(date_in_day,temp_tmp,pch=16,col="red",t="p",cex=1.5)
	lines(temp_model)
	print(mean(abs(diff(temp_model))))
}
dev.off()


#Checking that the thermal amplitude is not completely off
print("##############################")
daily_temp=read.table("../../Plankton/data/raw/Meteo/meteo_1987-2012.csv",sep=";",header=T,dec=",")
print(mean(abs(diff(daily_temp$T.eau.eyrac)),na.rm=T))


