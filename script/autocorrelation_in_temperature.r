#### CP 25/03/2020 : Looking for autocorrelation in temperature without seasonality

rm(list=ls())
graphics.off()
library(lubridate)

evt_tab=read.table(paste("param/Augerhydro.txt",sep=""),sep=";",header=T)
evt_tab$Date=as.Date(evt_tab$Date)
evt_tab=evt_tab[year(evt_tab$Date)>=1996,]

dates=evt_tab$Date
temp=evt_tab[,"TEMP"]
#We need to build a simulation for temperatures as we don't have enough data in the real dataset
min_temp=min(temp,na.rm=T)
max_temp=max(temp,na.rm=T)
sd_temp=sd(temp,na.rm=T)
mean_temp=mean(temp,na.rm=T)

#Interpolation
#We will do it year after year because it seems that there is too much NA from one year to another, which leads to high autocorrelation
yy=unique(year(dates))

y=yy[1]
dates_bis=seq(dates[year(dates)==y][1],dates[year(dates)==y][length(dates[year(dates)==y])],1) #Regular time grid
temp_bis=approx(temp,x=dates,xout=dates_bis)$y
t=1:length(dates_bis)
for(y in yy[2:length(yy)]){
last_date=dates_bis[length(dates_bis)]
first_date=dates[year(dates)==y][1]
dates_bis=c(dates_bis,seq(dates[year(dates)==y][1],dates[year(dates)==y][length(dates[year(dates)==y])],1)) #Regular time grid
id=which(year(dates)==y)
idi2=which(year(dates_bis)==y)
temp_bis=c(temp_bis,approx(temp[id],x=dates[id],xout=dates_bis[idi2])$y)
t=c(t,t[length(t)]+as.numeric(first_date-last_date)+seq(0,as.numeric(dates_bis[length(dates_bis)]-first_date),1))
}
#t=1:length(temp_bis)
omega=2*pi/(365.25)
sini=sin(omega*t)
cosi=cos(omega*t)
#gap_temp=!is.na(temp_bis)
season=lm(temp_bis~sini+cosi)

residuals=lm(temp_bis~season$fitted.values)$residuals

par(mfrow=c(1,2))
#id1=which(year(dates_bis)>1999&year(dates_bis)<2003)
id1=which(year(dates_bis)==2000)
plot(dates_bis[id1],season$coefficients[1]+season$coefficients[2]*sini[id1]+season$coefficients[3]*cosi[id1],t="l")
#id2=which(year(evt_tab$Date)>1999&year(evt_tab$Date)<2003)
id2=which(year(evt_tab$Date)==2000)
points(evt_tab$Date[id2],temp[id2],col="blue",pch=16)


plot(dates_bis[id2],residuals[id2])
