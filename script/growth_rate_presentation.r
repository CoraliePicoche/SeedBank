######## CP 29/04/2020 Explore different growth rate function for SI (based on two previous scripts: comparison_growth_rate_deprecated.r and variation_in_area.r).

rm(list=ls())
graphics.off()

growth_rate_SV=function(temp,T_opt,B){ #from Scranton and Vasseur 2016 
        #Define parameters that are fixed, not phytoplankton specific
        a_r=386/365.25
        E_r=0.467
        k=8.6173324*10^(-5)
        t_0=293

        #Compute r(temp)
        ftmp=rep(NA,length(T_opt))
        rtmp=rep(NA,length(T_opt))
        pp=a_r*exp(E_r*(temp-t_0)/(k*temp*t_0))
        for(i in 1:length(T_opt)){
                if(temp<=T_opt[i]){
                        ftmp[i]=exp(-(abs(temp-T_opt[i]))^3/B[i])
                }else{
                        ftmp[i]=exp(-5*(abs(temp-T_opt[i]))^3/B[i])
                }
                rtmp[i]=pp*ftmp[i]
        }
        return(list(ftmp,rtmp,pp))
}

growth_rate_Bissinger=function(temp,irradiance){
        tmp=irradiance*0.81*exp(0.0631*temp) 
        return(tmp)
}

growth_rate_Eppley=function(temp,irradiance){
        tmp=irradiance*0.59*exp(0.0633*temp)
        return(tmp)
}

growth_rate_Bissinger_complete=function(temp,T_opt,B,irradiance){ #from Scranton and Vasseur 2016 
        #Compute r(temp)
        ftmp=rep(NA,length(T_opt))
        rtmp=rep(NA,length(T_opt))
        pp=growth_rate_Bissinger(temp-273.15,irradiance)
        for(i in 1:length(T_opt)){
                if(temp<=T_opt[i]){
                        ftmp[i]=exp(-(abs(temp-T_opt[i]))^3/B[i])
                }else{
                        ftmp[i]=exp(-5*(abs(temp-T_opt[i]))^3/B[i])
                }
                rtmp[i]=pp*ftmp[i]
        }
        return(list(ftmp,rtmp,pp))
}

B=1300
T_opt=273.15+15
T_min=0+273.15
T_max=30+273.15
seq_temp=seq(T_min,T_max,length.out=1000)

val_growth_SV_init=c()
val_growth_SV_meta=c()
val_growth_SV_niche=c()
val_growth_Biss=c()
val_growth_Eppley=c()

for(t in seq_temp){
	plou=growth_rate_SV(t,T_opt,B)
        val_growth_SV_init=c(val_growth_SV_init,plou[[2]])
        val_growth_SV_meta=c(val_growth_SV_meta,plou[[3]])
        val_growth_SV_niche=c(val_growth_SV_niche,plou[[1]])
        val_growth_Biss=c(val_growth_Biss,growth_rate_Bissinger(t-273.15,0.5))
        val_growth_Eppley=c(val_growth_Eppley,growth_rate_Eppley(t-273.15,0.5))
}

#Val in the literature
gr_min_lit=0.3
gr_max_lit=1.78

pdf("output/compare_growth_rate_rawBissinger.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(4,4.5,1,1))
plot(0,0,xlim=range(seq_temp)-273.15,ylim=c(0,1.25),ylab=expression(paste("Growth rate (day"^"-1",")",sep="")),xlab="Temperature",t="n")
lines(seq_temp-273.15,val_growth_SV_init,col="black",lty=1,lwd=1.5)
lines(seq_temp-273.15,val_growth_SV_meta,col="blue",lty=1)
lines(seq_temp-273.15,val_growth_SV_niche,col="aquamarine3",lty=1,lwd=1.1)


legend("topleft",legend=c("Growth rate","Niche part","Metabolism part"),col=c("black","aquamarine3","blue"),lty=c(1),bty="n",lwd=1)


id=seq(1,1000,length.out=20)
plot(seq_temp[id]-273.15,val_growth_SV_meta[id],ylim=c(0.25,2.0),t="p",xlab="Temperature",ylab="",xlim=c(0,30),col="black",pch="+")

lines(seq_temp-273.15,val_growth_Biss,lty=1,col="orchid")
lines(seq_temp-273.15,val_growth_Eppley,lty=1,col="steelblue1")
abline(h=gr_min_lit,lty=2)
abline(h=gr_max_lit,lty=2)
legend(0,1.75,c("Eppley+MTE","Eppley","Bissinger","Literature"),bty="n",col=c("black","steelblue1","orchid","black"),lty=c(NA,1,1,2),pch=c("+",NA,NA,NA))
dev.off()


nspp=10
B=c(seq(10,50,length.out=5),seq(180,1500,length.out=5))
acol=rainbow(nspp,start=0.7,end=0.025)

pdf("output/further_notes_on_SV.pdf",width=10,height=5)
par(mfrow=c(1,2),mar=c(4,4.5,1,1))
plot(0,0,xlim=range(seq_temp)-273.15,ylim=c(0,1.25),ylab=expression(paste("Growth rate (day"^"-1",")",sep="")),xlab="Temperature",t="n")
abline(v=T_opt-273)
for(a in 1:nspp){
        rtmp_temp=rep(NA,length(seq_temp))
        for(t in 1:length(seq_temp)){
                plou=growth_rate_Bissinger_complete(seq_temp[t],T_opt,B[a],0.5)
                rtmp_temp[t]=plou[[2]]

        }
        lines(seq_temp-273.15,rtmp_temp,col=acol[a],lty=1)
}

legend("topleft",legend=B,col=acol,lty=1,bty="n",lwd=2)

id=seq(1,1000,length.out=75)
plot(0,0,xlim=range(seq_temp)-273.15,ylim=c(0,1.25),ylab="",xlab="Temperature",t="n")
abline(v=T_opt-273.15)
for(a in c(nspp,1)){
        ftmp_temp=rep(NA,length(seq_temp))
        rtmp_temp=rep(NA,length(seq_temp))
        temp_var=rep(NA,length(seq_temp))
        for(t in 1:length(seq_temp)){
                plou=growth_rate_Bissinger_complete(seq_temp[t],T_opt,B[a],0.5)
                ftmp_temp[t]=plou[[1]]
                rtmp_temp[t]=plou[[2]]
                temp_var[t]=plou[[3]]

        }
        points(seq_temp[id]-273.15,ftmp_temp[id],col=acol[a],pch=16)
        lines(seq_temp-273.15,rtmp_temp,col=acol[a],lty=1)
        lines(seq_temp-273.15,temp_var,col=acol[a],lty=4)
}


legend("topleft",legend=c("Growth rate","Niche part","Metabolism part"),col="black",lty=c(1,NA,4),bty="n",lwd=2,pch=c(NA,16,NA))
dev.off()


