####### 21/03/2020: CP checking the exact effect of width and height of the niche when varying A only for Scranton & Vasseur (2016)
#For the same temperature, different A


rm(list=ls())
graphics.off()

source("script/step_functions.r")
source("script/infer_interaction_matrix_growth_rate.r")
library(RColorBrewer)


growth_rate_debug=function(temp,T_opt,B){ #from Scranton and Vasseur 2016 
        #Define parameters that are fixed, not phytoplankton specific
        a_r=386/365.25
        E_r=0.467
        k=8.6173324*10^(-5)
        t_0=293

        #Compute r(temp)
        ftmp=rep(NA,length(T_opt))
        rtmp=rep(NA,length(T_opt))
        for(i in 1:length(T_opt)){
                if(temp<=T_opt[i]){
                        ftmp[i]=exp(-(abs(temp-T_opt[i]))^3/B[i])
                }else{
                        ftmp[i]=exp(-5*(abs(temp-T_opt[i]))^3/B[i])
                }
		pp=a_r*exp(E_r*(temp-t_0)/(k*temp*t_0))
                rtmp[i]=a_r*ftmp[i]*exp(E_r*(temp-t_0)/(k*temp*t_0))
        }
        return(list(ftmp,rtmp,pp))
}

growth_rate_noMTE_Bissinger=function(temp,T_opt,B){ #from Scranton and Vasseur 2016 
        #Define parameters that are fixed, not phytoplankton specific
        a_r=386/365.25
        E_r=0.467
        k=8.6173324*10^(-5)
        t_0=293

        #Compute r(temp)
        ftmp=rep(NA,length(T_opt))
        rtmp=rep(NA,length(T_opt))
	metabolism=0.81*exp(0.0631*(temp-273))
        for(i in 1:length(T_opt)){
                if(temp<=T_opt[i]){
                        ftmp[i]=exp(-(abs(temp-T_opt[i]))^3/B[i])
                }else{
                        ftmp[i]=exp(-5*(abs(temp-T_opt[i]))^3/B[i])
                }
                rtmp[i]=ftmp[i]*metabolism
        }
        return(list(ftmp,rtmp))
}


nspp=10
T_opt=15+273
T_min=0+273
T_max=30+273
seq_T=seq(T_min,T_max,length.out=100)

#A=seq(5,15,length.out=nspp)
#B=rep(NA,nspp)

B=seq(200,20000,length.out=nspp)
acol=rainbow(nspp,start=0.7,end=0.025)


pdf("output/further_notes_on_SV.pdf",width=12.5)
par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(0,0,xlim=range(seq_T),ylim=c(0,1.25),ylab="growth rate",xlab="temperature",t="n")
abline(v=T_opt)
for(a in 1:nspp){
        rtmp_temp=rep(NA,length(seq_T))
        for(t in 1:length(seq_T)){
                plou=growth_rate_debug(seq_T[t],T_opt,B[a])
                rtmp_temp[t]=plou[[2]]
                
        }
        lines(seq_T,rtmp_temp,col=acol[a],lty=1)
}

legend("topleft",legend=B,col=acol,lty=1,bty="n",lwd=2)

plot(0,0,xlim=range(seq_T),ylim=c(0,1.25),ylab="",xlab="temperature",t="n")
abline(v=T_opt)
for(a in c(1,nspp)){
	ftmp_temp=rep(NA,length(seq_T))
	rtmp_temp=rep(NA,length(seq_T))
	temp_var=rep(NA,length(seq_T))
	for(t in 1:length(seq_T)){
		plou=growth_rate_debug(seq_T[t],T_opt,B[a])
		ftmp_temp[t]=plou[[1]]
		rtmp_temp[t]=plou[[2]]
		temp_var[t]=plou[[3]]
		
	}
	lines(seq_T,ftmp_temp,col=acol[a],lty=2)
	lines(seq_T,rtmp_temp,col=acol[a],lty=1)
	lines(seq_T,temp_var,col=acol[a],lty=4)
}


legend("topleft",legend=c("Growth rate","Niche part","Metabolism part"),col="black",lty=c(1,2,4),bty="n",lwd=2)
dev.off()

gr_min_lit=0.3
gr_max_lit=1.78

pdf("output/Eppley_Bissinger.pdf")
par(mfrow=c(1,1))
plot(0,0,xlim=range(seq_T),ylim=c(0,3),ylab="growth rate",xlab="temperature",t="n")
abline(v=T_opt)
for(a in 1:nspp){
        rtmp_temp=rep(NA,length(seq_T))
        rtmp_temp_noMTE=rep(NA,length(seq_T))
        rtmp_temp_other=rep(NA,length(seq_T))
        for(t in 1:length(seq_T)){
                plou=growth_rate_noMTE_Bissinger(seq_T[t],T_opt,B[a])
                rtmp_temp_noMTE[t]=plou[[2]]
                plou=growth_rate_debug(seq_T[t],T_opt,B[a])
                rtmp_temp[t]=plou[[2]]
        }
        lines(seq_T,rtmp_temp,col=acol[a],lty=2)
        lines(seq_T,rtmp_temp_noMTE,col=acol[a],lty=1)
	if(a==nspp){
        rtmp_temp_other=rep(NA,length(seq_T))
        	for(t in 1:length(seq_T)){
                plou=growth_rate_noMTE_Bissinger(seq_T[t],T_opt+10,B[a])
                rtmp_temp_other[t]=plou[[2]]
		}
	}
}

#abline(h=gr_min_lit,lty=1)
#abline(h=gr_max_lit,lty=1)

legend("topleft",legend=c("Bissinger-noMTE","Original"),col="black",lty=c(1,2),bty="n",lwd=2)
dev.off()
