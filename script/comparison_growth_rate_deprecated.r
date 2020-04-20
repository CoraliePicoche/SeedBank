############
#06/02/2020: CP compares SV, Bissinger, and other values from the literature
#21/03/2020: Removed the half-correction of Bissinger (there is no reason, at least in their paper, to consider the daylength effect on the growth rate) and corrected the legend to make it more clear. Also added the Eppley curve in the figure.

rm(list=ls())
graphics.off()
source("script/step_functions.r")

correct=0.25

tab_temp=read.table("param/Augerhydro.txt",sep=";",header=T)
temp=tab_temp[,"TEMP"]
min_temp=min(temp,na.rm=T)
max_temp=max(temp,na.rm=T)

#To compute the b parameter in Scranton and Vasseur (2016) definition of the growth rate
f_to_optimize=function(b,T_min,T_max,T_opt,A){
        f1=integrate(Vectorize(growth_rate_SV),lower=T_min,upper=T_max,T_opt,b)$val
        tmp=abs(f1-A)
        return(tmp)
}

temp=273+15
A_init=10^(3.1)/365.25
A_correct=15
T_min=min_temp+273
T_max=max_temp+273

T_min=0+273
T_max=35+273

T_opt=temp #We are in the best possible conditions

B_init=optimize(f_to_optimize,T_min,T_max,T_opt,A_init,interval=c(0,100000000))$minimum
print(B_init)

B_correct=optimize(f_to_optimize,T_min,T_max,T_opt,A_correct,interval=c(0,100000000))$minimum


T_opt_low=min_temp+273
T_opt_high=max_temp+273
B_low=optimize(f_to_optimize,T_min,T_max,T_opt_low,A_correct,interval=c(0,100000000))$minimum
B_high=optimize(f_to_optimize,T_min,T_max,T_opt_high,A_correct,interval=c(0,100000000))$minimum


#gr_min_SV=growth_rate_SV(T_min,T_opt,B)
#gr_max_SV=growth_rate_SV(T_max,T_opt,B)
#print(paste("Scranton and Vasseur",gr_min_SV,gr_max_SV))

seq_temp=seq(T_min,T_max,length.out=1000)
val_growth_SV_init=c()
val_growth_SV_corrected=c()
val_growth_SV_high=c()
val_growth_SV_low=c()
val_growth_Biss=c()
val_growth_Eppley=c()

for(t in seq_temp){
	val_growth_SV_init=c(val_growth_SV_init,growth_rate_SV(t,T_opt,B_init))
	val_growth_SV_corrected=c(val_growth_SV_corrected,growth_rate_SV(t,T_opt,B_correct)+correct)
	val_growth_SV_low=c(val_growth_SV_low,growth_rate_SV(t,T_opt_low,B_low))
	val_growth_SV_high=c(val_growth_SV_high,growth_rate_SV(t,T_opt_high,B_high))
	val_growth_Biss=c(val_growth_Biss,growth_rate_Bissinger(t-273,1))
	val_growth_Eppley=c(val_growth_Eppley,growth_rate_Eppley(t-273,1))
}

##Bissinger
gr_min_Bissinger=growth_rate_Bissinger(min_temp,1)
gr_max_Bissinger=growth_rate_Bissinger(max_temp,1)
print(paste("Bissinger",gr_min_Bissinger,gr_max_Bissinger))

#Val in the literature
gr_min_lit=0.3
gr_max_lit=1.78

pdf("output/compare_growth_rate_rawBissinger.pdf")
plot(seq_temp-273,val_growth_SV_init,ylim=c(0,3),t="l",lwd=3,xlab="Temperature",ylab="Growth rate (d-1)",xlim=c(0,30),col="grey")
lines(seq_temp-273,val_growth_SV_corrected,lwd=3,col="black")
lines(seq_temp-273,val_growth_SV_low+correct,ylim=c(0,2),lwd=3,lty=3,col="blue")
lines(seq_temp-273,val_growth_SV_high+correct,ylim=c(0,2),lwd=3,lty=3,col="red")
#abline(h=gr_min_Bissinger,lty=2)
#abline(h=gr_max_Bissinger,lty=2)
lines(seq_temp-273,val_growth_Biss,lty=2)
lines(seq_temp-273,val_growth_Eppley,lty=2,col="orchid")
abline(h=gr_min_lit,lty=1)
abline(h=gr_max_lit,lty=1)
legend("topleft",c("SV","SV corrected, A=15","low T_opt SV corrected","high T_opt SV corrected","Bissinger","Eppley","Literature"),lty=c(1,1,3,3,2,2,1),lwd=c(3,3,3,3,2,1,1,1),bty="n",col=c("grey","black","blue","red","black","orchid","black"))
dev.off()
