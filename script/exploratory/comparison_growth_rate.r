############
#06/02/2020: CP compares SV, Bissinger, and other values from the literature

rm(list=ls())
graphics.off()
source("script/step_functions.r")

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
A=10^(3.1)/365.25
T_min=min_temp+273
T_max=max_temp+273

T_opt=temp #We are in the best possible conditions

B=optimize(f_to_optimize,T_min,T_max,T_opt,A,interval=c(0,100))$minimum

gr_min_SV=growth_rate_SV(T_min,T_opt,B)
gr_max_SV=growth_rate_SV(T_max,T_opt,B)
print(paste("Scranton and Vasseur",gr_min_SV,gr_max_SV))

##Bissinger

gr_min_Bissinger=growth_rate_Bissinger(min_temp,0.5)
gr_max_Bissinger=growth_rate_Bissinger(max_temp,0.5)
print(paste("Bissinger",gr_min_Bissinger,gr_max_Bissinger))
