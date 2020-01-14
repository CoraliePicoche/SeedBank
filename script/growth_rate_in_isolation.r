#####
#13/01/2020 CP: Simplified dynamics to model growth rates only and stick to values in the literature in the absence of competition from other species
#####

rm(list=ls())
graphics.off()
source("script/step_functions.r")

#To compute the b parameter in Scranton and Vasseur (2016) definition of the growth rate
f_to_optimize=function(b,T_min,T_max,T_opt,A){
        f1=integrate(Vectorize(growth_rate),lower=T_min-5,upper=T_max+5,T_opt,b)$val
        tmp=abs(f1-A)
        return(tmp)
}

temp=293
A=10^(3.1)/365.25
T_min=288
T_max=298

#Let's start with Chaetoceros species which are supposed to grow well
T_opt=temp #We are in the best possible conditions

B=optimize(f_to_optimize,T_min,T_max,T_opt,A,interval=c(0,100))$minimum


gr=growth_rate(temp,T_opt,B)

print("N1/N0")
print(exp(gr)) #This is the value we have in papers. Should be around 1 to 2 

##### Now, we add the density-dependance,focusing on one species
fixed_growth=seq(0.1,20,length.out=1000)
#intra_sp=-1*seq(10^(-8),10^(-5),length.out=1000)
intra_sp=-0.35

mean_val=rep(NA,length(intra_sp))

n_iter=10000

theta=1.3
mean_temp=273+20
sd_temp=2.5
temp=mean_temp+theta*sd_temp*sin(2*pi*1:n_iter/365.25)+rnorm(n_iter,0,sd_temp*sqrt(1-theta^2/2))

#for(i in 1:length(intra_sp)){
for(i in 1:length(fixed_growth)){
N=rep(NA,n_iter)
N[1]=10^6

for(t in 1:(n_iter-1)){
	N[t+1]=exp(fixed_growth[i])*N[t]/(1-intra_sp*N[t]) #The minus sign is there so that negative interactions do reduce growth rates
}


mean_val[i]=mean(log10(N[(n_iter-1000):n_iter]))

}

plot(fixed_growth,mean_val)
abline(h=6)
