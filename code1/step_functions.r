##########
# 17/12/2019 CP: define the functions we will use in the code
##########

growth_rate=function(temp,T_opt,B){
	#Define parameters that are fixed, not phytoplankton specific
	a_r=386
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
		rtmp[i]=a_r*ftmp[i]*exp(E_r*(temp-t_0)/(k*temp*t_0))
	}
	return(rtmp)
}

step1=function(n_t,list_inter,temp,T_opt,M,B){
	tmp=matrix(NA,dim(n_t)[1],dim(n_t)[2])
	for(i in 1:2){
		tmp[i,]=exp(growth_rate(temp,T_opt,B))*n_t[i,]/(1+list_inter[[i]]%*%n_t[i,])
	}
	tmp[3,]=n_t[3,]*(1-M)
	return(tmp)
}

step2=function(n_t,S,Gamma,e){
	tmp=matrix(NA,dim(n_t)[1],dim(n_t)[2])
	tmp[1,]=n_t[1,]*(1-S-e)+Gamma*n_t[3,]+e*n_t[2,]
	print(tmp[1,])
	tmp[2,]=n_t[2,]*(1-e)+e*n_t[1,]
	print(tmp[2,])
	tmp[3,]=n_t[3,]*(1-Gamma)+S*n_t[1,]
	print(tmp[3,])
	return(tmp)
}
