##########
# 17/12/2019 CP: define the functions we will use in the code
# 05/02/2020 CP: test for ssh
# 21/03/2020 CP: added the Eppley curve and cleaned useless variables that were not used anymore
##########

growth_rate_noMTE_Bissinger=function(temp,T_opt,B,d){ #from Scranton and Vasseur 2016 
        #Define parameters that are fixed, not phytoplankton specific
        
	#Compute r(temp)
        ftmp=rep(NA,length(T_opt))
        rtmp=rep(NA,length(T_opt))
        metabolism=0.81*exp(0.0631*(temp-273))*d
        for(i in 1:length(T_opt)){
                if(temp<=T_opt[i]){
                        ftmp[i]=exp(-(abs(temp-T_opt[i]))^3/B[i])
                }else{
                        ftmp[i]=exp(-5*(abs(temp-T_opt[i]))^3/B[i])
                }
                rtmp[i]=ftmp[i]*metabolism
        }
	#print(paste("ftmp",ftmp))
	#print(paste("metabolism",metabolism))
	#print(paste("rtmp",rtmp))
        return(rtmp)
}


step1=function(n_t,list_inter,temp,M,mort,irradiance,T_opt,B,threshold){
	tmp=matrix(NA,dim(n_t)[1],dim(n_t)[2])
	effect_compet=matrix(NA,2,dim(n_t)[2])
	colnames(tmp)=names(n_t)
	for(i in 1:2){
		growth=growth_rate_noMTE_Bissinger(temp,T_opt,B,irradiance)
		effect_compet[i,]=1+list_inter[[i]]%*%n_t[i,]
		tmp[i,]=exp(growth)*n_t[i,]/pmax(threshold,1+list_inter[[i]]%*%n_t[i,]) - mort*n_t[i,]
	}
	tmp[3,]=n_t[3,]*(1-M)
#	print(tmp)
	if(sum(c(tmp)<0)>0){
	print(tmp[1,])
	print(tmp[2,])
	print(tmp[3,])
	stop()}
	return(list(tmp,effect_compet))
}

step2=function(n_t,S,Gamma,e){
	tmp=matrix(NA,dim(n_t)[1],dim(n_t)[2])
	tmp[1,]=n_t[1,]*(1-S-e)+Gamma*n_t[3,]+e*n_t[2,]
	tmp[2,]=n_t[2,]*(1-S-e)+e*n_t[1,]
	tmp[3,]=n_t[3,]*(1-Gamma)+S*n_t[1,]
	return(tmp)
}
