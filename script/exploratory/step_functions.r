##########
# 17/12/2019 CP: define the functions we will use in the code
# 05/02/2020 CP: test for ssh
# 21/03/2020 CP: added the Eppley curve and cleaned useless variables that were not used anymore
##########

growth_rate_SV=function(temp,T_opt,B){ #from Scranton and Vasseur 2016 
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
		rtmp[i]=a_r*ftmp[i]*exp(E_r*(temp-t_0)/(k*temp*t_0))
	}
	print(paste("ftmp",ftmp))
	print(paste("metabolism",a_r*exp(E_r*(temp-t_0)/(k*temp*t_0))))
	print(paste("rtmp",rtmp))
	return(rtmp)
}

growth_rate_Bissinger=function(temp,irradiance){
	tmp=irradiance*0.81*exp(0.0631*temp) #This is supposed to be specific growth rate already
	return(tmp)	
}

growth_rate_Eppley=function(temp,irradiance){
	tmp=irradiance*0.59*exp(0.0633*temp) #This is supposed to be specific growth rate already
	return(tmp)	
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
        metabolism=0.81*exp(0.0631*(temp-273))*0.5
        for(i in 1:length(T_opt)){
                if(temp<=T_opt[i]){
                        ftmp[i]=exp(-(abs(temp-T_opt[i]))^3/B[i])
                }else{
                        ftmp[i]=exp(-5*(abs(temp-T_opt[i]))^3/B[i])
                }
                rtmp[i]=ftmp[i]*metabolism
        }
	print(paste("ftmp",ftmp))
	print(paste("metabolism",metabolism))
	print(paste("rtmp",rtmp))
        return(rtmp)
}


step1=function(n_t,list_inter,temp,M,mort,correct,model="BH",threshold=0.001,fixed_growth=NA,gr="Bissinger",irradiance=NA,T_opt=NA,B=NA){
	tmp=matrix(NA,dim(n_t)[1],dim(n_t)[2])
	colnames(tmp)=names(n_t)
	for(i in 1:2){
		mat_pos=mat_neg=list_inter[[i]]
		mat_pos[mat_pos<0]=0
		mat_neg[mat_neg>0]=0
		if(model=="BH"){
			growth=NA
			if(gr=="SV"){
				growth=growth_rate_SV(temp,T_opt,B)
			}else if(gr=="SV_Bissinger"){
				growth=growth_rate_noMTE_Bissinger(temp,T_opt,B)
			}else if(gr=="B"){
				growth=growth_rate_Bissinger(temp-273,irradiance)
			}else if(gr=="fixed"){
				growth=fixed_growth
			}
			tmp[i,]=exp(growth+correct)*n_t[i,]/pmax(threshold,1+list_inter[[i]]%*%n_t[i,]) - mort*n_t[i,]#We can also use the minus sign as 1-list_inter to make sure we interprete the interactions the right way.
#			tmp[i,]=exp(growth+correct)*n_t[i,]/pmax(threshold,1+mat_pos%*%n_t[i,]) - mort*n_t[i,]#We only take into account competitive interactions to see if it can compensate for the very high growth rates
			#print(paste("Growth",exp(growth+correct)*n_t[i,]))
			#print(paste("Inter",pmax(threshold,1+list_inter[[i]]%*%n_t[i,])))
			#print(paste("Loss",mort*n_t[i,]))
#			tmp[i,tmp[i,]<0]=0.0001
		}else if(model=="Martorell"){
################### This is the formula from Martorell
		for(j in 1:dim(tmp)[2]){
			prod_pos=1
			sum_neg=0
			for(k in 1:dim(tmp)[2]){
				#prod_pos=prod_pos*((1+mat_pos[j,k]*n_t[i,k])^0.1)
				if(mat_pos[j,k]>0){
				prod_pos=prod_pos*((1+n_t[i,k])^0.1)
				}
				sum_neg=sum_neg-mat_neg[j,k]*n_t[i,k]
			}
			tmp[i,j]=exp(growth_rate_SV(temp,T_opt[j],B[j]))*n_t[i,j]*prod_pos/pmax(threshold,1+sum_neg) 
		}
		}
################## End of the formula from Martorell
	}
	tmp[3,]=n_t[3,]*(1-M)
#	print(tmp[1,])
	if(sum(c(tmp)<0)>0){
	print(tmp[1,])
	print(tmp[2,])
	print(tmp[3,])
	stop()}
	return(tmp)
}

step2=function(n_t,S,Gamma,e){
	tmp=matrix(NA,dim(n_t)[1],dim(n_t)[2])
	tmp[1,]=n_t[1,]*(1-S-e)+Gamma*n_t[3,]+e*n_t[2,]
	tmp[2,]=n_t[2,]*(1-S-e)+e*n_t[1,]
	tmp[3,]=n_t[3,]*(1-Gamma)+S*n_t[1,]
	return(tmp)
}