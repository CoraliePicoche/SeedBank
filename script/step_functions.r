##########
# 17/12/2019 CP: define the functions we will use in the code
##########

growth_rate=function(temp,T_opt,B){
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
	return(rtmp)
}

step1=function(n_t,list_inter,temp,T_opt,M,B,model="BH",threshold=0.001,fixed_growth=NA){
	tmp=matrix(NA,dim(n_t)[1],dim(n_t)[2])
	colnames(tmp)=names(n_t)
	for(i in 1:2){
		mat_pos=mat_neg=list_inter[[i]]
		mat_pos[mat_pos<0]=0
		mat_neg[mat_neg>0]=0
		if(model=="BH"){
			tmp[i,]=exp(growth_rate(temp,T_opt,B))*n_t[i,]/pmax(threshold,1-list_inter[[i]]%*%n_t[i,]) #The minus sign is there so that negative interactions do reduce growth rates
		}else if(model=="fixed"){
			tmp[i,]=fixed_growth*n_t[i,]/pmax(threshold,1-list_inter[[i]]%*%n_t[i,])
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
			tmp[i,j]=exp(growth_rate(temp,T_opt[j],B[j]))*n_t[i,j]*prod_pos/pmax(threshold,1+sum_neg)
		}
		}
################## End of the formula from Martorell
	}
	tmp[3,]=n_t[3,]*(1-M)
	if(sum(c(tmp)<0)>0){
	print(tmp[1,])
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
