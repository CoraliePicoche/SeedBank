#### CP 08/04/2020 Compute saturating interactions 

library(lubridate)
source("../../script/infer_interaction_matrix_growth_rate.r")
source("step_functions.r")
source("../../script/matrix_MAR_clean.r")


load(paste("../../param/Auger_pencen_null_regular_common_MO.RData",sep=""))
name_spp=colnames(cis$call$model$B)
B_matrix=clean_matrix(cis,signif=F)
rownames(B_matrix)=colnames(B_matrix)=name_spp
nspp=length(name_spp)

abundances_tab=read.table(paste("../../param/corres_hernandez_Auger.txt",sep=""),sep=";",header=T)
#abundances_tab=read.table(paste("param/corres_hernandez_Auger.txt",sep=""),sep=";",header=T)
dates=as.Date(abundances_tab$Date)
abundances_tab=abundances_tab[year(dates)>=1996,name_spp]#Using data from 1996
N_max=apply(abundances_tab,2,max,na.rm=T)
N_mean=apply(abundances_tab,2,mean,na.rm=T)

max_g=growth_rate_noMTE_Bissinger(273+25,273+25,1500,0.5) #25 degree is the maximum temperature we observe in Auger; 1500 is the maximum value we have for the generalist. We place ourselves in the best case scenario
O_y=1

MAR2saturation=function(B,N_mean,N_max,max_g,O_y)
	#Step 1 translate MAR to classical BH
	A=MAR2BH(B,N_mean)

	#Step 2: compute maximum competition strength and maximum facilitation strength
	#Competition
	a_C=sum(abs(A)%*%N_max)

	#Facilitation
	#Compute max growth rate
	a_F=1/0.7*(O_y*exp(max_g)-1-0.3*a_C)

	#Step 3: compute corresponding half-saturation coefficients
	#+We will need a matrix of these values for the step_function
	half_saturation=matrix(NA,nspp,nspp)
	max_interactions=matrix(NA,nspp,nspp)
	for(s1 in 1:nspp){
		for(s2 in 1:nspp){
			max_interactions[s1,s2]=a_C*(A[s1,s2]>0)+a_F*(A[s1,s2]<0)
			half_saturation[s1,s2]=max_interactions[s1,s2]/A[s1,s2]
		}
	}
	return(max_interactions,half_saturation)

print(MAR2saturation(B_matrix,N_mean,N_max,max_g,O_y))
