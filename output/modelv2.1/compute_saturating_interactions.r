#### CP 08/04/2020 Compute saturating interactions 
#### CP 29/04/2020 Corrected a bug in computing the maximum competitive coefficients: before, took the sum of all competition effects on all species. Should actually be the mean of the maximum effect on any species
#### CP 25/05/2020 Removed the path to step functions

library(lubridate)
source("../../script/infer_interaction_matrix_growth_rate.r")
source("../../script/matrix_MAR_clean.r")


MAR2saturation=function(B,N_mean,N_max,O_y,ratio_pos){
	#Step 1 translate MAR to classical BH
	A=MAR2BH(B,N_mean)
	name_spp=colnames(A)
	#Step 2: compute maximum competition strength and maximum facilitation strength
	#Competition
#	a_C=sum(abs(A)%*%N_max)
	a_C=mean(abs(A)%*%N_max)

	#Facilitation
	#Compute max growth rate
	a_F=(1/ratio_pos)*(1/O_y-1-(1-ratio_pos)*a_C)

	#Step 3: compute corresponding half-saturation coefficients
	#+We will need a matrix of these values for the step_function
	half_saturation=matrix(0,nspp,nspp)
	max_interactions=matrix(0,nspp,nspp)
	colnames(half_saturation)=name_spp
	rownames(half_saturation)=name_spp
	colnames(max_interactions)=name_spp
	rownames(max_interactions)=name_spp
	for(s1 in 1:nspp){
		for(s2 in 1:nspp){
			if(A[s1,s2]!=0){
				max_interactions[s1,s2]=a_C*(A[s1,s2]>0)+a_F*(A[s1,s2]<0)
				half_saturation[s1,s2]=max_interactions[s1,s2]/A[s1,s2]
			}
		}
	}
	return(list(max_interactions,half_saturation))
}
