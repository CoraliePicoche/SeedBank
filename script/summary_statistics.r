#### CP 15/04/20 Diagnostics on different simulations / summary statistics
#### CP 16/04/20 Turns into a function
#### CP 23/04/20 From abs(error) to (error)^2 (sqrt(mean) is done in the main calibration file) + return a vector instead of a final number to keep track of the species-specific indicator
#### CP 08/05/20 Use error to classify model 
#### CP 20/05/20 Use log10(mean) for both obs and sim  + now outputs only a vector of error. The choice of absolute error or squared error is made out of this script.

summary_statistics=function(tab_mean,tab_pheno,tab_coast,nb_year){
#Simulation
name_spp=colnames(tab_coast)
n_iter=nrow(tab_coast)
#id=seq(1+floor(n_iter/365)*365-365*nb_year,floor(n_iter/365)*365)
id=(n_iter-365*nb_year):n_iter


#Average abundance
mean_sim=apply(tab_coast[id,],2,mean)
vec_diff_abundances=log10(mean_sim+10^(-5))-log10(tab_mean[,2]+10^(-5))
#diff_abundances=sum(vec_diff_abundances)

#table_summary_per_species[name_spp,"Abundances"]=vec_diff_abundances[name_spp]
#table_summary_all_models[m,"Abundances"]=diff_abundances

transfo_N_coast=log10(tab_coast[id,]+10^(-5))
#Amplitude
amplitude=apply(transfo_N_coast,2,max)-apply(transfo_N_coast,2,min)
vec_diff_amplitude=amplitude-tab_pheno[,"Mean_amplitude"]
#diff_amplitude=sum(vec_diff_amplitude)

#table_summary_per_species[name_spp,"Amplitude"]=vec_diff_abundances[name_spp]
#table_summary_all_models[m,"Amplitude"]=diff_amplitude

#Phenology
vec_diff_season=rep(NA,length(name_spp))
names(vec_diff_season)=name_spp
val_bloom=apply(transfo_N_coast,2,median)
tab_beginning_bloom=matrix(NA,length(name_spp),nb_year)
rownames(tab_beginning_bloom)=name_spp
for(s in 1:length(name_spp)){
	if(tab_coast[id[length(id)],s]!=0){
	season=NA
	for(n in 1:nb_year){
		id_y=seq((365*(n-1)+1),(365*n))
		date_bloom=which(transfo_N_coast[id_y,s]>val_bloom[s])
		tab_beginning_bloom[s,n]=date_bloom[1]
                }
	month_beginning=mean(tab_beginning_bloom[s,])%%30
	if((month_beginning>=3&month_beginning<5)|(month_beginning>=9&month_beginning<=11)){ #blooms in spring or autumn
        	season=1
      	}else if(month_beginning>=5&month_beginning<9){ #blooms in summer
        	season=2
        }else{ #blooms in winter
        season=0
	}
	if(season==tab_pheno[name_spp[s],"Season"]){
		vec_diff_season[s]=0
	}else{
		vec_diff_season[s]=1
	}
}else{
	vec_diff_season[s]=1
}
}

#table_summary_per_species[name_spp,"Phenology"]=vec_diff_season[name_spp]
#diff_season=sum(vec_diff_season)
#table_summary_all_models[m,"Phenology"]=diff_season


#Objective function
#final=diff_abundances+diff_amplitude+diff_season
#print(final)

return(list(vec_diff_abundances,vec_diff_amplitude,vec_diff_season))
#return(c(diff_abundances,diff_amplitude,diff_season))
      }

#Function which gives the rank of each model, starting with the one that minimizes the fit criteria as much as possible
classify_model=function(mat_diff_abundances,mat_diff_amplitude,mat_diff_season){
	nspp=ncol(mat_diff_abundances)
	error_abundances=translate(mat_diff_abundances)
	error_amplitude=translate(mat_diff_amplitude)
	error_pheno=translate(mat_diff_season)
	order_abundances=rank(error_abundances)
	order_amplitude=rank(error_amplitude)
	order_pheno=rank(error_pheno)
	order_total=rank(order_abundances+order_amplitude+order_pheno)

	return(order_total)
}

translate=function(mat){
	nspp=ncol(mat)
	tmp=sqrt(apply(mat^2,1,sum)/nspp)
	return(tmp)
}
