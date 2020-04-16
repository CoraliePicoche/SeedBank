#### CP 15/04/20 Diagnostics on different simulations / summary statistics
#### CP 16/04/20 Turns into a function

summary_statistics=function(tab_mean,tab_pheno,tab_coast,nb_year){
#Simulation
name_spp=colnames(tab_coast)
n_iter=nrow(tab_coast)
id=(n_iter-365*nb_year):n_iter

transfo_N_coast=log10(tab_coast[id,]+10^(-5))

#Average abundance
mean_sim=apply(transfo_N_coast,2,mean)
vec_diff_abundances=abs(mean_sim-log10(tab_mean[,2]+10^(-5)))
diff_abundances=sum(vec_diff_abundances)

#table_summary_per_species[name_spp,"Abundances"]=vec_diff_abundances[name_spp]
#table_summary_all_models[m,"Abundances"]=diff_abundances

#Amplitude
amplitude=apply(transfo_N_coast,2,max)-apply(transfo_N_coast,2,min)
vec_diff_amplitude=abs(amplitude-tab_pheno[,"Mean_amplitude"])
diff_amplitude=sum(vec_diff_amplitude)

#table_summary_per_species[name_spp,"Amplitude"]=vec_diff_abundances[name_spp]
#table_summary_all_models[m,"Amplitude"]=diff_amplitude

#Phenology
vec_diff_season=rep(NA,length(name_spp))
names(vec_diff_season)=name_spp
val_bloom=apply(transfo_N_coast,2,median)
tab_beginning_bloom=matrix(NA,length(name_spp),nb_year)
rownames(tab_beginning_bloom)=name_spp
for(s in 1:length(name_spp)){
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
}
#table_summary_per_species[name_spp,"Phenology"]=vec_diff_season[name_spp]
diff_season=sum(vec_diff_season)
#table_summary_all_models[m,"Phenology"]=diff_season


#Objective function
final=diff_abundances+diff_amplitude+diff_season
print(final)

#return(c(vec_diff_abundances,vec_diff_amplitude,vec_diff_season))
return(c(diff_abundances,diff_amplitude,diff_season))
}
