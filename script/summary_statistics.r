#### CP 15/04/20 DIagnostics on different simulations / summary statistics

rm(list=ls())
graphics.off()

nb_year=2
#Simulation
tab_coast=read.table(paste("output/modelv2.1/out_coast.csv",sep=""),sep=";",dec=".")
name_spp=colnames(tab_coast)
n_iter=nrow(tab_coast)
id=(n_iter-365*nb_year):n_iter

transfo_N_coast=log10(tab_coast[id,]+10^(-5))



#Observations
tab_mean=read.table("param/abundances_Auger.txt",sep=",",dec=".",header=T)
tab_pheno=read.table("param/generalist_specialist_spp_added_amplitude_season.csv",sep=";",dec=".",header=T,row.names=1)

#Average abundance
mean_sim=apply(transfo_N_coast,2,mean)
vec_diff_abundances=abs(mean_sim-log10(tab_mean[,2]+10^(-5)))
diff_abundances=sum(vec_diff_abundances)

#Amplitude
amplitude=apply(transfo_N_coast,2,max)-apply(transfo_N_coast,2,min)
vec_diff_amplitude=abs(amplitude-tab_pheno[,"Mean_amplitude"])
diff_amplitude=sum(vec_diff_amplitude)

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
diff_season=sum(vec_diff_season)


#Objective function
final=diff_abundances+diff_amplitude+diff_season
print(final)