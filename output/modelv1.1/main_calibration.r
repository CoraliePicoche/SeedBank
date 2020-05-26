#################
###CP 18/02/2020: Clean version of the main function
###CP 16/04/2020: Version for main calibration
###CP 22/04/2020: Corrected lots of bugs
###CP 23/04/2020: Use RMSE instead of abs(error) when doing summary statistics, use a vector instead of a final value to keep track of the species-specific values
###CP 12/05/2020: Implemented the rank-based choice of model instead of using direct values of summary statistics
###CP 25/05/2020: Corrected the use of summary statistics (squaring the error is not applied in this function anymore), added the original simulation in the list of simulation to take into account when choosing the right community matrix, corrected paths to scripts and added persistence in the ocean
#################

rm(list=ls())
graphics.off()
set.seed(40) #Warning: results and even feasibility can be highly seed-sensitive

source("../../script/matrix_MAR_clean.r")
source("../../script/common_scripts/step_functions.r")
source("../../script/infer_interaction_matrix_growth_rate.r")
source("../../script/summary_statistics.r")

####Calibration
value=c(0.5,0.9,1.1,2) #Parameter space
nb_simu=1000
nb_year=2
nb_persistence=6*30 #Persistence is defined as the fact that the species has disappeared for the last 6 months (taking into account extreme case in which species could get back from the deads from the seed bank in summer, when germination begins)

#Fixed parameters
tab=read.table("../../param/simu.csv",sep=";",dec=".",header=T)
cyst_mortality=as.numeric(as.character(tab[tab[,1]=="cyst_mortality",2]))
cyst_burial=as.numeric(as.character(tab[tab[,1]=="cyst_burial",2]))
M=cyst_mortality+cyst_burial
k_coast2ocean=as.numeric(as.character(tab[tab[,1]=="k_coast2ocean",2]))
morta=as.numeric(as.character(tab[tab[,1]=="loss_rate",2]))
e=as.numeric(as.character(tab[tab[,1]=="exchange",2]))
germination=as.numeric(as.character(tab[tab[,1]=="germination",2]))
resuspension=as.numeric(as.character(tab[tab[,1]=="resuspension",2]))
Gamma=resuspension*germination
S_max=as.numeric(as.character(tab[tab[,1]=="max_sinking",2]))
quad_prog=tab[tab[,1]=="quad_prog",2]
n_iter=as.numeric(as.character(tab[tab[,1]=="n_iter",2]))
temp_germin=as.numeric(as.character(tab[tab[,1]=="germin_threshold",2]))
a_d=as.numeric(as.character(tab[tab[,1]=="daylength",2]))
threshold=as.numeric(as.character(tab[tab[,1]=="threshold",2]))

#Data to use (Auger)
a=as.matrix(read.table("../../param/reconstructed_temperature_Auger.txt", row.names=1,header=T,sep=";",dec="."))
temp_model=a[1:n_iter]


B_matrix=read.table("../../param/community_matrix_B_Auger.csv",sep=";",dec=".",row.names=1,header=T)
name_spp=colnames(B_matrix)
nspp=length(name_spp)
#load(paste("../../param/Auger_pencen_null_regular_common_MO.RData",sep=""))
#name_spp=colnames(cis$call$model$B)
#B_matrix=clean_matrix(cis,signif=F)
#rownames(B_matrix)=colnames(B_matrix)=name_spp

####Two options are possible for mean values
#If we use raw values of corres_hernandez, we avoid the artefacts created by the interpolation and the random value used when gaps are over 2 points in the time series, but we increase the mean value artificially as cells are counted only when they are numerous. The inverse is true when using interpolated data. This is only a matter of choice.

#abundances_tab=read.table(paste("param/","raw_abundances_Auger.txt",sep=""),sep=";",header=T)
#dates=as.Date(abundances_tab$Date)
#abundances_tab=abundances_tab[year(dates)>=1996,name_spp]#Using data from 1996
#dates=dates[year(dates)>=1996]
#x_obs=apply(abundances_tab,mean,na.rm=T)

pop_table=read.table("../../param/mean_interpolated_abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
pop_table=pop_table[name_spp,]
x_obs=pop_table[name_spp,"Mean_abundances"]
names(x_obs)=name_spp

inter_mat=MAR2BH(B_matrix,x_obs)

#Sinking rates and T_opt
tab=read.table(paste("../../param/species_specific_parameters.txt",sep=""),sep=";",dec=".",header=T)
tab=tab[name_spp,]
S=S_max*tab$S
T_opt=tab$T_opt+273.15+5
names(S)=name_spp
names(T_opt)=name_spp

#First proxy for r
B=tab$Val_b
names(B)=name_spp
r_mean=growth_rate_noMTE_Bissinger(288,T_opt,B,a_d)

write.table(inter_mat,paste("interaction_matrix_before_calibration.csv",sep=""),sep=";",row.names=F,dec=".")


if(quad_prog==1){
	tmp=quadprog(inter_mat,x_obs,r_mean,tol=10000,r_calibrate=F)
	inter_mat=tmp[[1]]
	growth_rate=tmp[[2]]
}else{
	growth_rate=r_mean
}
#write.table(cbind(inter_mat,growth_rate),paste("matrix_A_after_quad.csv",sep=""),sep=";",row.names=F,dec=".")

#############Here, we begin the loop
links_to_explore=which(inter_mat!=0,arr.ind=T)
nb_link=nrow(links_to_explore)

name_links=c()
#Write links, just for the sake of clarity
for(l in 1:nb_link){
	name_links=c(name_links,paste(links_to_explore[l,1],links_to_explore[l,2],sep="-"))
}

tab_simu=matrix(NA,nrow=nb_simu+1,ncol=nb_link)
rownames(tab_simu)=1:(nb_simu+1)
colnames(tab_simu)=name_links

tab_summary=array(NA,dim=c(nb_simu+1,6,nspp),dimnames=list(1:(nb_simu+1),c("Abundance","Amplitude","Phenology","Diff","Persistence_coast","Persistence_ocean"),name_spp))

size_subset=floor(nb_link/(length(value)+1))

id_persistence=seq(n_iter-nb_persistence+1,n_iter)

tab_pheno=read.table("../../param/generalist_specialist_spp_added_amplitude_season.csv",sep=";",dec=".",header=T,row.names=1)

#Before calibration
list_inter=list(inter_mat,k_coast2ocean*inter_mat)
N=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
N[1,,]=rep(10^3,nspp*3)
try(
for(t in 1:(n_iter-1)){
                var_tmp=step1_modelI(N[t,,],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
#                effect_compet[t+1,,]=var_tmp[[2]]
                N[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
,silent=T)

if((any(apply(N[id_persistence,"coast",]==0,2,sum)==nb_persistence))|(any(apply(N[id_persistence,"ocean",]==0,2,sum)==nb_persistence))){
	stop("even in the first simulation we can't have all species")
}

tab_coast=N[,1,]
final_summary_tmp=summary_statistics(pop_table,tab_pheno,tab_coast,nb_year)
final_summary=matrix(unlist(final_summary_tmp),ncol=3)
rownames(final_summary)=names(final_summary_tmp[[1]])
colnames(final_summary)=colnames(tab_summary)[1:3]
write.table(final_summary,"summary_statistics_per_species_before_calibration.txt",dec=".",sep=";")

#Just in case the simulation before calibration is better than everything else
tab_simu[nb_simu+1,]=1.0
tab_summary[nb_simu+1,'Abundance',]=final_summary[[1]]
tab_summary[nb_simu+1,'Amplitude',]=final_summary[[2]]
tab_summary[nb_simu+1,'Phenology',]=final_summary[[3]]
tab_summary[nb_simu+1,"Diff",]=apply(tab_summary[nb_simu+1,1:3,],2,sum)
tab_summary[nb_simu+1,"Persistence_coast",]=apply(N[id_persistence,"coast",]==0,2,sum)<nb_persistence
tab_summary[nb_simu+1,"Persistence_ocean",]=apply(N[id_persistence,"ocean",]==0,2,sum)<nb_persistence


t1=Sys.time()
for(sim in 1:nb_simu){
#for(sim in 1:1){
	print(sim)
	tmp_inter=inter_mat
	set=1:nb_link
	for(v in 1:length(value)){
		inter_v=sample(set,size_subset,replace=F)
		set_tmp=setdiff(set,inter_v)
		for(i in 1:length(inter_v)){
			tmp_inter[links_to_explore[inter_v[i],1],links_to_explore[inter_v[i],2]]=inter_mat[links_to_explore[inter_v[i],1],links_to_explore[inter_v[i],2]]*value[v]
		}
		set=set_tmp
		tab_simu[sim,inter_v]=value[v]
	}
	inter_no_change=set
	tab_simu[sim,inter_no_change]=1.0

	list_inter=list(tmp_inter,k_coast2ocean*tmp_inter)

#Initialize
N=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
#effect_compet=array(NA,dim=c(n_iter,2,nspp),dimnames=list(NULL,c("coast","ocean"),name_spp))
N[1,,]=rep(10^3,nspp*3)

##Run
try(
for(t in 1:(n_iter-1)){
	        var_tmp=step1_modelI(N[t,,],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
		Ntmp=var_tmp[[1]]
#		effect_compet[t+1,,]=var_tmp[[2]]
        	N[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
,silent=T)

tab_coast=N[,1,]
if(sum(is.na(tab_coast)>0)){
	tab_summary[sim,1:3,]=100
}else{
	final_summary=summary_statistics(pop_table,tab_pheno,tab_coast,nb_year)
#	tab_summary[sim,1:3]=c(sqrt(sum(final_summary[[1]])/nspp),sqrt(sum(final_summary[[2]])/nspp),sqrt(sum(final_summary[[3]])/nspp))
	tab_summary[sim,'Abundance',]=final_summary[[1]]
	tab_summary[sim,'Amplitude',]=final_summary[[2]]
	tab_summary[sim,'Phenology',]=final_summary[[3]]
}
#tab_summary[sim,4]=sum(tab_summary[sim,1:3])
tab_summary[sim,"Diff",]=apply(tab_summary[sim,1:3,],2,sum)
tab_summary[sim,"Persistence_coast",]=apply(N[id_persistence,"coast",]==0,2,sum)<nb_persistence
tab_summary[sim,"Persistence_ocean",]=apply(N[id_persistence,"ocean",]==0,2,sum)<nb_persistence

}

print(Sys.time()-t1)
#Write simulations
write.table(tab_simu,paste("list_simulation_calibration.csv",sep=""),sep=";",dec=".")
tab_summary_sum_species=apply(tab_summary[,1:3,],2,translate)
persistence_coast=apply(tab_summary[,"Persistence_coast",],1,sum)
persistence_ocean=apply(tab_summary[,"Persistence_ocean",],1,sum)
diff_sum=apply(tab_summary_sum_species,1,sum)


#best=which(tab_summary[,4]==min(tab_summary[,4],na.rm=T)) Before: we added all RMSE for each interaction model and then the best model was the one which minimized the sum.

#Now we rank each model for each criterion and then choose the one which balances everything
compare_ranks=classify_model(tab_summary[,"Abundance",],tab_summary[,"Amplitude",],tab_summary[,'Phenology',])
#Now, we ensure that a model with at least one missing species is at the end of the ranking
max_rank=max(compare_ranks)+1
id_model_death=unique(c(which(apply(tab_summary[,"Persistence_coast",],1,sum)<nspp),which(apply(tab_summary[,"Persistence_ocean",],1,sum)<nspp)))
compare_ranks[id_model_death]=max_rank
best=which(compare_ranks==min(compare_ranks))[1] #Just in case two models have the same score. We arbitrarily choose the first one.


tab_summary_tmp=cbind(tab_summary_sum_species,diff_sum,persistence_coast,persistence_ocean,compare_ranks)
write.table(tab_summary_tmp,paste("list_statistics_calibration.csv",sep=""),sep=";",dec=".")

best_line_inter=tab_simu[best,]

tmp_inter=inter_mat
for(l in 1:nrow(links_to_explore)){
	tmp_inter[links_to_explore[l,1],links_to_explore[l,2]]=best_line_inter[l]*inter_mat[links_to_explore[l,1],links_to_explore[l,2]]
}

###Finally, we can run the model for the best simulation
list_inter=list(tmp_inter,k_coast2ocean*tmp_inter)

#Initialize
N=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
effect_compet=array(NA,dim=c(n_iter,2,nspp),dimnames=list(NULL,c("coast","ocean"),name_spp))
N[1,,]=rep(10^3,nspp*3)

##Run 
for(t in 1:(n_iter-1)){
                var_tmp=step1_modelI(N[t,,],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                effect_compet[t+1,,]=var_tmp[[2]]
                N[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
tab_coast=N[,1,]
#Statistics per species
final_summary_tmp=summary_statistics(pop_table,tab_pheno,tab_coast,nb_year)
final_summary=matrix(unlist(final_summary_tmp),ncol=3)
rownames(final_summary)=names(final_summary_tmp[[1]])
colnames(final_summary)=colnames(tab_summary)[1:3]
write.table(final_summary,"summary_statistics_per_species_after_calibration.txt",dec=".",sep=";")

#Matrix
write.table(tmp_inter,"interaction_matrix_after_calibration.csv",sep=";",dec=".")

#Dynamics
write.table(N[,1,],paste("out_coast.csv",sep=""),sep=";",dec=".")
write.table(N[,2,],paste("out_ocean.csv",sep=""),sep=";",dec=".")
write.table(N[,3,],paste("out_seed.csv",sep=""),sep=";",dec=".")
#write.table(effect_compet[,1,],paste("compet_coast.csv",sep=""),sep=";",dec=".")
#write.table(effect_compet[,2,],paste("compet_ocean.csv",sep=""),sep=";",dec=".")

source("../../script/diagnostics_single_simu.r")

