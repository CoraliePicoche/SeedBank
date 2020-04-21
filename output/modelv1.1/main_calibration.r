#################
###CP 18/02/2020: Clean version of the main function
###CP 16/04/2020: Version for main calibration
#################

rm(list=ls())
graphics.off()
set.seed(40) #Warning: results and even feasibility can be highly seed-sensitive

source("../../script/matrix_MAR_clean.r")
source("../../script/infer_interaction_matrix_growth_rate.r")
source("../../script/summary_statistics.r")
source("step_functions.r")

####Calibration
value=c(0.9,1.1) #Parameter space
nb_simu=1000
nb_year=2

dataset="Auger"
#Fixed parameters
tab=read.table("simu.csv",sep=";",dec=".",header=T)
cyst_mortality=as.numeric(as.character(tab[tab[,1]=="cyst_mortality",2]))
cyst_burial=as.numeric(as.character(tab[tab[,1]=="cyst_burial",2]))
M=cyst_mortality+cyst_burial
k_coast2ocean=as.numeric(as.character(tab[tab[,1]=="k_coast2ocean",2]))
morta=as.numeric(as.character(tab[tab[,1]=="loss_rate",2]))
e=as.numeric(as.character(tab[tab[,1]=="exchange",2]))
germination=as.numeric(as.character(tab[tab[,1]=="germination",2]))
resuspension=as.numeric(as.character(tab[tab[,1]=="resuspension",2]))
Gamma=resuspension*germination
correct=as.numeric(as.character(tab[tab[,1]=="gain",2]))
growth_model=as.character(tab[tab[,1]=="growth_model",2])
S_max=as.numeric(as.character(tab[tab[,1]=="max_sinking",2]))
quad_prog=tab[tab[,1]=="quad_prog",2]
n_iter=as.numeric(as.character(tab[tab[,1]=="n_iter",2]))
temp_germin=as.numeric(as.character(tab[tab[,1]=="germin_threshold",2]))
a_d=as.numeric(as.character(tab[tab[,1]=="daylength",2]))
threshold=as.numeric(as.character(tab[tab[,1]=="threshold",2]))

#Data to use (Auger)
a=as.matrix(read.table("../../param/reconstructed_temperature_Auger.txt", row.names=1,header=T,sep=";",dec="."))
temp_model=a[1:n_iter]

load(paste("../../param/",dataset,"_pencen_null_regular_common_MO.RData",sep=""))
name_spp=colnames(cis$call$model$B)
B_matrix=clean_matrix(cis,signif=F)
rownames(B_matrix)=colnames(B_matrix)=name_spp
nspp=length(name_spp)

####Two options are possible for mean values
#If we use raw values of corres_hernandez, we avoid the artefacts created by the interpolation and the random value used when gaps are over 2 points in the time series, but we increase the mean value artificially as cells are counted only when they are numerous. The inverse is true when using interpolated data. This is only a matter of choice.

#abundances_tab=read.table(paste("param/","corres_hernandez_Auger.txt",sep=""),sep=";",header=T)
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
tab=read.table(paste("species_specific_parameters.txt",sep=""),sep=";",dec=".",header=T)
tab=tab[name_spp,]
S=S_max*tab$S
T_opt=tab$T_opt+273+5
names(S)=name_spp
names(T_opt)=name_spp

#First proxy for r
B=tab$Val_b
names(B)=name_spp
r_mean=growth_rate_noMTE_Bissinger(288,T_opt,B,a_d)

#write.table(cbind(inter_mat,r_mean),paste("matrix_A_before_quad.csv",sep=""),sep=";",row.names=F,dec=".")


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

tab_simu=matrix(NA,nrow=nb_simu,ncol=nb_link)
rownames(tab_simu)=1:nb_simu
colnames(tab_simu)=name_links

tab_summary=matrix(NA,nrow=nb_simu,ncol=4)
rownames(tab_summary)=1:nb_simu
colnames(tab_summary)=c("Abundance","Amplitude","Phenology","Diff")

size_subset=floor(nb_link/(length(value)+1))


tab_pheno=read.table("../../param/generalist_specialist_spp_added_amplitude_season.csv",sep=";",dec=".",header=T,row.names=1)

t1=Sys.time()
for(sim in 1:nb_simu){
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
	inter_no_change=setdiff(set,nb_link-size_subset*length(value))
	tab_simu[sim,inter_no_change]=1.0

	list_inter=list(tmp_inter,k_coast2ocean*tmp_inter)

#Initialize
N=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
effect_compet=array(NA,dim=c(n_iter,2,nspp),dimnames=list(NULL,c("coast","ocean"),name_spp))
N[1,,]=rep(10^3,nspp*3)

##Run
for(t in 1:(n_iter-1)){
	        var_tmp=step1(N[t,,],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
		Ntmp=var_tmp[[1]]
		effect_compet[t+1,,]=var_tmp[[2]]
        	N[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}

tab_coast=N[,1,]
final_summary=summary_statistics(pop_table,tab_pheno,tab_coast,nb_year)
tab_summary[sim,1:3]=final_summary
tab_summary[sim,4]=sum(final_summary)
}
print(Sys.time()-t1)

#Write simulations
write.table(tab_simu,paste("list_simulation.csv",sep=""),sep=";",dec=".")
write.table(tab_summary,paste("list_statistics.csv",sep=""),sep=";",dec=".")


best=which(tab_summary[,4]==min(tab_summary[,4]))

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
                var_tmp=step1(N[t,,],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                effect_compet[t+1,,]=var_tmp[[2]]
                N[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}

#Matrix
write.table(tmp_inter,"matrix_A_after_calibration.csv",sep=";",dec=".")

#Dynamics
write.table(N[,1,],paste("out_coast.csv",sep=""),sep=";",dec=".")
write.table(N[,2,],paste("out_ocean.csv",sep=""),sep=";",dec=".")
write.table(N[,3,],paste("out_seed.csv",sep=""),sep=";",dec=".")
write.table(effect_compet[,1,],paste("compet_coast.csv",sep=""),sep=";",dec=".")
write.table(effect_compet[,2,],paste("compet_ocean.csv",sep=""),sep=";",dec=".")

source("../../script/diagnostics_single_simu.r")

