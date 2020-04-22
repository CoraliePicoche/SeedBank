#################
###CP 18/02/2020: Clean version of the main function 
###CP 08/04/2020: Main function for the saturating interaction model
###CP 16/04/20202: Added calibration on H_ij
#################

rm(list=ls())
graphics.off()
set.seed(40) #Warning: results and even feasibility can be highly seed-sensitive

source("../../script/matrix_MAR_clean.r")
source("../../script/summary_statistics.r")
source("compute_saturating_interactions.r")

####Calibration
value=c(0.5,0.9,1.1,2) #Parameter space
nb_simu=100
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
#threshold=as.numeric(as.character(tab[tab[,1]=="threshold",2]))
O_y=as.numeric(as.character(tab[tab[,1]=="overyielding",2]))
ratio_pos=as.numeric(as.character(tab[tab[,1]=="ratio_pos",2]))

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

abundances_tab=read.table(paste("../../param/corres_hernandez_Auger.txt",sep=""),sep=";",header=T)
#abundances_tab=read.table(paste("param/corres_hernandez_Auger.txt",sep=""),sep=";",header=T)
dates=as.Date(abundances_tab$Date)
abundances_tab=abundances_tab[year(dates)>=1996,name_spp]#Using data from 1996
N_max=apply(abundances_tab,2,max,na.rm=T)
#N_mean=apply(abundances_tab,2,mean,na.rm=T)

pop_table=read.table("../../param/mean_interpolated_abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
pop_table=pop_table[name_spp,]
N_mean=pop_table[name_spp,"Mean_abundances"]
names(N_mean)=name_spp


inter_mat=MAR2saturation(B_matrix,N_mean,N_max,O_y,ratio_pos)
print("For now, interactions on the coast and in the ocean are the same")
type_inter=list(inter_mat[[1]],inter_mat[[1]])
list_H=list(inter_mat[[2]],inter_mat[[2]])

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

print("No quadratic programming available for saturating interactions, for now")
write.table(inter_mat[[2]],paste("matrix_A_before_calibration.csv",sep=""),sep=";",row.names=F,dec=".")

#list_inter=list(inter_mat,k_coast2ocean*inter_mat)

#############Here, we begin the loop
links_to_explore=which(inter_mat[[2]]!=0,arr.ind=T)
nb_link=nrow(links_to_explore)

name_links=c()
#Write links, just for the sake of clarity
for(l in 1:nb_link){
        name_links=c(name_links,paste(links_to_explore[l,1],links_to_explore[l,2],sep="-"))
}

tab_simu=matrix(NA,nrow=nb_simu,ncol=nb_link)
rownames(tab_simu)=1:nb_simu
colnames(tab_simu)=name_links

tab_summary=matrix(NA,nrow=nb_simu,ncol=5)
rownames(tab_summary)=1:nb_simu
colnames(tab_summary)=c("Abundance","Amplitude","Phenology","Diff","Persistence")

size_subset=floor(nb_link/(length(value)+1))

pop_table=as.data.frame(cbind(names(N_mean),N_mean))
pop_table[,2]=as.numeric(as.character(N_mean))
tab_pheno=read.table("../../param/generalist_specialist_spp_added_amplitude_season.csv",sep=";",dec=".",header=T,row.names=1)

t1=Sys.time()
for(sim in 1:nb_simu){
        print(sim)
        tmp_inter=list_H[[1]]
        set=1:nb_link
        for(v in 1:length(value)){
                inter_v=sample(set,size_subset,replace=F)
                set_tmp=setdiff(set,inter_v)
                for(i in 1:length(inter_v)){
                        tmp_inter[links_to_explore[inter_v[i],1],links_to_explore[inter_v[i],2]]=list_H[[1]][links_to_explore[inter_v[i],1],links_to_explore[inter_v[i],2]]*value[v]
                }
                set=set_tmp
                tab_simu[sim,inter_v]=value[v]
        }
        inter_no_change=set
        tab_simu[sim,inter_no_change]=1.0

        list_H_tmp=list(tmp_inter,tmp_inter)


#Initialize
N=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
effect_compet=array(NA,dim=c(n_iter,2,nspp),dimnames=list(NULL,c("coast","ocean"),name_spp))
N[1,,]=rep(10^3,nspp*3)

##Run
try(
for(t in 1:(n_iter-1)){
		var_tmp=step1(N[t,,],list_H_tmp,type_inter,temp_model[t],M,morta,a_d,T_opt,B)
		Ntmp=var_tmp[[1]]
		effect_compet[t+1,,]=var_tmp[[2]]
        	N[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
,silent=T)
tab_coast=N[,1,]
if(sum(is.na(tab_coast)>0)){
tab_summary[sim,1:3]=100
}else{
final_summary=summary_statistics(pop_table,tab_pheno,tab_coast,nb_year)
tab_summary[sim,1:3]=c(sum(final_summary[[1]]),sum(final_summary[[2]]),sum(final_summary[[3]]))
}
tab_summary[sim,4]=sum(tab_summary[sim,1:3])
tab_summary[sim,5]=sum(tab_coast>0,na.rm=T)
if(tab_summary[sim,5]<ncol(tab_coast)){tab_summary[,4]=tab_summary[,4]*2} #We can't have a model with a missing species

}

print(Sys.time()-t1)

#Write simulations
write.table(tab_simu,paste("list_simulation.csv",sep=""),sep=";",dec=".")
write.table(tab_summary,paste("list_statistics.csv",sep=""),sep=";",dec=".")


best=which(tab_summary[,4]==min(tab_summary[,4],na.rm=T))

best_line_inter=tab_simu[best,]

tmp_inter=list_H[[1]]
for(l in 1:nrow(links_to_explore)){
        tmp_inter[links_to_explore[l,1],links_to_explore[l,2]]=best_line_inter[l]*list_H[[1]][links_to_explore[l,1],links_to_explore[l,2]]
}

###Finally, we can run the model for the best simulation
        list_H_tmp=list(tmp_inter,tmp_inter)

#Initialize
N=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
effect_compet=array(NA,dim=c(n_iter,2,nspp),dimnames=list(NULL,c("coast","ocean"),name_spp))
N[1,,]=rep(10^3,nspp*3)

##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1(N[t,,],list_H_tmp,type_inter,temp_model[t],M,morta,a_d,T_opt,B)
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
write.table(final_summary,"summary_statistics_per_species.txt",dec=".",sep=";")



colnames(list_H_tmp[[1]])=name_spp
rownames(list_H_tmp[[1]])=name_spp
colnames(type_inter[[1]])=name_spp
rownames(type_inter[[1]])=name_spp
write.table(list_H_tmp[[1]],"matrix_A_after_calibration.csv",sep=";",dec=".")
write.table(type_inter[[1]],"facilitation_and_competition.csv",sep=";",dec=".")


write.table(N[,1,],paste("out_coast.csv",sep=""),sep=";",dec=".")
write.table(N[,2,],paste("out_ocean.csv",sep=""),sep=";",dec=".")
write.table(N[,3,],paste("out_seed.csv",sep=""),sep=";",dec=".")
write.table(effect_compet[,1,],paste("compet_coast.csv",sep=""),sep=";",dec=".")
write.table(effect_compet[,2,],paste("compet_ocean.csv",sep=""),sep=";",dec=".")

source("../../script/diagnostics_single_simu.r")

