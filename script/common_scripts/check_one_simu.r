##########
rm(list=ls())
graphics.off()

source("step_functions.r")

nb_year=2
nb_persistence=6*30

#Fixed parameters: golden parameter set
tab=read.table("../../param/simu_to_explore.csv",sep=";",dec=".",header=T)
cyst_mortality=as.numeric(as.character(tab[tab[,1]=="cyst_mortality",2]))
cyst_burial=as.numeric(as.character(tab[tab[,1]=="cyst_burial",2]))
k_coast2ocean=as.numeric(as.character(tab[tab[,1]=="k_coast2ocean",2]))
morta=as.numeric(as.character(tab[tab[,1]=="loss_rate",2]))
e=as.numeric(as.character(tab[tab[,1]=="exchange",2]))
germination=as.numeric(as.character(tab[tab[,1]=="germination",2]))
resuspension=as.numeric(as.character(tab[tab[,1]=="resuspension",2]))
S_max=as.numeric(as.character(tab[tab[,1]=="max_sinking",2]))
quad_prog=tab[tab[,1]=="quad_prog",2]
temp_germin=as.numeric(as.character(tab[tab[,1]=="germin_threshold",2]))
a_d=as.numeric(as.character(tab[tab[,1]=="daylength",2]))
threshold=as.numeric(as.character(tab[tab[,1]=="threshold",2]))

n_iter=as.numeric(as.character(tab[tab[,1]=="n_iter",2]))
id_persistence=seq(n_iter-nb_persistence+1,n_iter)
id=seq(n_iter-365*nb_year+1,n_iter)


#Data to use (Auger)
a=as.matrix(read.table("../../param/reconstructed_temperature_Auger_with_corrected_phase_and_amplitude.txt", row.names=1,header=T,sep=";",dec="."))
temp_model=a[1:n_iter]
#Used for summary statistics
tab_pheno=read.table("../../param/generalist_specialist_spp_added_amplitude_season.csv",sep=";",dec=".",header=T,row.names=1)
name_spp=rownames(tab_pheno)
nspp=length(name_spp)
####Two options are possible for mean values
pop_table=read.table("../../param/mean_interpolated_abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
pop_table=pop_table[name_spp,]

#Sinking rates and T_opt
tab=read.table(paste("../../param/species_specific_parameters.txt",sep=""),sep=";",dec=".",header=T)
tab=tab[name_spp,]
S=S_max*tab$S
T_opt=tab$T_opt+273.15+5
names(S)=name_spp
names(T_opt)=name_spp
B=tab$Val_b
names(B)=name_spp

#Composite phenomena
M=cyst_mortality+cyst_burial
Gamma=resuspension*germination
S=S_max*tab$S

####Model I
#Interaction
tmp_inter_I=read.table("../../output/modelv1.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
tmp_inter_I=as.matrix(tmp_inter_I)
list_inter_I=list(tmp_inter_I,k_coast2ocean*tmp_inter_I)

###############Golden_set
#Initialize
N_original_set=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
N_original_set[1,,]=10^3

##Run
tt=c()
where=c()
who=c()
for(t in 1:(n_iter-1)){
                var_tmp=step1_modelI(N_original_set[t,,],list_inter_I,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_original_set[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
		if(sum(N_original_set[t+1,,]==0)>0){
			a=which(N_original_set[t+1,,]==0,arr.ind=T)
			tt=c(tt,t)
			where=c(where,dimnames(N_original_set)[[2]][a[,1]])
			who=c(who,name_spp[a[,2]])
		}
}
##Every year, DIT disappears and reappears
