#################
###CP 18/02/2020: Clean version of the main function
#################

rm(list=ls())
graphics.off()
set.seed(40) #Warning: results and even feasibility can be highly seed-sensitive

source("../../script/matrix_MAR_clean.r")
source("../../script/infer_interaction_matrix_growth_rate.r")
source("step_functions.r")

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
a=as.matrix(read.table("param/reconstructed_temperature_Auger.txt", row.names=1,header=T,sep=";",dec="."))
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

write.table(cbind(inter_mat,r_mean),paste("matrix_A_before_quad.csv",sep=""),sep=";",row.names=F,dec=".")


if(quad_prog==1){
	tmp=quadprog(inter_mat,x_obs,r_mean,tol=10000,r_calibrate=F)
	inter_mat=tmp[[1]]
	growth_rate=tmp[[2]]
}else{
	growth_rate=r_mean
}
write.table(cbind(inter_mat,growth_rate),paste("matrix_A_after_quad.csv",sep=""),sep=";",row.names=F,dec=".")

inter_ocean=inter_mat
inter_ocean[inter_mat<0]=0
list_inter=list(inter_mat,k_coast2ocean*inter_ocean)

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

write.table(N[,1,],paste("out_coast.csv",sep=""),sep=";",dec=".")
write.table(N[,2,],paste("out_ocean.csv",sep=""),sep=";",dec=".")
write.table(N[,3,],paste("out_seed.csv",sep=""),sep=";",dec=".")
write.table(effect_compet[,1,],paste("compet_coast.csv",sep=""),sep=";",dec=".")
write.table(effect_compet[,2,],paste("compet_ocean.csv",sep=""),sep=";",dec=".")

source("../../script/diagnostics_single_simu.r")

