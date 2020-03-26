#################
###CP 18/02/2020: Clean version of the main function
#################

rm(list=ls())
graphics.off()
set.seed(40) #Warning: results and even feasibility can be highly seed-sensitive

source("script/matrix_MAR_clean.r")
source("script/infer_interaction_matrix_growth_rate.r")

args = commandArgs(trailingOnly=TRUE)
n_simulation=args[1]

n_simulation=1

#Fixed parameters
tab=read.table(paste("param/simu",n_simulation,".csv",sep=""),sep=";",dec=".",header=T)
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

#Data to use (Auger)
dataset=as.character(tab[tab[,1]=="dataset",2])
evt_tab=read.table(paste("param/",dataset,"hydro.txt",sep=""),sep=";",header=T)
temp=evt_tab[,"TEMP"]
#We need to build a simulation for temperatures as we don't have enough data in the real dataset
min_temp=min(temp,na.rm=T)
max_temp=max(temp,na.rm=T)
sd_temp=sd(temp,na.rm=T)
mean_temp=mean(temp,na.rm=T)
theta=1.3
temp_model=273+mean_temp+theta*sd_temp*sin(2*pi*1:n_iter/365.25)+rnorm(n_iter,0,sd_temp*sqrt(1-theta^2/2))

load(paste("param/",dataset,"_pencen_null_regular_common_MO.RData",sep=""))
name_spp=colnames(cis$call$model$B)
B_matrix=clean_matrix(cis,signif=F)
rownames(B_matrix)=colnames(B_matrix)=name_spp
nspp=length(name_spp)

pop_table=read.table("param/abundances_Auger.txt",sep=",",dec=".",header=T)
name_spp=pop_table$sp
x_obs=pop_table[name_spp,"Mean_abundances"]
names(x_obs)=name_spp

inter_mat=MAR2BH(B_matrix,x_obs)

#Sinking rates and T_opt
tab=read.table(paste("param/species_specific_parameters_",n_simulation,".txt",sep=""),sep=";",dec=".",header=T)
S=S_max*tab$S
T_opt=tab$T_opt+273+5
names(S)=tab$sp
names(T_opt)=tab$sp

#First proxy for r
if(growth_model=="B"){ #B for Bissinger
	r_mean=rep(growth_rate_Bissinger(mean_temp,0.5),nspp)
}else if(growth_model=="SV"){ #SV for Scranton Vasseur, SV_Bissinger for SV model with the Bissinger/noMTE metabolism part
	#Niche area to compute growth rates + range of optimal temperatures
	B=tab$Val_b
	r_mean=growth_rate_SV(293,T_opt,B)
}else if(growth_model=="SV_Bissinger"){ #SV for Scranton Vasseur, SV_Bissinger for SV model with the Bissinger/noMTE metabolism part
	B=tab$Val_b
        r_mean=growth_rate_noMTE_Bissinger(293,T_opt,B)
}else if(growth_model=="fixed"){
	r_mean=0.1
	print("Warning: for now, this option is not implemented")
}else{
	stop("This growth model does not exist")
}

write.table(cbind(inter_mat,r_mean),paste("output/matrix_A_before_quad_",n_simulation,".csv",sep=""),sep=";",row.names=F,dec=".")


if(quad_prog==1){
	if(growth_model=="fixed"){
		tmp=quadprog(inter_mat,x_obs,r_mean,tol=0.1,r_calibrate=T)
	}else{
		tmp=quadprog(inter_mat,x_obs,r_mean,tol=10000,r_calibrate=F)
	}
	inter_mat=tmp[[1]]
	growth_rate=tmp[[2]]
}else{
	growth_rate=r_mean
}
write.table(cbind(inter_mat,growth_rate),paste("output/matrix_A_after_quad_",n_simulation,".csv",sep=""),sep=";",row.names=F,dec=".")

list_inter=list(inter_mat,k_coast2ocean*inter_mat)

#Initialize
N=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
N[1,,]=rep(10^3,nspp*3)


##Run
for(t in 1:(n_iter-1)){
	if(growth_model=="B"){
	        Ntmp=step1(N[t,,],list_inter,temp_model[t],M,morta,model="BH",gr=growth_model,threshold=0.001,irradiance=mean_irr)
	}else if(growth_model=="SV"|growth_model=="SV_Bissinger"){
	        Ntmp=step1(N[t,,],list_inter,temp_model[t],M,morta,correct,model="BH",gr=growth_model,threshold=0.001,T_opt=T_opt,B=B)
	}else if(growth_model=="fixed"){
	        Ntmp=step1(N[t,,],list_inter,temp_model[t],M,morta,model="BH",gr=growth_model,threshold=0.001,fixed_growth=r_mean)
	}else{
		stop("This growth model does not exist")
	}
        N[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}

write.table(N[,1,],paste("output/out_coast",n_simulation,".csv",sep=""),sep=";",dec=".")
write.table(N[,2,],paste("output/out_ocean",n_simulation,".csv",sep=""),sep=";",dec=".")
write.table(N[,3,],paste("output/out_seed",n_simulation,".csv",sep=""),sep=";",dec=".")

#source("script/diagnostics.r")

