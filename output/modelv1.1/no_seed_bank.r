#####
## 16/05/2020 CP: Simulations without seed bank. 
#####

rm(list=ls())
graphics.off()

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

#Interaction
tmp_inter=read.table("interaction_matrix_after_calibration.csv",sep=";",dec=".")
name_spp=rownames(tmp_inter)
nspp=length(name_spp)
list_inter=list(as.matrix(tmp_inter),k_coast2ocean*as.matrix(tmp_inter)) #Interaction list contains coastal and oceanic interactions

#Data to use (Auger)
a=as.matrix(read.table("../../param/reconstructed_temperature_Auger.txt", row.names=1,header=T,sep=";",dec="."))
temp_model=a[1:n_iter]
####Two options are possible for mean values
pop_table=read.table("../../param/mean_interpolated_abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
pop_table=pop_table[name_spp,]

#Sinking rates and T_opt
tab=read.table(paste("../../param/species_specific_parameters.txt",sep=""),sep=";",dec=".",header=T)
tab=tab[name_spp,]
S=S_max*tab$S
T_opt=tab$T_opt+273+5
names(S)=name_spp
names(T_opt)=name_spp
B=tab$Val_b
names(B)=name_spp

#Composite phenomena
M=cyst_mortality+cyst_burial
Gamma=resuspension*germination
S=S_max*tab$S


N_original_test_coast=read.table("out_coast.csv",sep=";",row.names=1,header=T)
N_original_test_ocean=read.table("out_ocean.csv",sep=";",row.names=1,header=T)
n_iter=nrow(N_original_test_coast)
id=(n_iter-365):n_iter
####First test, only removing seed bank
M=1
#Initialize
N_no_seed=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
#effect_compet=array(NA,dim=c(n_iter,2,nspp),dimnames=list(NULL,c("coast","ocean"),name_spp))
N_no_seed[1,,]=rep(10^3,nspp*3)


##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1(N_no_seed[t,,],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_no_seed[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}

pdf("seed_bank.pdf",width=10)
for(s in name_spp){
print(s)

if(N_no_seed[dim(N_no_seed)[1],"coast",s]>0){
	id_tmp=id
}else{
	first_0_coast=which(N_no_seed[,"coast",s]==0)[1]
	first_0_ocean=which(N_no_seed[,"ocean",s]==0)[1]
	first_0=max(c(first_0_coast,first_0_ocean))
	id_tmp=seq(1,first_0)
}


par(mfrow=c(1,2))
plot(id_tmp,log10(N_original_test_coast[id_tmp,s]+10^(-5)),col="cyan",pch=16,t="o",lty=1,main=s)
lines(id_tmp,log10(N_original_test_ocean[id_tmp,s]+10^(-5)),col="darkblue")

plot(id_tmp,log10(N_no_seed[id_tmp,"coast",s]+10^(-5)),col="cyan",pch=16,t="o",lty=1,main=s)
lines(id_tmp,log10(N_no_seed[id_tmp,"ocean",s]+10^(-5)),col="darkblue")
}
dev.off()


##No, what happens if I lower competition
tmp_inter_lower=tmp_inter
tmp_inter_lower[which(tmp_inter>0,arr.ind=T)]=tmp_inter[which(tmp_inter>0,arr.ind=T)]*0.1
list_inter_lower=list(as.matrix(tmp_inter_lower),as.matrix(tmp_inter_lower))
#Initialize
N_no_seed_lower_compet=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
#effect_compet=array(NA,dim=c(n_iter,2,nspp),dimnames=list(NULL,c("coast","ocean"),name_spp))
N_no_seed_lower_compet[1,,]=rep(10^3,nspp*3)

print("Lower")
##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1(N_no_seed_lower_compet[t,,],list_inter_lower,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_no_seed_lower_compet[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
print(sum(N_no_seed_lower_compet[1000,"ocean",]>0))

##No, what happens if I increase competition
tmp_inter_higher=tmp_inter
tmp_inter_higher[which(tmp_inter>0,arr.ind=T)]=tmp_inter[which(tmp_inter>0,arr.ind=T)]*5
list_inter_higher=list(as.matrix(tmp_inter_higher),as.matrix(tmp_inter_higher))
#Initialize
N_no_seed_higher_compet=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
#effect_compet=array(NA,dim=c(n_iter,2,nspp),dimnames=list(NULL,c("coast","ocean"),name_spp))
N_no_seed_higher_compet[1,,]=rep(10^3,nspp*3)

print("higher")
##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1(N_no_seed_higher_compet[t,,],list_inter_higher,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_no_seed_higher_compet[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
print(sum(N_no_seed_higher_compet[1000,"ocean",]>0))


