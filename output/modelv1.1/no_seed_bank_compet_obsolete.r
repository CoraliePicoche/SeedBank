2#####
## 16/05/2020 CP: Simulations without seed bank. 
#####

rm(list=ls())
graphics.off()
set.seed(42)

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
tmp_inter=as.matrix(tmp_inter)
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
mean_val=log10(apply(N_original_test_coast[id,],2,mean))
min_val=log10(apply(N_original_test_coast[id,],2,min))
range_val=log10(apply(apply(N_original_test_coast[id,],2,range),2,diff))
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


pdf("no_seed_bank_compet.pdf")
par(mfrow=c(2,2))
surviv=rep(0,nspp)
names(surviv)=name_spp

nb_simu=100
lower_simu=seq(0.1,1,length.out=nb_simu)
higher_simu=seq(1,10,length.out=nb_simu)

N_no_seed_lower_compet=array(NA,dim=c(n_iter,3,nspp,nb_simu,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(lower_simu,digits=1),c("Inc. self-reg","Not inc. self-reg")))
N_no_seed_lower_compet[1,,,,]=rep(10^3,nspp*3*nb_simu*2)
N_no_seed_higher_compet=array(NA,dim=c(n_iter,3,nspp,nb_simu,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(higher_simu,digits=1),c("Inc. self-reg","Not inc. self-reg")))
N_no_seed_higher_compet[1,,,,]=rep(10^3,nspp*3*nb_simu*2)
N_no_seed_lower_mutual=array(NA,dim=c(n_iter,3,nspp,nb_simu),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(lower_simu,digits=1)))
N_no_seed_lower_mutual[1,,,]=rep(10^3,nspp*3*nb_simu)
N_no_seed_higher_mutual=array(NA,dim=c(n_iter,3,nspp,nb_simu),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(higher_simu,digits=1)))
N_no_seed_higher_mutual[1,,,]=rep(10^3,nspp*3*nb_simu)

##Now, what happens if I lower competition?
richness_lower=array(NA,dim=c(length(lower_simu),2,2),dimnames=list(format(lower_simu,digits=1),c("Inc. self-reg","Not inc. self-reg"),c("Coast","Ocean")))
plot(0,0,t="n",ylim=c(0,11),xlim=range(lower_simu))
for(sr in 1:2){
for(fac in 1:length(lower_simu)){
#Changing all competititive interactions, including density-dependence
tmp_inter_lower=tmp_inter
if(sr==2){
	diag(tmp_inter_lower)=NA
}
tmp_inter_lower[which(tmp_inter>0,arr.ind=T)]=tmp_inter[which(tmp_inter>0,arr.ind=T)]*lower_simu[fac]
if(sr==2){
	diag(tmp_inter_lower)=diag(tmp_inter)
}
list_inter_lower=list(as.matrix(tmp_inter_lower),as.matrix(tmp_inter_lower))

#Initialize
print("Lower")
##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1(N_no_seed_lower_compet[t,,,fac,sr],list_inter_lower,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_no_seed_lower_compet[t+1,,,fac,sr]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
#print(sum(N_no_seed_lower_compet[1000,"ocean",]>0))
richness_lower[fac,sr,"Coast"]=sum(N_no_seed_lower_compet[1000,"coast",,fac,sr]>0)
richness_lower[fac,sr,"Ocean"]=sum(N_no_seed_lower_compet[1000,"ocean",,fac,sr]>0)
surviv=surviv+(N_no_seed_lower_compet[1000,"ocean",,fac,sr]>0)
}
if(sr==1){
	apch=1
}else{
	apch=16
}
points(lower_simu,richness_lower[,sr,"Coast"],col="cyan",pch=apch)
points(lower_simu,richness_lower[,sr,"Ocean"],col="darkblue",pch=apch,cex=0.5)
}

##Now, what happens if I increase competition
richness_higher=array(NA,dim=c(length(higher_simu),2,2),dimnames=list(format(higher_simu,digits=1),c("Inc. self-reg","Not inc. self-reg"),c("Coast","Ocean")))
plot(0,0,t="n",xlim=range(higher_simu),ylim=c(0,11))
for(sr in 1:2){
for(fac in 1:length(higher_simu)){
tmp_inter_higher=tmp_inter
if(sr==2){
	diag(tmp_inter_lower)=NA
}
tmp_inter_higher[which(tmp_inter>0,arr.ind=T)]=tmp_inter[which(tmp_inter>0,arr.ind=T)]*higher_simu[fac]
if(sr==2){
	diag(tmp_inter_higher)=diag(tmp_inter)
}
list_inter_higher=list(as.matrix(tmp_inter_higher),as.matrix(tmp_inter_higher))

print("higher")
##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1(N_no_seed_higher_compet[t,,,fac,sr],list_inter_higher,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_no_seed_higher_compet[t+1,,,fac,sr]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
#print(sum(N_no_seed_higher_compet[1000,"ocean",]>0))
richness_higher[fac,sr,"Coast"]=sum(N_no_seed_higher_compet[1000,"coast",,fac,sr]>0)
richness_higher[fac,sr,"Ocean"]=sum(N_no_seed_higher_compet[1000,"ocean",,fac,sr]>0)
}
if(sr==1){
	apch=1
}else{
	apch=16
}
points(higher_simu,richness_higher[,sr,"Coast"],col="cyan",pch=apch)
points(higher_simu,richness_higher[,sr,"Ocean"],col="darkblue",pch=apch,cex=0.8)
surviv=surviv+(N_no_seed_higher_compet[1000,"ocean",,fac,sr]>0)
}

#And what happens if I decrease mutualism?
richness_lower=matrix(NA,length(lower_simu),2,dimnames=list(format(lower_simu,digits=1),c("Coast","Ocean")))
for(fac in 1:length(lower_simu)){
tmp_inter_lower=tmp_inter
tmp_inter_lower[which(tmp_inter<0,arr.ind=T)]=tmp_inter[which(tmp_inter<0,arr.ind=T)]*lower_simu[fac]
list_inter_lower=list(as.matrix(tmp_inter_lower),as.matrix(tmp_inter_lower))

print("Lower")
##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1(N_no_seed_lower_mutual[t,,,fac],list_inter_lower,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_no_seed_lower_mutual[t+1,,,fac]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
#print(sum(N_no_seed_lower_compet[1000,"ocean",]>0))
richness_lower[fac,"Coast"]=sum(N_no_seed_lower_mutual[1000,"coast",,fac]>0)
richness_lower[fac,"Ocean"]=sum(N_no_seed_lower_mutual[1000,"ocean",,fac]>0)
surviv=surviv+(N_no_seed_lower_mutual[1000,"ocean",,fac]>0)
}
plot(lower_simu,richness_lower[,"Coast"],col="cyan",pch=16)
points(lower_simu,richness_lower[,"Ocean"],col="darkblue",pch=16)

##Now, what happens if I increase mutualism?
richness_higher=matrix(NA,length(higher_simu),2,dimnames=list(format(higher_simu,digits=1),c("Coast","Ocean")))
for(fac in 1:length(higher_simu)){
tmp_inter_higher=tmp_inter
tmp_inter_higher[which(tmp_inter<0,arr.ind=T)]=tmp_inter[which(tmp_inter<0,arr.ind=T)]*higher_simu[fac]
list_inter_higher=list(as.matrix(tmp_inter_higher),as.matrix(tmp_inter_higher))

print("higher")
##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1(N_no_seed_higher_mutual[t,,,fac],list_inter_higher,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_no_seed_higher_mutual[t+1,,,fac]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
#print(sum(N_no_seed_higher_compet[1000,"ocean",]>0))
richness_higher[fac,"Coast"]=sum(N_no_seed_higher_mutual[1000,"coast",,fac]>0)
richness_higher[fac,"Ocean"]=sum(N_no_seed_higher_mutual[1000,"ocean",,fac]>0)
surviv=surviv+(N_no_seed_higher_mutual[1000,"ocean",,fac]>0)
}

plot(higher_simu,richness_higher[,"Coast"],col="cyan",pch=16)
points(higher_simu,richness_higher[,"Ocean"],col="darkblue",pch=16)
dev.off()

pdf("species_no_seed_bank.pdf",width=10)
par(mfrow=c(3,4))

for(s in 1:nspp){
	plot(1:nb_simu,apply(log10(N_no_seed_lower_compet[id,"ocean",s,,1]+10^-5),2,mean),pch=1,ylim=c(-2,7),main=name_spp[s],col="lightblue",cex=2)
	points(1:nb_simu,apply(log10(N_no_seed_lower_compet[id,"ocean",s,,2]+10^-5),2,mean),pch=16,col="lightblue",cex=2)
	points(1:nb_simu,apply(log10(N_no_seed_higher_compet[id,"ocean",s,,1]+10^-5),2,mean),pch=1,col="darkblue",cex=2)
	points(1:nb_simu,apply(log10(N_no_seed_higher_compet[id,"ocean",s,,2]+10^-5),2,mean),pch=16,col="darkblue",cex=2)
	points(1:nb_simu,apply(log10(N_no_seed_lower_mutual[id,"ocean",s,]+10^-5),2,mean),pch=16,col="orange",cex=2)
	points(1:nb_simu,apply(log10(N_no_seed_higher_mutual[id,"ocean",s,]+10^-5),2,mean),pch=16,col="darkred",cex=2)
	abline(h=log10(mean_val[s]))
}
dev.off()

pdf("dynamics_vs_survival_no_seed_bank.pdf")
par(mfrow=c(1,3))
plot(log10(surviv+0.1),mean_val)
text(log10(surviv+0.1),mean_val,name_spp,pos=2)

plot(log10(surviv+0.1),range_val)
text(log10(surviv+0.1),range_val,name_spp,pos=2)

plot(log10(surviv+0.1),min_val)
text(log10(surviv+0.1),min_val,name_spp,pos=2)
dev.off()
