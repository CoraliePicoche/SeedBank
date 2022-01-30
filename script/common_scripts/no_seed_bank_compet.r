####
# 18/05/2020 CP: Effect of competition on communities with and without seed bank
# 21/05/2020 CP: Wrote the code so that it could be applied to both models at the same time
# 25/02/2021 CP: Same code, but adding the cyst bank in the figure
# 01/03/2021 CP: Bray Curtis seed bank vs ocean coast
# 25/01/2022 CP: Added a threshold on the denominator on model II
###

rm(list=ls())
graphics.off()
set.seed(42)

source("step_functions.r")
source("Bray_Curtis.r")
nb_year=2
cpt="ocean"
nb_persistence=6*30

doyouload=TRUE #If TRUE, load previous results to plot them ; if FALSE, relaunch simulations

#Fixed parameters: golden parameter set
tab=read.table("../../param/simu.csv",sep=";",dec=".",header=T)
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
M_orig=M
Gamma=resuspension*germination
S=S_max*tab$S

####Model I
#Interaction
tmp_inter_I=read.table("../../output/modelv1.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
tmp_inter_I=as.matrix(tmp_inter_I)
list_inter_I=list(tmp_inter_I,k_coast2ocean*tmp_inter_I) #Interaction list contains coastal and oceanic interactions

####ModelII
list_H_tmp=read.table("../../output/modelv2.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
list_H=list(list_H_tmp,list_H_tmp)
type_inter_tmp=read.table("../../output/modelv2.1/facilitation_and_competition.csv",sep=";",dec=".")
type_inter_II=list(type_inter_tmp,k_coast2ocean*type_inter_tmp)

if(!doyouload){

###############Golden_set
#Initialize
N_original_set=array(NA,dim=c(n_iter,3,nspp,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,c("modelI","modelII")))
N_original_set[1,,,]=10^3

##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1_modelI(N_original_set[t,,,"modelI"],list_inter_I,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_original_set[t+1,,,'modelI']=step2(Ntmp,S,Gamma,e)

                var_tmp=step1_modelII(N_original_set[t,,,'modelII'],list_H,type_inter_II,temp_model[t],M,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
                N_original_set[t+1,,,'modelII']=step2(Ntmp,S,Gamma,e)
}

mean_val=log10(apply(N_original_set[id,"ocean",,],c(2,3),mean))
min_val=log10(apply(N_original_set[id,"ocean",,],c(2,3),min))
range_val=apply(log10(apply(N_original_set[id,"ocean",,],c(2,3),range)),c(2,3),diff)


######################Param changed when using seed_bank
M_orig=M
M=1
nb_div=10
simu=seq(1,10,length.out=nb_div)
fac_simu=unique(c(rev(1/simu),1,simu))


###Initialize
N_array=array(NA,dim=c(length(id),3,nspp,length(fac_simu),2,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(fac_simu,digits=3),c("compet","facil"),c("ModelI","ModelII")))
N_orig=array(NA,dim=c(length(id),3,nspp,length(fac_simu),2,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(fac_simu,digits=3),c("compet","facil"),c("ModelI","ModelII")))

gr_array=array(NA,dim=c(length(id),2,nspp,length(fac_simu),2,2),dimnames=list(NULL,c("coast","ocean"),name_spp,format(fac_simu,digits=3),c("compet","facil"),c("ModelI","ModelII")))
gr_orig=array(NA,dim=c(length(id),2,nspp,length(fac_simu),2,2),dimnames=list(NULL,c("coast","ocean"),name_spp,format(fac_simu,digits=3),c("compet","facil"),c("ModelI","ModelII")))

for(fac in 1:length(fac_simu)){
	print(fac_simu[fac])
	####ModelI
	#Compet
	tmp_inter_compet=tmp_inter_I
	tmp_inter_compet[which(tmp_inter_I>0,arr.ind=T)]=tmp_inter_I[which(tmp_inter_I>0,arr.ind=T)]*fac_simu[fac]
	list_inter_compet=list(tmp_inter_compet,tmp_inter_compet)

	#Facilitation
	tmp_inter_facil=tmp_inter_I
	tmp_inter_facil[which(tmp_inter_I<0,arr.ind=T)]=tmp_inter_I[which(tmp_inter_I<0,arr.ind=T)]*fac_simu[fac]
	list_inter_facil=list(tmp_inter_facil,tmp_inter_facil)


	###ModelII
	tmp_inter_compet=type_inter_tmp
        tmp_inter_compet[which(type_inter_tmp>0,arr.ind=T)]=type_inter_tmp[which(type_inter_tmp>0,arr.ind=T)]*fac_simu[fac]
        type_inter_compet=list(tmp_inter_compet,tmp_inter_compet)

        #Facilitation
        tmp_inter_facil=type_inter_tmp
        tmp_inter_facil[which(type_inter_tmp<0,arr.ind=T)]=type_inter_tmp[which(type_inter_tmp<0,arr.ind=T)]*fac_simu[fac]
        type_inter_facil=list(tmp_inter_facil,tmp_inter_facil)


	N_array_simu=array(NA,dim=c(n_iter,3,nspp,2,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,c("compet","facil"),c("ModelI","ModelII")))
	gr_array_simu=array(NA,dim=c(n_iter,2,nspp,2,2),dimnames=list(NULL,c("coast","ocean"),name_spp,c("compet","facil"),c("ModelI","ModelII")))
	N_array_simu[1,,,,]=10^3
	N_orig_simu=array(NA,dim=c(n_iter,3,nspp,2,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,c("compet","facil"),c("ModelI","ModelII")))
	gr_orig_simu=array(NA,dim=c(n_iter,2,nspp,2,2),dimnames=list(NULL,c("coast","ocean"),name_spp,c("compet","facil"),c("ModelI","ModelII")))
	N_orig_simu[1,,,,]=10^3
	for(t in 1:(n_iter-1)){
		#ModelI		
		#Without Seed Bank
		#Compet
                var_tmp=step1_modelI(N_array_simu[t,,,"compet","ModelI"],list_inter_compet,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
		gr_array_simu[t+1,,,"compet","ModelI"]=var_tmp[[2]]
                N_array_simu[t+1,,,"compet","ModelI"]=step2(Ntmp,S,Gamma,e)
                
		#Facil
		var_tmp=step1_modelI(N_array_simu[t,,,"facil","ModelI"],list_inter_facil,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
		gr_array_simu[t+1,,,"facil","ModelI"]=var_tmp[[2]]
                N_array_simu[t+1,,,"facil","ModelI"]=step2(Ntmp,S,Gamma,e)
		
		#With Seed Bank
		#Compet
                var_tmp=step1_modelI(N_orig_simu[t,,,"compet","ModelI"],list_inter_compet,temp_model[t],M_orig,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
		gr_orig_simu[t+1,,,"compet","ModelI"]=var_tmp[[2]]
                N_orig_simu[t+1,,,"compet","ModelI"]=step2(Ntmp,S,Gamma,e)
                
		#Facil
		var_tmp=step1_modelI(N_orig_simu[t,,,"facil","ModelI"],list_inter_facil,temp_model[t],M_orig,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
		gr_orig_simu[t+1,,,"facil","ModelI"]=var_tmp[[2]]
                N_orig_simu[t+1,,,"facil","ModelI"]=step2(Ntmp,S,Gamma,e)

		##ModelII
                #Without Seed Bank
                #Compet
                var_tmp=step1_modelII(N_array_simu[t,,,"compet","ModelII"],list_H,type_inter_compet,temp_model[t],M,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
		gr_array_simu[t+1,,,"compet","ModelII"]=var_tmp[[2]]
                N_array_simu[t+1,,,"compet","ModelII"]=step2(Ntmp,S,Gamma,e)

                #Facil
                var_tmp=step1_modelII(N_array_simu[t,,,"facil","ModelII"],list_H,type_inter_facil,temp_model[t],M,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
		gr_array_simu[t+1,,,"facil","ModelII"]=var_tmp[[2]]
                N_array_simu[t+1,,,"facil","ModelII"]=step2(Ntmp,S,Gamma,e)

                #With Seed Bank
                #Compet
                var_tmp=step1_modelII(N_orig_simu[t,,,"compet","ModelII"],list_H,type_inter_compet,temp_model[t],M_orig,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
		gr_orig_simu[t+1,,,"compet","ModelII"]=var_tmp[[2]]
                N_orig_simu[t+1,,,"compet","ModelII"]=step2(Ntmp,S,Gamma,e)

                #Facil
                var_tmp=step1_modelII(N_orig_simu[t,,,"facil","ModelII"],list_H,type_inter_facil,temp_model[t],M_orig,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
		gr_orig_simu[t+1,,,"facil","ModelII"]=var_tmp[[2]]
                N_orig_simu[t+1,,,"facil","ModelII"]=step2(Ntmp,S,Gamma,e)



	}
	N_array[,,,fac,,]=N_array_simu[id,,,,]
	N_orig[,,,fac,,]=N_orig_simu[id,,,,]
	gr_array[,,,fac,,]=gr_array_simu[id,,,,]
	gr_orig[,,,fac,,]=gr_orig_simu[id,,,,]
}

save(N_array, N_orig, mean_val,min_val,range_val,fac_simu,nb_persistence,gr_array,gr_orig,id,file = "no_seed_compet_simu_light_threshold_modelII_whole_range.RData", envir = .GlobalEnv)
}else{ #end !nodoyouload
load("no_seed_compet_simu_light_threshold_modelII_whole_range.RData")
}

id_persistence=seq(length(id)-nb_persistence+1,length(id)) 

pdf("no_seed_bank_compet.pdf")
par(mfrow=c(2,2),mar=c(1.75,4.5,3,0))
 #Final richness, taking into account the last months, instead of only the end of the simulation
plot(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'compet','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),ylim=c(0,11.5),ylab="Richness",xlab="",main="Competition",log="x",col="black",lty=1,t="l",xaxt="n",cex.axis=1.5,cex.lab=1.5,lwd=2,cex.main=1.75) 
mtext("a",2,line=-0.6,at=12.6,cex=1.25,las=2,font=2)
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'compet','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),col="black",lty=2,lwd=2)
lines(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'compet','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="grey",lty=1,lwd=2)
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'compet','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="grey",lty=2,lwd=2)

par(mar=c(1.75,2.5,3,2))
plot(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'facil','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),ylim=c(0,11.5),ylab="",xlab="",main="Facilitation",log="x",col="black",lty=1,t="l",xaxt="n",yaxt="n",lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.75)
mtext("b",2,line=-0.6,at=12.6,cex=1.25,las=2,font=2)
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'facil','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),col="black",lty=2,lwd=2)
lines(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'facil','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="grey",lty=1,lwd=2)
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'facil','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="grey",lty=2,lwd=2)

par(mar=c(4,4.5,1.25,0))

N_tot=matrix(NA,length(fac_simu),2)
N_tot_orig=matrix(NA,length(fac_simu),2)

#SIMPSON
simpson=matrix(NA,length(fac_simu),2)
simpson_orig=matrix(NA,length(fac_simu),2)
for(f in 1:length(fac_simu)){
for(m in 1:2){
	abundance_tot=mean(apply(N_array[,'ocean',,f,'compet',m],1,sum,na.rm=T))
	abundance_per_species=apply(N_array[,'ocean',,f,'compet',m],2,mean,na.rm=T)
	proportion=abundance_per_species/abundance_tot
	simpson[f,m]=sum(proportion^2)
	abundance_tot_orig=mean(apply(N_orig[,'ocean',,f,'compet',m],1,sum,na.rm=T))
	abundance_per_species_orig=apply(N_orig[,'ocean',,f,'compet',m],2,mean,na.rm=T)
	proportion_orig=abundance_per_species_orig/abundance_tot_orig
	simpson_orig[f,m]=sum(proportion_orig^2)
}
}
plot(fac_simu,1/simpson[,1],ylim=c(0,5),ylab="1/Simpson",log="x",xlab=expression("a"["ij,sim"]/"a"["ij,ref"]),col="black",t="l",lty=1,cex.axis=1.5,cex.lab=1.5,lwd=2)

print("Competition, modelI, no cyst bank")
print(range(1/simpson[,1]))
lines(fac_simu,1/simpson_orig[,1],col="black",lty=2,lwd=2)
print("Competition, modelI, WITH cyst bank")
print(range(1/simpson_orig[,1]))
lines(fac_simu,1/simpson[,2],col="grey",lty=1,lwd=2)
print("Competition, modelII, no cyst bank")
print(range(1/simpson[,2]))
lines(fac_simu,1/simpson_orig[,2],col="grey",lty=2,lwd=2)
print("Competition, modelII, WITH cyst bank")
print(range(1/simpson_orig[,2]))
mtext("c",2,line=-0.6,at=5.5,cex=1.25,las=2,font=2)

#SIMPSON
par(mar=c(4,2.5,1.25,2))
simpson=matrix(NA,length(fac_simu),2)
simpson_orig=matrix(NA,length(fac_simu),2)
for(f in 1:length(fac_simu)){
for(m in 1:2){
        abundance_tot=mean(apply(N_array[,'ocean',,f,'facil',m],1,sum,na.rm=T))
        abundance_per_species=apply(N_array[,'ocean',,f,'facil',m],2,mean,na.rm=T)
        proportion=abundance_per_species/abundance_tot
        simpson[f,m]=sum(proportion^2)
        abundance_tot_orig=mean(apply(N_orig[,'ocean',,f,'facil',m],1,sum,na.rm=T))
        abundance_per_species_orig=apply(N_orig[,'ocean',,f,'facil',m],2,mean,na.rm=T)
        proportion_orig=abundance_per_species_orig/abundance_tot_orig
        simpson_orig[f,m]=sum(proportion_orig^2)
}
}
plot(fac_simu,1/simpson[,1],ylim=c(0,5),ylab="",log="x",xlab=expression("a"["ij,sim"]/"a"["ij,ref"]),col="black",t="l",lty=1,cex.axis=1.5,cex.lab=1.5,lwd=2,yaxt="n")
text(0.16,5.5,"d",3,las=2,xpd=NA,cex=1.5,font=2)
lines(fac_simu,1/simpson_orig[,1],col="black",lty=2,lwd=2)
lines(fac_simu,1/simpson[,2],col="grey",lty=1,lwd=2)
lines(fac_simu,1/simpson_orig[,2],col="grey",lty=2,lwd=2)

print("Facilitation, modelI, no cyst bank")
print(range(1/simpson[,1]))
print("Facilitation, modelI, WITH cyst bank")
print(range(1/simpson_orig[,1]))
print("Facilitation, modelII, no cyst bank")
print(range(1/simpson[,2],na.rm=T))
print("Facilitation, modelII, WITH cyst bank")
print(range(1/simpson_orig[,2]))



legend("topleft",c("W bank","W/o bank","Model I","Model II"),lty=c(2,1,1,1),col=c("black","black","black","grey"),bty="n",cex=1.25,lwd=2)
dev.off()


surviv_compet=apply(apply(N_array[id_persistence,"ocean",,,'compet',]==0,c(2,3,4),sum)<nb_persistence,c(1,3),sum)/(length(fac_simu))
surviv_facil=apply(apply(N_array[id_persistence,"ocean",,,'facil',]==0,c(2,3,4),sum)<nb_persistence,c(1,3),sum)/length(fac_simu)

pdf("dynamics_vs_survival_no_seed_bank_minab_amplitude_niche.pdf",width=9,height=3)
set.seed(42)
par(mfrow=c(1,3),oma=c(1,1,1,1),mar=c(4.5,4.5,1,1))
a=jitter(min_val[,1],amount=0.1)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Log10(min abundance)",pch=16,ylim=c(0,1),ylab="Prob survival",cex=1.5,xlim=range(c(min_val))+c(-0.1,0.1),cex.lab=1.5,cex.axis=1.5)
text(-1.415,1.08,"a",las=2,xpd=NA,cex=1.5,font=2)
b=jitter(min_val[,2],amount=0.1)
points(b,surviv_compet[,2],col="grey",pch=17,cex=1.75)

a=jitter(range_val[,1],amount=0.5)
plot(a,surviv_compet[,1],t="p",col="black",pch=16,xlab="Log. amplitude",ylab="",ylim=c(0,1),cex=1.5,cex.axis=1.5,cex.lab=1.5)
text(0.85,1.09,"b",las=2,xpd=NA,cex=1.5,font=2)
b=jitter(range_val[,2],amount=0.5)
points(b,surviv_compet[,2],col="grey",pch=17,cex=1.5)

a=jitter(log10(B),amount=0.5)
plot(a,surviv_compet[,1],t="p",col="black",pch=16,xlab="Proxy niche width",ylab="",ylim=c(0,1),cex=1.5,xlim=range(c(log10(B)))+c(-0.1,0.5),cex.lab=1.5,cex.axis=1.5)
text(0.7,1.08,"c",las=2,xpd=NA,cex=1.5,font=2)
b=jitter(log10(B),amount=0.5)
points(b,surviv_compet[,2],col="grey",pch=17,cex=1.5)

legend("bottomright",c("Model I","Model II"),pch=c(16,17),col=c("black","grey"),pt.lwd=c(1,1,1.5,1.5),bty="n",cex=1.25)
dev.off()

pdf("dynamics_vs_survival_no_seed_bank_meanab_sinking_temperature.pdf",width=9,height=3)
set.seed(42)
par(mfrow=c(1,3),oma=c(1,1,1,1),mar=c(4.5,4.5,1,1))

a=jitter(mean_val[,1],amount=0.1)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Log10(mean abundance)",pch=16,ylim=c(0,1),ylab="Prob survival",cex=1.5,xlim=range(c(mean_val))+c(-0.1,0.1),cex.lab=1.5,cex.axis=1.5)
text(2.2,1.08,"a",las=2,xpd=NA,cex=1.5,font=2)
b=jitter(mean_val[,2],amount=0.1)
points(b,surviv_compet[,2],col="grey",pch=17,cex=1.75)

a=jitter(S,amount=0.01)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Sinking rate",pch=16,ylim=c(0,1),ylab="Prob survival",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlim=range(c(S))+c(-0.05,0.05))
text(-0.5,1.08,"b",las=2,xpd=NA,cex=1.5,font=2)
points(a,surviv_compet[,2],col="grey",pch=17,cex=1.75)

a=jitter(T_opt,amount=0.5)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Optimal temperature",pch=16,ylim=c(0,1),ylab="Prob survival",cex=1.5,cex.lab=1.5,cex.axis=1.5,xlim=range(c(T_opt))+c(-0.5,0.5))
text(286.75,1.08,"c",las=2,xpd=NA,cex=1.5,font=2)
points(a,surviv_compet[,2],col="grey",pch=17,cex=1.75)

legend("topright",c("Model I","Model II"),pch=c(16,17),col=c("black","grey"),pt.lwd=c(1,1,1.5,1.5),bty="n")
dev.off()


##Assessing diversity in the cyst bank
#pdf("cyst_bank_effect_compet.pdf")
par(mfrow=c(2,2),mar=c(1.75,4.5,3,0))
 #Final richness, taking into account the last months, instead of only the end of the simulation
plot(fac_simu,apply(apply(N_orig[id_persistence,"seed",,,'compet','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),ylim=c(0,11.5),ylab="Richness",xlab="",main="Competition",log="x",col="black",lty=2,t="l",xaxt="n",cex.axis=1.5,cex.lab=1.5,lwd=2,cex.main=1.75)
mtext("a",2,line=-0.6,at=12.6,cex=1.25,las=2,font=2)
lines(fac_simu,apply(apply(N_orig[id_persistence,"seed",,,'compet','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="grey",lty=2,lwd=2)

lines(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'compet','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),col="darkblue",lty=1,t="l",xaxt="n",cex.axis=1.5,cex.lab=1.5,lwd=2,cex.main=1.75)
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'compet','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),col="darkblue",lty=2,lwd=2)
lines(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'compet','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="lightblue",lty=1,lwd=2)
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'compet','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="lightblue",lty=2,lwd=2)


par(mar=c(1.75,2.5,3,2))
plot(fac_simu,apply(apply(N_orig[id_persistence,"seed",,,'facil','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),ylim=c(0,11.5),ylab="",xlab="",main="Facilitation",log="x",col="black",lty=2,t="l",xaxt="n",yaxt="n",lwd=2,cex.axis=1.5,cex.lab=1.5,cex.main=1.75)
mtext("b",2,line=-0.6,at=12.6,cex=1.25,las=2,font=2)
lines(fac_simu,apply(apply(N_orig[id_persistence,"seed",,,'facil','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="grey",lty=2,lwd=2)

lines(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'facil','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),col="darkblue",lty=1,t="l",xaxt="n",cex.axis=1.5,cex.lab=1.5,lwd=2,cex.main=1.75)
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'facil','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),col="darkblue",lty=2,lwd=2)
lines(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'facil','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="lightblue",lty=1,lwd=2)
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'facil','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="lightblue",lty=2,lwd=2)

par(mar=c(4,4.5,1.25,0))

N_tot=matrix(NA,length(fac_simu),2)
N_tot_orig=matrix(NA,length(fac_simu),2)

#SIMPSON
simpson=matrix(NA,length(fac_simu),2)
simpson_orig=matrix(NA,length(fac_simu),2)
for(f in 1:length(fac_simu)){
for(m in 1:2){
        abundance_tot_orig=mean(apply(N_orig[,'seed',,f,'compet',m],1,sum,na.rm=T))
        abundance_per_species_orig=apply(N_orig[,'seed',,f,'compet',m],2,mean,na.rm=T)
        proportion_orig=abundance_per_species_orig/abundance_tot_orig
        simpson_orig[f,m]=sum(proportion_orig^2)
}
}
plot(fac_simu,1/simpson_orig[,1],ylim=c(0,5),ylab="1/Simpson",log="x",xlab=expression("a"["ij,sim"]/"a"["ij,ref"]),col="black",t="l",lty=2,cex.axis=1.5,cex.lab=1.5,lwd=2)

lines(fac_simu,1/simpson_orig[,1],col="black",lty=2,lwd=2)
lines(fac_simu,1/simpson_orig[,2],col="grey",lty=2,lwd=2)
mtext("c",2,line=-0.6,at=5.5,cex=1.25,las=2,font=2)

for(f in 1:length(fac_simu)){
for(m in 1:2){
        abundance_tot=mean(apply(N_array[,'ocean',,f,'compet',m],1,sum,na.rm=T))
        abundance_per_species=apply(N_array[,'ocean',,f,'compet',m],2,mean,na.rm=T)
        proportion=abundance_per_species/abundance_tot
        simpson[f,m]=sum(proportion^2)
        abundance_tot_orig=mean(apply(N_orig[,'ocean',,f,'compet',m],1,sum,na.rm=T))
        abundance_per_species_orig=apply(N_orig[,'ocean',,f,'compet',m],2,mean,na.rm=T)
        proportion_orig=abundance_per_species_orig/abundance_tot_orig
        simpson_orig[f,m]=sum(proportion_orig^2)
}
}
lines(fac_simu,1/simpson[,1],col="darkblue",t="l",lty=1,cex.axis=1.5,cex.lab=1.5,lwd=2)
lines(fac_simu,1/simpson_orig[,1],col="darkblue",lty=2,lwd=2)
lines(fac_simu,1/simpson[,2],col="lightblue",lty=1,lwd=2)
lines(fac_simu,1/simpson_orig[,2],col="lightblue",lty=2,lwd=2)

for(f in 1:length(fac_simu)){
for(m in 1:2){
        abundance_tot=mean(apply(N_array[,'coast',,f,'compet',m],1,sum,na.rm=T))
        abundance_per_species=apply(N_array[,'coast',,f,'compet',m],2,mean,na.rm=T)
        proportion=abundance_per_species/abundance_tot
        simpson[f,m]=sum(proportion^2)
        abundance_tot_orig=mean(apply(N_orig[,'coast',,f,'compet',m],1,sum,na.rm=T))
        abundance_per_species_orig=apply(N_orig[,'coast',,f,'compet',m],2,mean,na.rm=T)
        proportion_orig=abundance_per_species_orig/abundance_tot_orig
        simpson_orig[f,m]=sum(proportion_orig^2)
}
}
legend("topleft",c("Cyst","Ocean"),lty=c(1,1),col=c("black","darkblue"),bty="n",cex=1.25,lwd=2)



#SIMPSON
par(mar=c(4,2.5,1.25,2))
simpson=matrix(NA,length(fac_simu),2)
simpson_orig=matrix(NA,length(fac_simu),2)
for(f in 1:length(fac_simu)){
for(m in 1:2){
        abundance_tot_orig=mean(apply(N_orig[,'seed',,f,'facil',m],1,sum,na.rm=T))
        abundance_per_species_orig=apply(N_orig[,'seed',,f,'facil',m],2,mean,na.rm=T)
        proportion_orig=abundance_per_species_orig/abundance_tot_orig
        simpson_orig[f,m]=sum(proportion_orig^2)
}
}
plot(fac_simu,1/simpson_orig[,1],ylim=c(0,5),ylab="",log="x",xlab=expression("a"["ij,sim"]/"a"["ij,ref"]),col="black",t="l",lty=2,cex.axis=1.5,cex.lab=1.5,lwd=2,yaxt="n")
text(0.16,5.5,"d",3,las=2,xpd=NA,cex=1.5,font=2)
lines(fac_simu,1/simpson_orig[,2],col="grey",lty=2,lwd=2)




print("Facilitation, modelI, no cyst bank")
print(range(1/simpson[,1]))
print("Facilitation, modelI, WITH cyst bank")
print(range(1/simpson_orig[,1]))
print("Facilitation, modelII, no cyst bank")
print(range(1/simpson[,2],na.rm=T))
print("Facilitation, modelII, WITH cyst bank")
print(range(1/simpson_orig[,2]))


for(f in 1:length(fac_simu)){
for(m in 1:2){
        abundance_tot=mean(apply(N_array[,'ocean',,f,'facil',m],1,sum,na.rm=T))
        abundance_per_species=apply(N_array[,'ocean',,f,'facil',m],2,mean,na.rm=T)
        proportion=abundance_per_species/abundance_tot
        simpson[f,m]=sum(proportion^2)
        abundance_tot_orig=mean(apply(N_orig[,'ocean',,f,'facil',m],1,sum,na.rm=T))
        abundance_per_species_orig=apply(N_orig[,'ocean',,f,'facil',m],2,mean,na.rm=T)
        proportion_orig=abundance_per_species_orig/abundance_tot_orig
        simpson_orig[f,m]=sum(proportion_orig^2)
}
}
lines(fac_simu,1/simpson[,1],col="darkblue",t="l",lty=1,cex.axis=1.5,cex.lab=1.5,lwd=2)
lines(fac_simu,1/simpson_orig[,1],col="darkblue",lty=2,lwd=2)
lines(fac_simu,1/simpson[,2],col="lightblue",lty=1,lwd=2)
lines(fac_simu,1/simpson_orig[,2],col="lightblue",lty=2,lwd=2)

for(f in 1:length(fac_simu)){
for(m in 1:2){
        abundance_tot=mean(apply(N_array[,'coast',,f,'facil',m],1,sum,na.rm=T))
        abundance_per_species=apply(N_array[,'coast',,f,'facil',m],2,mean,na.rm=T)
        proportion=abundance_per_species/abundance_tot
        simpson[f,m]=sum(proportion^2)
        abundance_tot_orig=mean(apply(N_orig[,'coast',,f,'facil',m],1,sum,na.rm=T))
        abundance_per_species_orig=apply(N_orig[,'coast',,f,'facil',m],2,mean,na.rm=T)
        proportion_orig=abundance_per_species_orig/abundance_tot_orig
        simpson_orig[f,m]=sum(proportion_orig^2)
}
}


legend("topleft",c("W bank","W/o bank","Model I","Model II"),lty=c(2,1,1,1),col=c("black","black","black","grey"),bty="n",cex=1.25,lwd=2)
#dev.off()


#pdf("species_no_seed_bank.pdf",width=10)
par(mfrow=c(3,4))
for(m in 1:2){
for(s in 1:nspp){
	ydelim=c(-3,max(log10(c(N_array[,"ocean",s,,,m],N_orig[,"ocean",s,,,m])),na.rm=T))
        plot(fac_simu,log10(apply(N_array[,"ocean",s,,'compet',m],2,mean)+10^(-5)),pch=1,ylim=ydelim,main=name_spp[s],col="red",log="x",ylab="Log10(abundance)",xlab="")
        points(fac_simu,log10(apply(N_array[,"ocean",s,,'facil',m],2,mean)+10^(-5)),pch=1,col="blue")
        lines(fac_simu,log10(apply(N_orig[,"ocean",s,,'compet',m],2,mean)+10^(-5)),col="red")
        lines(fac_simu,log10(apply(N_orig[,"ocean",s,,'facil',m],2,mean)+10^(-5)),col="blue")
        abline(h=mean_val[s])
}
plot(0,0,t="n")
legend("top",c("Compet w/o bank","Facil w/o bank","Compet w bank","Facil w bank"),col=c("red","blue"),pch=c(1,1,NA,NA),lty=c(NA,NA,1,1),bty="n")
}
#dev.off()


##### Bray Curtis for cyst/ ocean or coast
matrix_bc=array(NA,dim=c(2,2,length(fac_simu),2,2),dimnames=list(c("origin","disturbed"),c("ocean","coast"),fac_simu,c("compet","facil"),1:2))
for (cpt in c("ocean","coast")){
	for (f in 1:length(fac_simu)){
		for (type in c("compet","facil")){
			for(m in 1:2){
				matrix_bc["origin",cpt,f,type,m]=BC_index(N_orig[,"seed",,f,type,m],N_orig[,cpt,,f,type,m])
				matrix_bc["disturbed",cpt,f,type,m]=BC_index(N_array[,"seed",,f,type,m],N_array[,cpt,,f,type,m])
			}
		}
	}
}

#pdf("bray_curtis_index_compet.pdf")
par(mfrow=c(2,2))
plot(fac_simu,rep(1,length(fac_simu)),t="n",ylim=c(0.6,1),xlab="",ylab="Bray-Curtis index",main="Compet ocean",log="x")
type="compet"
cpt="ocean"
lines(fac_simu,matrix_bc["origin",cpt,,type,1],lty=2,col="black")
lines(fac_simu,matrix_bc["origin",cpt,,type,2],lty=2,col="grey")
lines(fac_simu,matrix_bc["disturbed",cpt,,type,1],lty=1,col="black")
lines(fac_simu,matrix_bc["disturbed",cpt,,type,2],lty=1,col="grey")
legend("topright",c("W bank","W/o bank","Model I","Model II"),lty=c(2,1,1,1),col=c("black","black","black","grey"),bty="n",cex=1.,lwd=2)


plot(fac_simu,rep(1,length(fac_simu)),t="n",ylim=c(0.6,1),xlab="",ylab="",main="Facil ocean",log="x")
type="facil"
cpt="ocean"
lines(fac_simu,matrix_bc["origin",cpt,,type,1],lty=2,col="black")
lines(fac_simu,matrix_bc["origin",cpt,,type,2],lty=2,col="grey")
lines(fac_simu,matrix_bc["disturbed",cpt,,type,1],lty=1,col="black")
lines(fac_simu,matrix_bc["disturbed",cpt,,type,2],lty=1,col="grey")

plot(fac_simu,rep(1,length(fac_simu)),t="n",ylim=c(0.6,1),xlab=expression("a"["ij,sim"]/"a"["ij,ref"]),ylab="Bray-Curtis index",main="Compet coast",log="x")
type="compet"
cpt="coast"
lines(fac_simu,matrix_bc["origin",cpt,,type,1],lty=2,col="black")
lines(fac_simu,matrix_bc["origin",cpt,,type,2],lty=2,col="grey")
lines(fac_simu,matrix_bc["disturbed",cpt,,type,1],lty=1,col="black")
lines(fac_simu,matrix_bc["disturbed",cpt,,type,2],lty=1,col="grey")


plot(fac_simu,rep(1,length(fac_simu)),t="n",ylim=c(0.6,1),xlab=expression("a"["ij,sim"]/"a"["ij,ref"]),main="Facil coast",log="x")
type="facil"
cpt="coast"
lines(fac_simu,matrix_bc["origin",cpt,,type,1],lty=2,col="black")
lines(fac_simu,matrix_bc["origin",cpt,,type,2],lty=2,col="grey")
lines(fac_simu,matrix_bc["disturbed",cpt,,type,1],lty=1,col="black")
lines(fac_simu,matrix_bc["disturbed",cpt,,type,2],lty=1,col="grey")
#dev.off()
