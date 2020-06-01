####
# 18/05/2020 CP: Effect of competition on communities with and without seed bank
# 21/05/2020 CP: Wrote the code so that it could be applied to both models at the same time
###

#if(1==0){
rm(list=ls())
graphics.off()
set.seed(42)

source("step_functions.r")
nb_year=2
cpt="ocean"
nb_persistence=6*30

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


###############Golden_set
#Initialize
N_original_set=array(NA,dim=c(n_iter,3,nspp,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,c("modelI","modelII")))
#effect_compet=array(NA,dim=c(n_iter,2,nspp),dimnames=list(NULL,c("coast","ocean"),name_spp))
N_original_set[1,,,]=10^3

##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1_modelI(N_original_set[t,,,"modelI"],list_inter_I,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_original_set[t+1,,,'modelI']=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)

                var_tmp=step1_modelII(N_original_set[t,,,'modelII'],list_H,type_inter_II,temp_model[t],M,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
                N_original_set[t+1,,,'modelII']=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}

mean_val=log10(apply(N_original_set[id,"ocean",,],c(2,3),mean))
min_val=log10(apply(N_original_set[id,"ocean",,],c(2,3),min))
range_val=apply(log10(apply(N_original_set[id,"ocean",,],c(2,3),range)),c(2,3),diff)


######################Param changed when using seed_bank
M_orig=M
M=1
nb_div=20
simu=seq(1,10,length.out=nb_div)
fac_simu=unique(c(rev(1/simu),simu))

###Initialize
N_array=array(NA,dim=c(length(id),3,nspp,length(fac_simu),2,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(fac_simu,digits=3),c("compet","facil"),c("ModelI","ModelII")))
N_orig=array(NA,dim=c(length(id),3,nspp,length(fac_simu),2,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(fac_simu,digits=3),c("compet","facil"),c("ModelI","ModelII")))

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
	N_array_simu[1,,,,]=10^3
	N_orig_simu=array(NA,dim=c(n_iter,3,nspp,2,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,c("compet","facil"),c("ModelI","ModelII")))
	N_orig_simu[1,,,,]=10^3
	for(t in 1:(n_iter-1)){
		#ModelI		
		#Without Seed Bank
		#Compet
                var_tmp=step1_modelI(N_array_simu[t,,,"compet","ModelI"],list_inter_compet,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_array_simu[t+1,,,"compet","ModelI"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
                
		#Facil
		var_tmp=step1_modelI(N_array_simu[t,,,"facil","ModelI"],list_inter_facil,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_array_simu[t+1,,,"facil","ModelI"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
		
		#With Seed Bank
		#Compet
                var_tmp=step1_modelI(N_orig_simu[t,,,"compet","ModelI"],list_inter_compet,temp_model[t],M_orig,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_orig_simu[t+1,,,"compet","ModelI"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
                
		#Facil
		var_tmp=step1_modelI(N_orig_simu[t,,,"facil","ModelI"],list_inter_facil,temp_model[t],M_orig,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_orig_simu[t+1,,,"facil","ModelI"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)

		##ModelII
                #Without Seed Bank
                #Compet
                var_tmp=step1_modelII(N_array_simu[t,,,"compet","ModelII"],list_H,type_inter_compet,temp_model[t],M,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
                N_array_simu[t+1,,,"compet","ModelII"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)

                #Facil
                var_tmp=step1_modelII(N_array_simu[t,,,"facil","ModelII"],list_H,type_inter_facil,temp_model[t],M,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
                N_array_simu[t+1,,,"facil","ModelII"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)

                #With Seed Bank
                #Compet
                var_tmp=step1_modelII(N_orig_simu[t,,,"compet","ModelII"],list_H,type_inter_compet,temp_model[t],M_orig,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
                N_orig_simu[t+1,,,"compet","ModelII"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)

                #Facil
                var_tmp=step1_modelII(N_orig_simu[t,,,"facil","ModelII"],list_H,type_inter_facil,temp_model[t],M_orig,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
                N_orig_simu[t+1,,,"facil","ModelII"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)



	}
	N_array[,,,fac,,]=N_array_simu[id,,,,]
	N_orig[,,,fac,,]=N_orig_simu[id,,,,]
}

save(N_array, N_orig, mean_val,min_val,range_val,fac_simu,nb_persistence,id,file = "no_seed_compet_simu_light.RData", envir = .GlobalEnv)
#}#1==0

pdf("no_seed_bank_compet.pdf")
par(mfrow=c(2,2),mar=c(4,4,1,1))
id_persistence=seq(length(id)-nb_persistence+1,length(id)) 
plot(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'compet','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),ylim=c(0,11.5),ylab="Richness",xlab="",main="Competition",log="x",col="black",lty=1,t="l")
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'compet','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),col="black",lty=2)
lines(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'compet','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="grey",lty=1)
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'compet','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="grey",lty=2)
#abline(v=1)
#axis(1,at=floor(seq(1,length(fac_simu),length.out=10)),format(fac_simu[floor(seq(1,length(fac_simu),length.out=10))],digits=2))

plot(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'facil','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),ylim=c(0,11.5),ylab="",xlab="",main="Facilitation",log="x",col="black",lty=1,t="l")
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'facil','ModelI']==0,c(2,3),sum)<nb_persistence,2,sum),col="black",lty=2)
lines(fac_simu,apply(apply(N_array[id_persistence,"ocean",,,'facil','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="grey",lty=1)
lines(fac_simu,apply(apply(N_orig[id_persistence,"ocean",,,'facil','ModelII']==0,c(2,3),sum)<nb_persistence,2,sum),col="grey",lty=2)

#abline(v=which(fac_simu==1))
#axis(1,at=floor(seq(1,length(fac_simu),length.out=10)),format(fac_simu[floor(seq(1,length(fac_simu),length.out=10))],digits=2))

N_tot=matrix(NA,length(fac_simu),2)
N_tot_orig=matrix(NA,length(fac_simu),2)
for(f in 1:length(fac_simu)){
#Log(sum abundance) - arithmetic mean
#	tmp=apply(N_array[id,'ocean',,f,'compet'],1,sum)
#	N_tot[f]=log10(mean(tmp))	
#	tmp=apply(N_orig[id,'ocean',,f,'compet'],1,sum)
#	N_tot_orig[f]=log10(mean(tmp))	
#Sum(log abundance) -geometric mean
for(m in 1:2){
        tmp=apply(log10(N_array[,'ocean',,f,'compet',m]+10^(-5)),1,sum,na.rm=T)/nspp
        N_tot[f,m]=mean(tmp)
        tmp=apply(log10(N_orig[,'ocean',,f,'compet',m]+10^(-5)),1,sum,na.rm=T)/nspp
	N_tot_orig[f,m]=mean(tmp)
}
}
plot(fac_simu,N_tot[,1],ylim=range(c(N_tot,N_tot_orig)),ylab="Log10(mean(abundance))",xlab="Multiplying factor",log="x",col="black",t="l",lty=1)
lines(fac_simu,N_tot_orig[,1],col="black",lty=2)
lines(fac_simu,N_tot[,2],col="grey",lty=1)
lines(fac_simu,N_tot_orig[,2],col="grey",lty=2)
#abline(v=which(fac_simu==1))
#axis(1,at=floor(seq(1,length(fac_simu),length.out=10)),format(fac_simu[floor(seq(1,length(fac_simu),length.out=10))],digits=2))


N_tot=matrix(NA,length(fac_simu),2)
N_tot_orig=matrix(NA,length(fac_simu),2)
for(f in 1:length(fac_simu)){
#Log(sum abundance)
#        tmp=apply(N_array[id,'ocean',,f,'facil'],1,sum)
#        N_tot[f]=log10(mean(tmp))
#        tmp=apply(N_orig[id,'ocean',,f,'facil'],1,sum)
#        N_tot_orig[f]=log10(mean(tmp))
#Sum(log abundance)
	for(m in 1:2){
        tmp=apply(log10(N_array[,'ocean',,f,'facil',m]+10^(-5)),1,sum,na.rm=T)/nspp
        N_tot[f,m]=mean(tmp)
        tmp=apply(log10(N_orig[,'ocean',,f,'facil',m]+10^(-5)),1,sum,na.rm=T)/nspp
        N_tot_orig[f,m]=mean(tmp)
	}
}
plot(fac_simu,N_tot[,1],ylim=range(c(N_tot,N_tot_orig)),log="x",ylab="",xlab="Multiplying factor",lty=1,t="l",col="black")
lines(fac_simu,N_tot_orig[,1],lty=2,col="black")
lines(fac_simu,N_tot[,2],lty=1,col="grey")
lines(fac_simu,N_tot_orig[,2],lty=2,col="grey")
legend("bottomleft",c("W bank","W/o bank","Model I","Model II"),lty=c(2,1,1,1),col=c("black","black","black","grey"),bty="n")
#abline(v=which(fac_simu==1))
#axis(1,at=floor(seq(1,length(fac_simu),length.out=10)),format(fac_simu[floor(seq(1,length(fac_simu),length.out=10))],digits=2))
dev.off()

surviv_compet=apply(apply(N_array[id_persistence,"ocean",,,'compet',]==0,c(2,3,4),sum)<nb_persistence,c(1,3),sum)/(length(fac_simu))
surviv_facil=apply(apply(N_array[id_persistence,"ocean",,,'facil',]==0,c(2,3,4),sum)<nb_persistence,c(1,3),sum)/length(fac_simu)
pdf("dynamics_vs_survival_no_seed_bank.pdf",width=12)
par(mfrow=c(2,3))
plot(mean_val[,1],surviv_compet[,1],t="p",col="red",xlab="Log10(mean abundance)",ylab="Prob survival",pch=16,ylim=c(0,1))
points(mean_val[,1],surviv_facil[,1],col="blue",pch=16)
points(mean_val[,2],surviv_compet[,2],col="red",pch=17)
points(mean_val[,2],surviv_facil[,2],col="blue",pch=17)

plot(min_val[,1],surviv_compet[,1],t="p",col="red",xlab="Log10(min abundance)",ylab="",pch=16,ylim=c(0,1))
points(min_val[,1],surviv_facil[,1],col="blue",pch=16)
points(min_val[,2],surviv_compet[,2],col="red",pch=17)
points(min_val[,2],surviv_facil[,2],col="blue",pch=17)

plot(range_val[,1],surviv_compet[,1],t="p",col="red",pch=16,xlab="Log. amplitude",ylab="",ylim=c(0,1))
points(range_val[,1],surviv_facil[,1],col="blue",pch=16)
points(range_val[,2],surviv_compet[,2],col="red",pch=17)
points(range_val[,2],surviv_facil[,2],col="blue",pch=17)

plot(S,surviv_compet[,1],t="p",col="red",pch=16,xlab="Sinking rate",ylab="Prob survival",ylim=c(0,1))
points(S,surviv_facil[,1],col="blue",pch=16)
points(S,surviv_compet[,2],col="red",pch=17)
points(S,surviv_facil[,2],col="blue",pch=17)

plot(T_opt,surviv_compet[,1],t="p",col="red",pch=16,xlab="T opt",ylab="",ylim=c(0,1))
points(T_opt,surviv_facil[,1],col="blue",pch=16)
points(T_opt,surviv_compet[,2],col="red",pch=17)
points(T_opt,surviv_facil[,2],col="blue",pch=17)

plot(B,surviv_compet[,1],t="p",col="red",pch=16,xlab="Proxy niche width",ylab="",ylim=c(0,1))
points(B,surviv_facil[,1],col="blue",pch=16)
points(B,surviv_compet[,2],col="red",pch=17)
points(B,surviv_facil[,2],col="blue",pch=17)

legend("bottomright",c("Model I","Model II","Compet","Facilitation"),pch=c(16,17,1,1),col=c("black","black","red","blue"),pt.lwd=c(1,1,1.5,1.5),bty="n")
dev.off()


pdf("species_no_seed_bank.pdf",width=10)
par(mfrow=c(3,4))
for(m in 1:2){
for(s in 1:nspp){
	ydelim=c(-3,max(log10(c(N_array[,"ocean",s,,,m],N_orig[,"ocean",s,,,m]))))
        plot(fac_simu,log10(apply(N_array[,"ocean",s,,'compet',m],2,mean)+10^(-5)),pch=1,ylim=ydelim,main=name_spp[s],col="red",log="x",ylab="Log10(abundance)",xlab="")
        points(fac_simu,log10(apply(N_array[,"ocean",s,,'facil',m],2,mean)+10^(-5)),pch=1,col="blue")
        lines(fac_simu,log10(apply(N_orig[,"ocean",s,,'compet',m],2,mean)+10^(-5)),col="red")
        lines(fac_simu,log10(apply(N_orig[,"ocean",s,,'facil',m],2,mean)+10^(-5)),col="blue")
        abline(h=mean_val[s])
}
plot(0,0,t="n")
legend("top",c("Compet w/o bank","Facil w/o bank","Compet w bank","Facil w bank"),col=c("red","blue"),pch=c(1,1,NA,NA),lty=c(NA,NA,1,1),bty="n")
}
dev.off()

