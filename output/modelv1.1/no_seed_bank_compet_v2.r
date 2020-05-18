####
# 18/05/2020 CP: Effect of competition on communities with and without seed bank
###

rm(list=ls())
graphics.off()
set.seed(42)

source("step_functions.r")


#######################Param used in every simulations
#Fixed parameters
tab=read.table("simu.csv",sep=";",dec=".",header=T)
cyst_mortality=as.numeric(as.character(tab[tab[,1]=="cyst_mortality",2]))
cyst_burial=as.numeric(as.character(tab[tab[,1]=="cyst_burial",2]))
M=cyst_mortality+cyst_burial
#k_coast2ocean=as.numeric(as.character(tab[tab[,1]=="k_coast2ocean",2])) #Not used anymore
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
list_inter=list(tmp_inter,tmp_inter) #Interaction list contains coastal and oceanic interactions

#Data to use (Auger)
a=as.matrix(read.table("../../param/reconstructed_temperature_Auger.txt", row.names=1,header=T,sep=";",dec="."))
temp_model=a[1:n_iter]

#Sinking rates and T_opt
tab=read.table(paste("../../param/species_specific_parameters.txt",sep=""),sep=";",dec=".",header=T)
tab=tab[name_spp,]
S=S_max*tab$S
T_opt=tab$T_opt+273+5
names(S)=name_spp
names(T_opt)=name_spp
B=tab$Val_b
names(B)=name_spp

#Original dataset
N_original_test_coast=read.table("out_coast.csv",sep=";",row.names=1,header=T)
N_original_test_ocean=read.table("out_ocean.csv",sep=";",row.names=1,header=T)
n_iter=nrow(N_original_test_coast)
id=(n_iter-365):n_iter
mean_val=log10(apply(N_original_test_coast[id,],2,mean))
min_val=log10(apply(N_original_test_coast[id,],2,min))
range_val=log10(apply(apply(N_original_test_coast[id,],2,range),2,diff))


######################Param changed when using seed_bank
M_orig=M
M=1
nb_div=100
simu=seq(1,10,length.out=nb_div)
fac_simu=unique(c(rev(1/simu),simu))

###Initialize
N_array=array(NA,dim=c(n_iter,3,nspp,length(fac_simu),2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(fac_simu,digits=3),c("compet","facil")))
N_array[1,,,,]=10^3
N_orig=array(NA,dim=c(n_iter,3,nspp,length(fac_simu),2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(fac_simu,digits=3),c("compet","facil")))
N_orig[1,,,,]=10^3

for(fac in 1:length(fac_simu)){
	print(fac_simu[fac])
	#Compet
	tmp_inter_compet=tmp_inter
	tmp_inter_compet[which(tmp_inter>0,arr.ind=T)]=tmp_inter[which(tmp_inter>0,arr.ind=T)]*fac_simu[fac]
	list_inter_compet=list(tmp_inter_compet,tmp_inter_compet)

	#Facilitation
	tmp_inter_facil=tmp_inter
	tmp_inter_facil[which(tmp_inter<0,arr.ind=T)]=tmp_inter[which(tmp_inter<0,arr.ind=T)]*fac_simu[fac]
	list_inter_facil=list(tmp_inter_facil,tmp_inter_facil)

	for(t in 1:(n_iter-1)){

		#Without Seed Bank
		#Compet
                var_tmp=step1(N_array[t,,,fac,"compet"],list_inter_compet,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_array[t+1,,,fac,"compet"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
                
		#Facil
		var_tmp=step1(N_array[t,,,fac,"facil"],list_inter_facil,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_array[t+1,,,fac,"facil"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
		
		#With Seed Bank
		#Compet
                var_tmp=step1(N_orig[t,,,fac,"compet"],list_inter_compet,temp_model[t],M_orig,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_orig[t+1,,,fac,"compet"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
                
		#Facil
		var_tmp=step1(N_orig[t,,,fac,"facil"],list_inter_facil,temp_model[t],M_orig,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_orig[t+1,,,fac,"facil"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)


	}
		if(fac_simu[fac]==1){
			print(N_original_test_ocean[n_iter,])
			print(N_orig[n_iter,"ocean",,fac,"compet"])
			print(N_orig[n_iter,"ocean",,fac,"facil"])

		}
}

pdf("no_seed_bank_compet.pdf")
par(mfrow=c(2,2))
plot(fac_simu,apply(N_array[n_iter,'ocean',,,'compet']>0,2,sum),ylim=c(0,11.5),ylab="Richness",xlab="",main="Competition",log="x")
lines(fac_simu,apply(N_orig[n_iter,'ocean',,,'compet']>0,2,sum))
#abline(v=1)
#axis(1,at=floor(seq(1,length(fac_simu),length.out=10)),format(fac_simu[floor(seq(1,length(fac_simu),length.out=10))],digits=2))

plot(fac_simu,apply(N_array[n_iter,'ocean',,,'facil']>0,2,sum),ylim=c(0,11.5),main="Facilitation",log="x",xlab="",ylab="")
lines(fac_simu,apply(N_orig[n_iter,'ocean',,,'facil']>0,2,sum))
#abline(v=which(fac_simu==1))
#axis(1,at=floor(seq(1,length(fac_simu),length.out=10)),format(fac_simu[floor(seq(1,length(fac_simu),length.out=10))],digits=2))

N_tot=rep(NA,length(fac_simu))
N_tot_orig=rep(NA,length(fac_simu))
for(f in 1:length(fac_simu)){
	tmp=apply(N_array[id,'ocean',,f,'compet'],1,sum)
	N_tot[f]=log10(mean(tmp))	
	tmp=apply(N_orig[id,'ocean',,f,'compet'],1,sum)
	N_tot_orig[f]=log10(mean(tmp))	
}
plot(fac_simu,N_tot,ylim=range(c(N_tot,N_tot_orig)),ylab="Log10(abundance)",xlab="Multiplying factor",log="x")
lines(fac_simu,N_tot_orig)
#abline(v=which(fac_simu==1))
#axis(1,at=floor(seq(1,length(fac_simu),length.out=10)),format(fac_simu[floor(seq(1,length(fac_simu),length.out=10))],digits=2))

N_tot=rep(NA,length(fac_simu))
N_tot_orig=rep(NA,length(fac_simu))
for(f in 1:length(fac_simu)){
        tmp=apply(N_array[id,'ocean',,f,'facil'],1,sum)
        N_tot[f]=log10(mean(tmp))
        tmp=apply(N_orig[id,'ocean',,f,'facil'],1,sum)
        N_tot_orig[f]=log10(mean(tmp))
}
plot(fac_simu,N_tot,ylim=range(c(N_tot,N_tot_orig)),log="x",ylab="",xlab="Multiplying factor")
lines(fac_simu,N_tot_orig)
legend("topleft",c("W bank","W/o bank"),pch=c(NA,1),lty=c(1,NA),bty="n")
#abline(v=which(fac_simu==1))
#axis(1,at=floor(seq(1,length(fac_simu),length.out=10)),format(fac_simu[floor(seq(1,length(fac_simu),length.out=10))],digits=2))
dev.off()

surviv_compet=apply(N_array[1000,"ocean",,,"compet"]>0,1,sum)/200
surviv_facil=apply(N_array[1000,"ocean",,,"facil"]>0,1,sum)/200
pdf("dynamics_vs_survival_no_seed_bank.pdf",width=12)
par(mfrow=c(1,3))
plot(surviv_compet,mean_val,pch=16,col="red",xlim=c(-0.2,1.2),ylab="Average abundance",xlab="Nb survival")
points(surviv_facil,mean_val,pch=16,col="blue")
text(surviv_compet,mean_val,name_spp,pos=2)
legend("topleft",c("Competition","Facilitation"),col=c("red","blue"),pch=16,bty="n")

plot(surviv_compet,range_val,pch=16,col="red",xlim=c(-0.2,1.2),ylab="Average amplitude",xlab="Nb survival")
points(surviv_facil,range_val,pch=16,col="blue")
text(surviv_compet,range_val,name_spp,pos=2)

plot(surviv_compet,min_val,pch=16,col="red",xlim=c(-0.2,1.2),ylab="Min amplitude",xlab="Nb survival")
points(surviv_facil,min_val,pch=16,col="blue")
text(surviv_compet,min_val,name_spp,pos=2)
dev.off()

pdf("species_no_seed_bank.pdf",width=10)
par(mfrow=c(3,4))

for(s in 1:nspp){
	ydelim=c(-3,max(log10(c(N_array[id,"ocean",s,,],N_orig[id,"ocean",s,,]))))
        plot(fac_simu,log10(apply(N_array[id,"ocean",s,,'compet'],2,mean)+10^(-5)),pch=1,ylim=ydelim,main=name_spp[s],col="red",log="x",ylab="Log10(abundance)",xlab="")
        points(fac_simu,log10(apply(N_array[id,"ocean",s,,'facil'],2,mean)+10^(-5)),pch=1,col="blue")
        lines(fac_simu,log10(apply(N_orig[id,"ocean",s,,'compet'],2,mean)+10^(-5)),col="red")
        lines(fac_simu,log10(apply(N_orig[id,"ocean",s,,'facil'],2,mean)+10^(-5)),col="blue")
        abline(h=mean_val[s])
}
plot(0,0,t="n")
legend("top",c("Compet w/o bank","Facil w/o bank","Compet w bank","Facil w bank"),col=c("red","blue"),pch=c(1,1,NA,NA),lty=c(NA,NA,1,1),bty="n")
dev.off()

