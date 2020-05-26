####
# 18/05/2020 CP: Effect of competition on communities with and without seed bank
# 21/05/2020 CP: General script for both models
###

rm(list=ls())
graphics.off()
source("step_functions.r")

n_iter=1000
nb_year=2
cpt="ocean"
#Data to use (Auger)
a=as.matrix(read.table("../../param/reconstructed_temperature_Auger.txt", row.names=1,header=T,sep=";",dec="."))
temp_model=a[1:n_iter]
mean_model=mean(temp_model)
var_model=var(temp_model)
#Used for summary statistics
tab_pheno=read.table("../../param/generalist_specialist_spp_added_amplitude_season.csv",sep=";",dec=".",header=T,row.names=1)
name_spp=rownames(tab_pheno)
nspp=length(name_spp)
####Two options are possible for mean values
pop_table=read.table("../../param/mean_interpolated_abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
pop_table=pop_table[name_spp,]

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
n_iter=as.numeric(as.character(tab[tab[,1]=="n_iter",2]))
temp_germin=as.numeric(as.character(tab[tab[,1]=="germin_threshold",2]))
a_d=as.numeric(as.character(tab[tab[,1]=="daylength",2]))
threshold=as.numeric(as.character(tab[tab[,1]=="threshold",2]))

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


####Model I
#Interaction
tmp_inter_I=read.table("../../output/modelv1.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
tmp_inter=as.matrix(tmp_inter_I)
list_inter=list(tmp_inter,k_coast2ocean*tmp_inter) #Interaction list contains coastal and oceanic interactions

####ModelII
list_H_tmp=read.table("../../output/modelv2.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
list_H=list(list_H_tmp,list_H_tmp)
type_inter_tmp=read.table("../../output/modelv2.1/facilitation_and_competition.csv",sep=";",dec=".")
type_inter=list(type_inter_tmp,k_coast2ocean*type_inter_tmp)

######################Param changed when using seed_bank
M_orig=M
M=1
nb_simu=100
mean_tmp=c(0,2,5,7) #Temperature increases
vari_tmp=c(1,0.75,1.25) #Not sure variance increases in all climate models
id=(n_iter-365):n_iter

###Initialize
N_array=array(NA,dim=c(n_iter,3,nspp,nb_simu,length(mean_tmp),length(vari_tmp),2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,1:nb_simu,mean_tmp,vari_tmp,c("Model I","Model II")))
N_array[1,,,,,,]=10^3
N_orig=array(NA,dim=c(n_iter,3,nspp,nb_simu,length(mean_tmp),length(vari_tmp),2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,1:nb_simu,mean_tmp,vari_tmp,c("Model I","Model II")))
N_orig[1,,,,,,]=10^3

theta=1.3
for(m in 1:length(mean_tmp)){
	print(paste("mean temp",mean_tmp[m]))
	mean_simu=mean_model+mean_tmp[m]
	for(v in 1:length(vari_tmp)){
		sd_simu=sqrt(var_model*vari_tmp[v])
		for(s in 1:nb_simu){
			temp_model=mean_simu+theta*sd_simu*sin(2*pi*1:n_iter/365.25)+rnorm(n_iter,0,sd_simu*sqrt(1-theta^2/2))
			for(t in 1:(n_iter-1)){
				#ModelI
				#Without Seed Bank
                		var_tmp=step1_modelI(N_array[t,,,s,m,v,"Model I"],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
  		              	Ntmp=var_tmp[[1]]
                		N_array[t+1,,,s,m,v,"Model I"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
		
				#With Seed Bank
		                var_tmp=step1_modelI(N_orig[t,,,s,m,v,"Model I"],list_inter,temp_model[t],M_orig,morta,a_d,T_opt,B,threshold)
                		Ntmp=var_tmp[[1]]
               			N_orig[t+1,,,s,m,v,"Model I"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)

				#Model II
				#Without Seed Bank
                                var_tmp=step1_modelII(N_array[t,,,s,m,v,"Model II"],list_H,type_inter,temp_model[t],M,morta,a_d,T_opt,B)
                                Ntmp=var_tmp[[1]]
                                N_array[t+1,,,s,m,v,"Model II"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)

                                #With Seed Bank
                                var_tmp=step1_modelII(N_orig[t,,,s,m,v,"Model II"],list_H,type_inter,temp_model[t],M_orig,morta,a_d,T_opt,B)
                                Ntmp=var_tmp[[1]]
                                N_orig[t+1,,,s,m,v,"Model II"]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)

			}
		}
	}
}
save(list = ls(all.names = TRUE), file = "no_seed_temp_simu.RData", envir = .GlobalEnv)

pdf("no_seed_bank_temp.pdf")
par(mfrow=c(2,2),mar=c(4,4,2,1))
plot(1:length(mean_tmp),rep(NA,length(mean_tmp)),ylim=c(0,11.5),main="Increase mean(T)",xaxt="n",ylab="Richness",xlab="",xlim=c(1,length(mean_tmp)+0.3))
axis(1,at=1:length(mean_tmp),labels=NA)
for(m in 1:2){
	if(m==1){
		acol="black"
		shift=0
	}else{
		acol="grey"
		shift=0.1
	}
for(i in 1:length(mean_tmp)){
	tmp_value=rep(NA,length(nb_simu))
	tmp_value_orig=rep(NA,length(nb_simu))
	for(j in 1:nb_simu){
		tmp_value[j]=sum(N_array[n_iter,'ocean',,j,i,1,m]>0)
		tmp_value_orig[j]=sum(N_orig[n_iter,'ocean',,j,i,1,m]>0)
	}
	points(i+shift,mean(tmp_value),pch=16,col=acol)
	arrows(i+shift,mean(tmp_value)-sd(tmp_value)*1.96,i+shift,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0,col=acol)
	points(i+shift,mean(tmp_value_orig),pch=17,col=acol)
	arrows(i+shift,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i+shift,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0,lty=2,col=acol)
}
}

plot(1:length(vari_tmp),rep(NA,length(vari_tmp)),ylim=c(0,11.5),main="Increase var(T)",xaxt="n",xlab="",ylab="",xlim=c(1,length(vari_tmp)+0.3))
axis(1,at=1:length(vari_tmp),labels=NA)
for(m in 1:2){
	if(m==1){
		acol="black"
		shift=0	
	}else{
		acol="grey"
		shift=0.1
	}
for(i in 1:length(vari_tmp)){
        tmp_value=rep(NA,length(nb_simu))
        tmp_valuei_orig=rep(NA,length(nb_simu))
        for(j in 1:nb_simu){
                tmp_value[j]=sum(N_array[n_iter,'ocean',,j,1,i,m]>0)
                tmp_value_orig[j]=sum(N_orig[n_iter,'ocean',,j,1,i,m]>0)
        }
        points(i+shift,mean(tmp_value),pch=16,col=acol)
        arrows(i+shift,mean(tmp_value)-sd(tmp_value)*1.96,i+shift,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0,col=acol)
	points(i+shift,mean(tmp_value_orig),pch=17,col=acol)
	arrows(i+shift,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i+shift,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0,col=acol,lty=2)
}
}
legend("bottom",c("W/o seed bank","W seed bank","Model I","Model II"),col=c("black","black","black","grey"),pch=c(16,17,16,16),lty=c(1,2,1,1),bty="n")

plot(1:length(mean_tmp),rep(NA,length(mean_tmp)),ylim=c(4,6),ylab="Log10(abundance)",xlab="Added tremp",xaxt="n",xlim=c(1,length(mean_tmp)+0.3))
axis(1,at=1:length(mean_tmp),labels=mean_tmp)
for(m in 1:2){
	if(m==1){
		acol="black"
		shift=0
	}else{
		acol="grey"
		shift=0.2
	}
for(i in 1:length(mean_tmp)){
        tmp_value=rep(NA,length(nb_simu))
        tmp_value_orig=rep(NA,length(nb_simu))
        for(j in 1:nb_simu){
                tmp_val=apply(N_array[id,'ocean',,j,i,1,m],1,sum)
		tmp_value[j]=log10(mean(tmp_val))
                tmp_val=apply(N_orig[id,'ocean',,j,i,1,m],1,sum)
		tmp_value_orig[j]=log10(mean(tmp_val))
        }
        points(i+shift,mean(tmp_value),pch=16,col=acol)
        arrows(i+shift,mean(tmp_value)-sd(tmp_value)*1.96,i+shift,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0,col=acol)
	points(i+shift+0.05,mean(tmp_value_orig),pch=17,col=acol)
	arrows(i+shift+0.05,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i+shift+0.05,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0,col=acol,lty=2)
}
}

plot(1:length(vari_tmp),rep(NA,length(vari_tmp)),ylim=c(4,6),ylab="",xlab="Multiplying variance",xaxt="n",xlim=c(1,length(vari_tmp)+0.3))
axis(1,at=1:length(vari_tmp),labels=vari_tmp)
for(m in 1:2){
	if(m==1){
		acol="black"
		shift=0
	}else{
		acol="grey"
		shift=0.2
	}
for(i in 1:length(vari_tmp)){
        tmp_value=rep(NA,length(nb_simu))
        tmp_value_orig=rep(NA,length(nb_simu))
        for(j in 1:nb_simu){
                tmp_val=apply(N_array[id,'ocean',,j,1,i,m],1,sum)
                tmp_value[j]=log10(mean(tmp_val))
                tmp_val=apply(N_orig[id,'ocean',,j,1,i,m],1,sum)
		tmp_value_orig[j]=log10(mean(tmp_val))
        }
        points(i+shift,mean(tmp_value),pch=16,col=acol)
        arrows(i+shift,mean(tmp_value)-sd(tmp_value)*1.96,i+shift,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0,col=acol)
	points(i+shift+0.05,mean(tmp_value_orig),pch=17,col=acol)
	arrows(i+shift+0.05,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i+shift+0.05,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0,col=acol,lty=2)
}
}
dev.off()