####
# 18/05/2020 CP: Effect of competition on communities with and without seed bank
# 21/05/2020 CP: General script for both models
###

doyouload=T #TRUE if you want to reload previous results to plot graphs; FALSE if you want to relaunch the analyses (and then plot the graphs)

if(!doyouload){
rm(list=ls())
graphics.off()
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
nb_simu=5
mean_tmp=c(0,2,5,7) #Temperature increases
vari_tmp=c(1,0.75,1.25) #Not sure variance increases in all climate models
id=(n_iter-365):n_iter

###Initialize
N_array=array(NA,dim=c(length(id),3,nspp,nb_simu,length(mean_tmp),length(vari_tmp),2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,1:nb_simu,mean_tmp,vari_tmp,c("Model I","Model II")))
N_orig=array(NA,dim=c(length(id),3,nspp,nb_simu,length(mean_tmp),length(vari_tmp),2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,1:nb_simu,mean_tmp,vari_tmp,c("Model I","Model II")))

theta=1.412
for(m in 1:length(mean_tmp)){
	print(paste("mean temp",mean_tmp[m]))
	mean_simu=mean_model+mean_tmp[m]
	for(v in 1:length(vari_tmp)){
		sd_simu=sqrt(var_model*vari_tmp[v])
		for(s in 1:nb_simu){
			temp_model=mean_simu+theta*sd_simu*sin(2*pi*(1:n_iter)/365+pi+pi/3)+rnorm(n_iter,0,sd_simu*sqrt(1-theta^2/2))
			N_array_simu=array(NA,dim=c(n_iter,3,nspp,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,c("Model I","Model II")))
			N_array_simu[1,,,]=10^3
			N_orig_simu=array(NA,dim=c(n_iter,3,nspp,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,c("Model I","Model II")))
			N_orig_simu[1,,,]=10^3
			for(t in 1:(n_iter-1)){
				#ModelI
				#Without Seed Bank
                		var_tmp=step1_modelI(N_array_simu[t,,,"Model I"],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
  		              	Ntmp=var_tmp[[1]]
                		N_array_simu[t+1,,,"Model I"]=step2(Ntmp,S,Gamma,e)
		
				#With Seed Bank
		                var_tmp=step1_modelI(N_orig_simu[t,,,"Model I"],list_inter,temp_model[t],M_orig,morta,a_d,T_opt,B,threshold)
                		Ntmp=var_tmp[[1]]
               			N_orig_simu[t+1,,,"Model I"]=step2(Ntmp,S,Gamma,e)

				#Model II
				#Without Seed Bank
                                var_tmp=step1_modelII(N_array_simu[t,,,"Model II"],list_H,type_inter,temp_model[t],M,morta,a_d,T_opt,B)
                                Ntmp=var_tmp[[1]]
                                N_array_simu[t+1,,,"Model II"]=step2(Ntmp,S,Gamma,e)

                                #With Seed Bank
                                var_tmp=step1_modelII(N_orig_simu[t,,,"Model II"],list_H,type_inter,temp_model[t],M_orig,morta,a_d,T_opt,B)
                                Ntmp=var_tmp[[1]]
                                N_orig_simu[t+1,,,"Model II"]=step2(Ntmp,S,Gamma,e)

			}
                        N_array[,,,s,m,v,]=N_array_simu[id,,,]
                        N_orig[,,,s,m,v,]=N_orig_simu[id,,,]
		}
	}
}
save(list = ls(all.names = TRUE), file = "no_seed_temp_simu.RData", envir = .GlobalEnv)
}else{##end !doyouload
load("no_seed_temp_simu.RData")
}
id_persistence=seq(length(id)-nb_persistence+1,length(id))

pdf("no_seed_bank_temp.pdf")
par(mfrow=c(2,2),mar=c(2,4.5,2,0))
plot(1:length(mean_tmp),rep(NA,length(mean_tmp)),ylim=c(0,11.5),xaxt="n",ylab="Richness",xlab="",xlim=c(0.8,length(mean_tmp)+0.5),cex.axis=1.5,cex.lab=1.5)
mtext("a",line=0.2,at=.75,cex=1.3,font=2)
#axis(1,at=1:length(mean_tmp),labels=NA)
for(m in 1:2){
	if(m==1){
		acol="black"
		shift=0
	}else{
		acol="grey"
		shift=0.3
	}
for(i in 1:length(mean_tmp)){
	tmp_value=rep(NA,length(nb_simu))
	tmp_value_orig=rep(NA,length(nb_simu))
	for(j in 1:nb_simu){
		#tmp_value[j]=apply(apply(N_array[id_persistence,"ocean",,j,i,1,m]==0,c(2,3),sum)<nb_persistence,2,sum)
		#tmp_value_orig[j]=apply(apply(N_orig[id_persistence,"ocean",,j,i,1,m]==0,c(2,3),sum)<nb_persistence,2,sum)
		tmp_value[j]=sum(apply(N_array[id_persistence,"ocean",,j,i,1,m]==0,2,sum)<nb_persistence)
		tmp_value_orig[j]=sum(apply(N_orig[id_persistence,"ocean",,j,i,1,m]==0,c(2),sum)<nb_persistence)
	}
	points(i+shift,mean(tmp_value),pch=16,col=acol,cex=2)
	arrows(i+shift,mean(tmp_value)-sd(tmp_value)*1.96,i+shift,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0,col=acol,lwd=2)
	points(i+shift+0.1,mean(tmp_value_orig),pch=17,col=acol,cex=2)
	arrows(i+shift+0.1,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i+shift+0.1,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0,lty=2,col=acol,lwd=2)
}
}

par(mar=c(2,2.5,2,2))

plot(1:length(vari_tmp),rep(NA,length(vari_tmp)),ylim=c(0,11.5),xaxt="n",xlab="",ylab="",xlim=c(0.8,length(vari_tmp)+0.5),cex.axis=1.5,cex.lab=1.5,yaxt="n")
mtext("b",line=0.2,at=0.75,cex=1.3,font=2)

#axis(1,at=1:length(vari_tmp),labels=NA)
for(m in 1:2){
	if(m==1){
		acol="black"
		shift=0	
	}else{
		acol="grey"
		shift=0.3
	}
for(i in 1:length(vari_tmp)){
        tmp_value=rep(NA,length(nb_simu))
        tmp_valuei_orig=rep(NA,length(nb_simu))
        for(j in 1:nb_simu){
		#tmp_value[j]=apply(apply(N_array[id_persistence,"ocean",,j,1,i,m]==0,c(2,3),sum)<nb_persistence,2,sum)
		#tmp_value_orig[j]=apply(apply(N_orig[id_persistence,"ocean",,j,1,i,m]==0,c(2,3),sum)<nb_persistence,2,sum)
		tmp_value[j]=sum(apply(N_array[id_persistence,"ocean",,j,1,i,m]==0,2,sum)<nb_persistence)
		tmp_value_orig[j]=sum(apply(N_orig[id_persistence,"ocean",,j,1,i,m]==0,2,sum)<nb_persistence)
        }
        points(i+shift,mean(tmp_value),pch=16,col=acol,cex=2)
        arrows(i+shift,mean(tmp_value)-sd(tmp_value)*1.96,i+shift,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0,col=acol,lwd=1.5)
	points(i+shift+0.1,mean(tmp_value_orig),pch=17,col=acol,cex=2)
	arrows(i+shift+0.1,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i+shift+0.1,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0,col=acol,lty=2,lwd=1.5)
}
}
legend("bottomleft",c("W/o seed bank","W seed bank","Model I","Model II"),col=c("black","black","black","grey"),pch=c(16,17,16,16),lty=c(1,2,1,1),bty="n",cex=1.25)

par(mar=c(4,4.5,0,0))
plot(1:length(mean_tmp),rep(NA,length(mean_tmp)),ylim=c(4,6),ylab="Log10(total abundance)",xlab="Increase in mean(T)",xaxt="n",xlim=c(0.8,length(mean_tmp)+0.5),cex.axis=1.5,cex.lab=1.5)
mtext("c",line=0.2,at=.75,cex=1.3,font=2)
axis(1,at=1:length(mean_tmp),labels=mean_tmp,cex.axis=1.5,cex.lab=1.5)
for(m in 1:2){
	if(m==1){
		acol="black"
		shift=0
	}else{
		acol="grey"
		shift=0.3
	}
for(i in 1:length(mean_tmp)){
        tmp_value=rep(NA,length(nb_simu))
        tmp_value_orig=rep(NA,length(nb_simu))
        for(j in 1:nb_simu){
                tmp_val=apply(N_array[,'ocean',,j,i,1,m],1,sum,na.rm=T)
		tmp_value[j]=log10(mean(tmp_val))
                tmp_val=apply(N_orig[,'ocean',,j,i,1,m],1,sum,na.rm=T)
		tmp_value_orig[j]=log10(mean(tmp_val))
        }
        points(i+shift,mean(tmp_value),pch=16,col=acol,cex=2)
        arrows(i+shift,mean(tmp_value)-sd(tmp_value)*1.96,i+shift,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0,col=acol,lwd=1.5)
	points(i+shift+0.1,mean(tmp_value_orig),pch=17,col=acol,cex=2)
	arrows(i+shift+0.1,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i+shift+0.1,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0,col=acol,lty=2,lwd=1.5)
}
}

par(mar=c(4,2.5,0,2))

plot(1:length(vari_tmp),rep(NA,length(vari_tmp)),ylim=c(4,6),ylab="",xlab=expression("Var(T)"["sim"]/"Var(T)"["init"]),xaxt="n",xlim=c(0.8,length(vari_tmp)+0.5),yaxt="n",cex.lab=1.5,cex.axis=1.5)
mtext("d",line=0.2,at=0.75,cex=1.3,font=2)
axis(1,at=1:length(vari_tmp),labels=vari_tmp,cex.axis=1.5,cex.lab=1.5)
for(m in 1:2){
	if(m==1){
		acol="black"
		shift=0
	}else{
		acol="grey"
		shift=0.3
	}
for(i in 1:length(vari_tmp)){
        tmp_value=rep(NA,length(nb_simu))
        tmp_value_orig=rep(NA,length(nb_simu))
        for(j in 1:nb_simu){
                tmp_val=apply(N_array[,'ocean',,j,1,i,m],1,sum,na.rm=T)
                tmp_value[j]=log10(mean(tmp_val))
                tmp_val=apply(N_orig[,'ocean',,j,1,i,m],1,sum,na.rm=T)
		tmp_value_orig[j]=log10(mean(tmp_val))
        }
	if(i==1){
		print("Total abundance with seed bank)")
		print(mean(tmp_value_orig))
		print("Total abundance without seed bank")
		print(mean(tmp_value))
	}
        points(i+shift,mean(tmp_value),pch=16,col=acol,cex=2)
        arrows(i+shift,mean(tmp_value)-sd(tmp_value)*1.96,i+shift,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0,col=acol,lwd=1.5)
	points(i+shift+0.1,mean(tmp_value_orig),pch=17,col=acol,cex=2)
	arrows(i+shift+0.1,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i+shift+0.1,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0,col=acol,lty=2,cex=2,lwd=1.5)
}
}
dev.off()
