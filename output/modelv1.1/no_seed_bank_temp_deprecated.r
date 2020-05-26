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
mean_model=mean(temp_model)
var_model=var(temp_model)

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
mean_val=apply(N_original_test_coast[id,],2,mean)

######################Param changed when using seed_bank
M_orig=M
M=1
nb_simu=100
mean_tmp=c(0,2,5,7) #Temperature increases
vari_tmp=c(1,0.75,1.25) #Not sure variance increases in all climate models

###Initialize
N_array=array(NA,dim=c(n_iter,3,nspp,nb_simu,length(mean_tmp),length(vari_tmp)),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,1:nb_simu,mean_tmp,vari_tmp))
N_array[1,,,,,]=10^3
N_orig=array(NA,dim=c(n_iter,3,nspp,nb_simu,length(mean_tmp),length(vari_tmp)),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,1:nb_simu,mean_tmp,vari_tmp))
N_orig[1,,,,,]=10^3

theta=1.3
for(m in 1:length(mean_tmp)){
	print(paste("mean temp",mean_tmp[m]))
	mean_simu=mean_model+mean_tmp[m]
	for(v in 1:length(vari_tmp)){
		sd_simu=sqrt(var_model*vari_tmp[v])
		for(s in 1:nb_simu){
			temp_model=mean_simu+theta*sd_simu*sin(2*pi*1:n_iter/365.25)+rnorm(n_iter,0,sd_simu*sqrt(1-theta^2/2))
			for(t in 1:(n_iter-1)){

				#Without Seed Bank
                		var_tmp=step1(N_array[t,,,s,m,v],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
  		              	Ntmp=var_tmp[[1]]
                		N_array[t+1,,,s,m,v]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
		
				#With Seed Bank
		                var_tmp=step1(N_orig[t,,,s,m,v],list_inter,temp_model[t],M_orig,morta,a_d,T_opt,B,threshold)
                		Ntmp=var_tmp[[1]]
               			N_orig[t+1,,,s,m,v]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
			}
		}
	}
}

pdf("no_seed_bank_temp.pdf")
par(mfrow=c(2,2))
plot(1:length(mean_tmp),rep(NA,length(mean_tmp)),ylim=c(0,11.5),main="mean",xaxt="n",ylab="Richness",xlab="")
axis(1,at=1:length(mean_tmp),labels=NA)
for(i in 1:length(mean_tmp)){
	tmp_value=rep(NA,length(nb_simu))
	tmp_value_orig=rep(NA,length(nb_simu))
	for(j in 1:nb_simu){
		tmp_value[j]=sum(N_array[n_iter,'ocean',,j,i,1]>0)
		tmp_value_orig[j]=sum(N_orig[n_iter,'ocean',,j,i,1]>0)
	}
	points(i,mean(tmp_value),pch=16,col="grey")
	arrows(i,mean(tmp_value)-sd(tmp_value)*1.96,i,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0,col="grey")
	points(i,mean(tmp_value_orig),pch=16,col="black")
	arrows(i,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0)
}
legend("bottomright",c("W/o seed bank","W seed bank"),col=c("grey","black",pch=16))

plot(1:length(vari_tmp),rep(NA,length(vari_tmp)),ylim=c(0,11.5),main="var",xaxt="n",xlab="",ylab="")
axis(1,at=1:length(vari_tmp),labels=NA)
for(i in 1:length(vari_tmp)){
        tmp_value=rep(NA,length(nb_simu))
        tmp_valuei_orig=rep(NA,length(nb_simu))
        for(j in 1:nb_simu){
                tmp_value[j]=sum(N_array[n_iter,'ocean',,j,1,i]>0)
                tmp_value_orig[j]=sum(N_orig[n_iter,'ocean',,j,1,i]>0)
        }
        points(i,mean(tmp_value),pch=16,col="grey")
        arrows(i,mean(tmp_value)-sd(tmp_value)*1.96,i,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0)
	points(i,mean(tmp_value_orig),pch=16,col="black")
	arrows(i,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0)
}

plot(1:length(mean_tmp),rep(NA,length(mean_tmp)),ylim=c(-3,8),ylab="Log10(abundance)",xlab="Increase in mean temp",xaxt="n")
axis(1,at=1:length(mean_tmp),labels=mean_tmp)
for(i in 1:length(mean_tmp)){
        tmp_value=rep(NA,length(nb_simu))
        tmp_value_orig=rep(NA,length(nb_simu))
        for(j in 1:nb_simu){
                tmp_val=apply(N_array[id,'ocean',,j,i,1],1,sum)
		tmp_value[j]=log10(mean(tmp_val))
                tmp_val=apply(N_orig[id,'ocean',,j,i,1],1,sum)
		tmp_value_orig[j]=log10(mean(tmp_val))
        }
        points(i,mean(tmp_value),pch=16,col="grey")
        arrows(i,mean(tmp_value)-sd(tmp_value)*1.96,i,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0)
	points(i,mean(tmp_value_orig),pch=16,col="black")
	arrows(i,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0)
}

plot(1:length(vari_tmp),rep(NA,length(vari_tmp)),ylim=c(-3,8),ylab="",xlab="Multiplying variance",xaxt="n")
axis(1,at=1:length(vari_tmp),labels=vari_tmp)
for(i in 1:length(vari_tmp)){
        tmp_value=rep(NA,length(nb_simu))
        tmp_value_orig=rep(NA,length(nb_simu))
        for(j in 1:nb_simu){
                tmp_val=apply(N_array[id,'ocean',,j,1,i],1,sum)
                tmp_value[j]=log10(mean(tmp_val))
                tmp_val=apply(N_orig[id,'ocean',,j,1,i],1,sum)
		tmp_value_orig[j]=log10(mean(tmp_val))
        }
        points(i,mean(tmp_value))
        arrows(i,mean(tmp_value)-sd(tmp_value)*1.96,i,mean(tmp_value)+sd(tmp_value)*1.96,code=3,angle=0)
	points(i,mean(tmp_value_orig),pch=16,col="black")
	arrows(i,mean(tmp_value_orig)-sd(tmp_value_orig)*1.96,i,mean(tmp_value_orig)+sd(tmp_value_orig)*1.96,code=3,angle=0)
}
dev.off()
