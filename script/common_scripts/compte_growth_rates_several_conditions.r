####CP 24/11/2020: This code to check the behaviour of the growth rates in several conditions to try and understand how can growth rates can relate to extinctions

rm(list=ls())
graphics.off()

source("step_functions.r")

load("no_seed_compet_simu_light.RData")
name_spp=rownames(mean_val)
nspp=length(name_spp)

#Sinking rates and T_opt
tab=read.table(paste("../../param/species_specific_parameters.txt",sep=""),sep=";",dec=".",header=T)
tab=tab[name_spp,]
T_opt=tab$T_opt+273.15+5
names(T_opt)=name_spp
B=tab$Val_b
names(B)=name_spp

####Model I
#InteractionG
tmp_inter_I=read.table("../../output/modelv1.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
tmp_inter_I=as.matrix(tmp_inter_I)
list_inter_I=list(tmp_inter_I,tmp_inter_I) #Interaction list contains coastal and oceanic interactions

####ModelII
list_H_tmp=read.table("../../output/modelv2.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
list_H=list(list_H_tmp,list_H_tmp)
type_inter_tmp=read.table("../../output/modelv2.1/facilitation_and_competition.csv",sep=";",dec=".")
type_inter_II=list(type_inter_tmp,type_inter_tmp)

save(N_array, N_orig, mean_val,min_val,range_val,fac_simu,nb_persistence,id,file = "no_seed_compet_simu_light.RData", envir = .GlobalEnv)

temp=15+273.15

growth_alone_mean_temp=rep(NA,nspp)
growth_alone_opt_temp=rep(NA,nspp)

#Growth at low density
growth_low_density_modelI=rep(NA,nspp)
growth_low_density_modelII=rep(NA,nspp)
for (s in 1:nspp){
	n_t=rep(1,nspp)
	n_t[s]=1

	growth=growth_rate_noMTE_Bissinger(temp,T_opt,B,0.5)
	growth_alone_mean_temp=growth        

	#Model I
	effect_compet=1+list_inter_I[[1]]%*%n_t
	tmp_growth=exp(growth)*n_t/pmax(1,1+list_inter_I[[1]]%*%n_t)
        growth_low_density_modelI[s]=tmp_growth[s] #We don't look at mortality

	#Model II
	for (sp1 in 1:nspp){
        	val1=0
                for(sp2 in 1:nspp){
                	if(list_H[[1]][sp1,sp2]>0){
                        	val1=val1+type_inter_II[[1]][sp1,sp2]*n_t[sp2]/(list_H[[1]][sp1,sp2]+n_t[sp2])
                        }
                }
       		effect_compet[sp1]=1+val1
	}
	tmp_growth=exp(growth)*n_t/effect_compet
        growth_low_density_modelII[s]=tmp_growth[s]
}

#Growth at average density and average temperature
growth=growth_rate_noMTE_Bissinger(temp,T_opt,B,0.5)
#Model I
n_t=mean_val[,1]
effect_compet=1+list_inter_I[[1]]%*%n_t
tmp_growth=exp(growth)*n_t/pmax(1,1+list_inter_I[[1]]%*%n_t)
growth_average_density_modelI=tmp_growth #We don't look at mortality

#Model II
n_t=mean_val[,2]
for (sp1 in 1:nspp){
	val1=0
	for(sp2 in 1:nspp){
        	if(list_H[[1]][sp1,sp2]>0){
                                val1=val1+type_inter_II[[1]][sp1,sp2]*n_t[sp2]/(list_H[[1]][sp1,sp2]+n_t[sp2])
                        }
                }
                effect_compet[sp1]=1+val1
        }
        tmp_growth=exp(growth)*n_t/effect_compet
growth_average_density_modelII=tmp_growth

#Growth at optimal temperature and average density
growth_optimal_temperature_modelI=rep(NA,nspp)
growth_optimal_temperature_modelII=rep(NA,nspp)
for (s in 1:nspp){
	temp=T_opt[s]
        growth=growth_rate_noMTE_Bissinger(temp,T_opt,B,0.5)
	growth_alone_opt_temp[s]=growth[s]

        #Model I
	n_t=mean_val[,1]
        effect_compet=1+list_inter_I[[1]]%*%n_t
        tmp_growth=exp(growth)*n_t/pmax(1,1+list_inter_I[[1]]%*%n_t)
        growth_optimal_temperature_modelI[s]=tmp_growth[s] #We don't look at mortality

        #Model II
	n_t=mean_val[,2]
        for (sp1 in 1:nspp){
                val1=0
                for(sp2 in 1:nspp){
                        if(list_H[[1]][sp1,sp2]>0){
                                val1=val1+type_inter_II[[1]][sp1,sp2]*n_t[sp2]/(list_H[[1]][sp1,sp2]+n_t[sp2])
                        }
                }
                effect_compet[sp1]=1+val1
        }
        tmp_growth=exp(growth)*n_t/effect_compet
        growth_optimal_temperature_modelII[s]=tmp_growth[s]
}

#Growth at average density and bad temperature
temp=30+273.15
growth=growth_rate_noMTE_Bissinger(temp,T_opt,B,0.5)
#Model I
n_t=mean_val[,1]
effect_compet=1+list_inter_I[[1]]%*%n_t
tmp_growth=exp(growth)*n_t/pmax(1,1+list_inter_I[[1]]%*%n_t)
growth_bad_temperature_modelI=tmp_growth #We don't look at mortality

#Model II
n_t=mean_val[,2]
for (sp1 in 1:nspp){
        val1=0
        for(sp2 in 1:nspp){
                if(list_H[[1]][sp1,sp2]>0){
                                val1=val1+type_inter_II[[1]][sp1,sp2]*n_t[sp2]/(list_H[[1]][sp1,sp2]+n_t[sp2])
                        }
                }
                effect_compet[sp1]=1+val1
        }
        tmp_growth=exp(growth)*n_t/effect_compet
growth_bad_temperature_modelII=tmp_growth

id_persistence=seq(length(id)-nb_persistence+1,length(id))
surviv_compet=apply(apply(N_array[id_persistence,"ocean",,,'compet',]==0,c(2,3,4),sum)<nb_persistence,c(1,3),sum)/(length(fac_simu))
surviv_facil=apply(apply(N_array[id_persistence,"ocean",,,'facil',]==0,c(2,3,4),sum)<nb_persistence,c(1,3),sum)/length(fac_simu)

pdf("survival_vs_growth_rates.pdf",width=6,height=6)
set.seed(42)
questionning=c(2,3,9,11) #These species have a high growth rate at low densities but still go extinct
par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(4,4,2,1))
limits=range(c(growth_low_density_modelI,growth_low_density_modelII))
a=jitter(growth_low_density_modelI,amount=0.1)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Realized PCGR at low ab.",pch=16,ylim=c(0,1),ylab="Prob survival",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,1],col="red",pch=16)
text(min(limits)-0.1,1.1,"a",las=2,xpd=NA,font=2)
a=jitter(growth_low_density_modelII,amount=0.1)
points(a,surviv_compet[,2],col="grey",pch=17,cex=1.5)
points(a[questionning],surviv_compet[questionning,2],col="red",pch=16)
legend("topleft",c("Model I","Model II"),pch=c(16,17),col=c("black","grey"),pt.lwd=c(1,1,1.5,1.5),bty="n")

limits=range(c(growth_average_density_modelI,growth_average_density_modelII))
a=jitter(growth_average_density_modelI,amount=0.1)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Realized PCGR at mean ab.",pch=16,ylim=c(0,1),ylab="",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,1],col="red",pch=16)
text(min(limits)-0.1,1.1,"b",las=2,xpd=NA,font=2)
a=jitter(growth_average_density_modelII,amount=0.1)
points(a,surviv_compet[,2],col="grey",pch=17,cex=1.5)
points(a[questionning],surviv_compet[questionning,2],col="red",pch=16)

limits=range(c(growth_optimal_temperature_modelI,growth_optimal_temperature_modelII))
a=jitter(growth_optimal_temperature_modelI,amount=0.1)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Realized PCGR at optimal temp.",pch=16,ylim=c(0,1),ylab="Prob survival",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,1],col="red",pch=16)
text(min(limits)-0.1,1.1,"c",las=2,xpd=NA,font=2)
a=jitter(growth_optimal_temperature_modelII,amount=0.1)
points(a,surviv_compet[,2],col="grey",pch=17,cex=1.5)
points(a[questionning],surviv_compet[questionning,2],col="red",pch=16)

limits=range(c(growth_bad_temperature_modelI,growth_bad_temperature_modelII))
a=jitter(growth_bad_temperature_modelI,amount=0.1)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Realized PCGR at suboptimal temp.",pch=16,ylim=c(0,1),ylab="",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,1],col="red",pch=16)
text(min(limits)-0.1,1.1,"d",las=2,xpd=NA,font=2)
a=jitter(growth_bad_temperature_modelII,amount=0.1)
points(a,surviv_compet[,2],col="grey",pch=17,cex=1.5)
points(a[questionning],surviv_compet[questionning,2],col="red",pch=16)


dev.off()


###Now differetentiating between species and models to make the results clearer
questionning=c(2,3,9,11)
pdf("survival_vs_growth_rates_focus_on_species.pdf",width=9,height=9)
set.seed(42)
par(mfrow=c(2,2),oma=c(1,1,1,1),mar=c(4,4,1,1))
##Model I
limits=range(c(growth_low_density_modelI))
a=jitter(growth_low_density_modelI,amount=0.1)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Growth low density",pch=16,ylim=c(-0.2,1.1),ylab="Prob survival",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,1],col="red",pch=16)
text(a,surviv_compet[,1],name_spp,pos=c(1,3))
text(min(limits)-0.1,1.2,"a",las=2,xpd=NA,font=2)

limits=range(c(growth_average_density_modelI))
a=jitter(growth_average_density_modelI,amount=0.1)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Growth average density",pch=16,ylim=c(-0.2,1),ylab="",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,1],col="red",pch=16)
text(min(limits)-0.1,1.2,"b",las=2,xpd=NA,font=2)
text(a,surviv_compet[,1],name_spp,pos=c(1,3))

limits=range(c(growth_optimal_temperature_modelI))
a=jitter(growth_optimal_temperature_modelI,amount=0.1)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Growth optimal temperature",pch=16,ylim=c(-0.2,1.1),ylab="Prob survival",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,1],col="red",pch=16)
text(min(limits)-0.1,1.2,"c",las=2,xpd=NA,font=2)
text(a,surviv_compet[,1],name_spp,pos=c(1,3))

limits=range(c(growth_bad_temperature_modelI))
a=jitter(growth_bad_temperature_modelI,amount=0.1)
plot(a,surviv_compet[,1],t="p",col="black",xlab="Growth bad temperature",pch=16,ylim=c(-0.2,1.1),ylab="",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,1],col="red",pch=16)
text(min(limits)-0.1,1.2,"d",las=2,xpd=NA,font=2)
text(a,surviv_compet[,1],name_spp,pos=c(1,3))

##Model II
limits=range(c(growth_low_density_modelII))
a=jitter(growth_low_density_modelII,amount=0.1)
plot(a,surviv_compet[,2],t="p",col="black",xlab="Growth low density",pch=16,ylim=c(-0.2,1.1),ylab="Prob survival",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,2],col="red",pch=16)
text(a,surviv_compet[,2],name_spp,pos=c(1,3))
text(min(limits)-0.1,1.2,"a",las=2,xpd=NA,font=2)

limits=range(c(growth_average_density_modelII))
a=jitter(growth_average_density_modelII,amount=0.1)
plot(a,surviv_compet[,2],t="p",col="black",xlab="Growth average density",pch=16,ylim=c(-0.2,1),ylab="",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,2],col="red",pch=16)
text(min(limits)-0.1,1.2,"b",las=2,xpd=NA,font=2)
text(a,surviv_compet[,2],name_spp,pos=c(1,3))

limits=range(c(growth_optimal_temperature_modelII))
a=jitter(growth_optimal_temperature_modelII,amount=0.1)
plot(a,surviv_compet[,2],t="p",col="black",xlab="Growth optimal temperature",pch=16,ylim=c(-0.2,1.1),ylab="Prob survival",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,2],col="red",pch=16)
text(min(limits)-0.1,1.2,"c",las=2,xpd=NA,font=2)
text(a,surviv_compet[,2],name_spp,pos=c(1,3))

limits=range(c(growth_bad_temperature_modelII))
a=jitter(growth_bad_temperature_modelII,amount=0.1)
plot(a,surviv_compet[,2],t="p",col="black",xlab="Growth bad temperature",pch=16,ylim=c(-0.2,1),ylab="",cex=1.5,xlim=limits+c(-0.1,0.1))
points(a[questionning],surviv_compet[questionning,2],col="red",pch=16)
text(min(limits)-0.1,1.2,"d",las=2,xpd=NA,font=2)
text(a,surviv_compet[,2],name_spp,pos=c(1,3))


dev.off()

##Doing survival=f(r_i), i.e. without taking the interactions into account
pdf("intrinsic_growth_rate.pdf",width=15,height=10)
par(mfrow=c(1,2))
a=jitter(growth_alone_mean_temp,amount=0.1)
plot(a,surviv_compet[,2],t="p",col="black",xlab="Intrinsic 15",pch=16,ylim=c(-0.2,1.1),ylab="Prob survival",cex=1.5)
points(a[questionning],surviv_compet[questionning,2],col="red",pch=16)
points(a,surviv_compet[,1],t="p",col="blue",pch=1,cex=1.25)
text(a,surviv_compet[,1],name_spp,pos=c(1),col="blue")
text(a,surviv_compet[,2],name_spp,pos=c(2))
points(a[questionning],surviv_compet[questionning,1],col="red",pch=16)
legend("topleft",c("Model I","Model II"),col=c("blue","black"),pch=16)

a=jitter(growth_alone_opt_temp,amount=0.1)
plot(a,surviv_compet[,2],t="p",col="black",xlab="Intrinsic opt",pch=16,ylim=c(-0.2,1.1),ylab="Prob survival",cex=1.5)
points(a[questionning],surviv_compet[questionning,2],col="red",pch=16)
points(a,surviv_compet[,1],t="p",col="blue",pch=1,cex=1.25)
text(a,surviv_compet[,2],name_spp,pos=c(1))
text(a,surviv_compet[,1],name_spp,pos=c(2),col="blue")
points(a[questionning],surviv_compet[questionning,1],col="red",pch=16)

dev.off()

#questionning=c(2,3,9,11)
#####GUI(3) PRP(11) PSE(9) DIT(2) are the ones with the high growth rate at low density but low proba of survival
