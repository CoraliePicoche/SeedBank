#################
###CP 18/02/2020: Clean version of the main function
###CP 21/04/2020: Version for main sensitivity
###CP 15/05/2020: Corrected a bug in amplitude plot (plotted abundance instead) + corrected the way average and amplitude were taken into account (before, was highly dependent on only the most abundant species)
#################

rm(list=ls())
graphics.off()
set.seed(40) #Warning: results and even feasibility can be highly seed-sensitive

source("step_functions.r")
source("../../script/summary_statistics.r")

n_iter=1000
nb_year=2

#Data to use (Auger)
a=as.matrix(read.table("../../param/reconstructed_temperature_Auger.txt", row.names=1,header=T,sep=";",dec="."))
temp_model=a[1:n_iter]
#Used for summary statistics
tab_pheno=read.table("../../param/generalist_specialist_spp_added_amplitude_season.csv",sep=";",dec=".",header=T,row.names=1)
name_spp=rownames(tab_pheno)
nspp=length(name_spp)
####Two options are possible for mean values
pop_table=read.table("../../param/mean_interpolated_abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
pop_table=pop_table[name_spp,]

#Interaction
tmp_inter=read.table("interaction_matrix_after_calibration.csv",sep=";",dec=".")

#Fixed parameters: golden parameter set
tab=read.table("simu.csv",sep=";",dec=".",header=T)
cyst_mortality=as.numeric(as.character(tab[tab[,1]=="cyst_mortality",2]))
cyst_burial=as.numeric(as.character(tab[tab[,1]=="cyst_burial",2]))
k_coast2ocean=as.numeric(as.character(tab[tab[,1]=="k_coast2ocean",2]))
list_inter=list(as.matrix(tmp_inter),k_coast2ocean*as.matrix(tmp_inter)) #Interaction list contains coastal and oceanic interactions
morta=as.numeric(as.character(tab[tab[,1]=="loss_rate",2]))
e=as.numeric(as.character(tab[tab[,1]=="exchange",2]))
germination=as.numeric(as.character(tab[tab[,1]=="germination",2]))
resuspension=as.numeric(as.character(tab[tab[,1]=="resuspension",2]))
correct=as.numeric(as.character(tab[tab[,1]=="gain",2]))
growth_model=as.character(tab[tab[,1]=="growth_model",2])
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


###############Golden_set
#Initialize
N_original_set=array(NA,dim=c(n_iter,3,nspp),dimnames=list(NULL,c("coast","ocean","seed"),name_spp))
#effect_compet=array(NA,dim=c(n_iter,2,nspp),dimnames=list(NULL,c("coast","ocean"),name_spp))
N_original_set[1,,]=rep(10^3,nspp*3)

##Run
for(t in 1:(n_iter-1)){
                var_tmp=step1(N_original_set[t,,],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_original_set[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}

#Parameters to move
free_param=read.table("sensitivity_parameters_woutcystmortality.txt",sep=";",dec=".",header=T,row.names=1)
list_simulation=matrix(NA,nrow=(ncol(free_param)-1)*nrow(free_param),ncol=nrow(free_param))
colnames(list_simulation)=rownames(free_param)
rownames(list_simulation)=1:nrow(list_simulation)

###############Store statistics
tab_summary=matrix(NA,nrow=nrow(list_simulation)+1,ncol=5)
rownames(tab_summary)=1:(nrow(list_simulation)+1)
colnames(tab_summary)=c("Abundance","Amplitude","Phenology","Diff","Persistence")

tab_coast=N_original_set[,1,]
final_summary=summary_statistics(pop_table,tab_pheno,tab_coast,nb_year)
tab_summary[nrow(list_simulation)+1,1:3]=c(sqrt(sum(final_summary[[1]])/nspp),sqrt(sum(final_summary[[2]])/nspp),sqrt(sum(final_summary[[3]])/nspp))
tab_summary[nrow(list_simulation)+1,4]=sum(tab_summary[nrow(list_simulation)+1,1:3])
tab_summary[nrow(list_simulation)+1,5]=sum(tab_coast[nrow(tab_coast),]>0)

#Initialize
N_sensitivity=array(NA,dim=c(n_iter,3,nspp,nrow(list_simulation)),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,1:nrow(list_simulation)))
N_sensitivity[1,,,]=rep(10^3,nspp*3)

nb_simu=0
for(param_to_move in rownames(free_param)){
	for(value in free_param[param_to_move,c("Value_min","Value_max")]){
		nb_simu=nb_simu+1
		#Keep other parameters
		all_others=as.numeric(as.character(free_param[,"Value_used"]))
		names(all_others)=rownames(free_param)
		for(a in 1:length(all_others)){
			assign(names(all_others)[a],all_others[a])
		}
		list_simulation[nb_simu,names(all_others)]=all_others
		value=as.numeric(as.character(value))
		names(value)=param_to_move
		assign(param_to_move,value)
		list_simulation[nb_simu,param_to_move]=value

		rownames(list_simulation)[nb_simu]=paste(param_to_move,value,sep="_")
		
		#Composite phenomena
		M=cyst_mortality+cyst_burial
		Gamma=resuspension*germination
		S=S_max*tab$S
		
		print(list_simulation[nb_simu,])

##Run
N_simu=N_sensitivity[,,,nb_simu]
for(t in 1:(n_iter-1)){
	        var_tmp=step1(N_simu[t,,],list_inter,temp_model[t],M,morta,a_d,T_opt,B,threshold)
		Ntmp=var_tmp[[1]]
        	N_simu[t+1,,]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
		N_sensitivity[,,,nb_simu]=N_simu[,,]
		tab_coast=N_simu[,1,]
		final_summary=summary_statistics(pop_table,tab_pheno,tab_coast,nb_year)
		tab_summary[nb_simu,1:3]=c(sqrt(sum(final_summary[[1]])/nspp),sqrt(sum(final_summary[[2]])/nspp),sqrt(sum(final_summary[[3]])/nspp))
		if(sum(tab_coast==0)>0){ #One species has died
			tab_summary[nb_simu,4]=0
		}else{
			tab_summary[nb_simu,4]=sum(tab_summary[nb_simu,1:3])
		}
		tab_summary[nb_simu,5]=sum(tab_coast[nrow(tab_coast),]>0)
}
}

id=seq(n_iter-365,n_iter)
mean_tot_original=mean(log10(apply(N_original_set[id,"coast",],1,sum)))

pdf("mean_abundance_sensitivity_v1.pdf",width=15)
analyses=rownames(list_simulation)
plot(0,0,t="n",xlim=c(1,nrow(free_param)),ylim=c(3.5,5.5),xaxt="n",ylab="Average total abundance",xlab="")
axis(1,labels=rownames(free_param),at=1:nrow(free_param))
mtext(c(all_others),1,line=3,at=1:nrow(free_param),cex=0.8)
abline(h=mean_tot_original)
l=0
val_text=c()
at_val_text=c()
for(param_to_move in rownames(free_param)){
	l=l+1
	id_param=grep(paste("^",param_to_move,sep=""),analyses)
	space=0.5/(2*length(id_param))
	seq_space=seq(-space,space,length.out=length(id_param))
	for(i in 1:length(id_param)){
		mean_tot=mean(log10(apply(N_sensitivity[id,"coast",,id_param[i]],1,sum)))
		print(mean_tot)
		points(l+seq_space[i],mean_tot,t="p",pch=16,col="black")
		if(tab_summary[id_param[i],4]==0){
			points(l+seq_space[i],mean_tot,t="p",pch=16,col="red")
		}
		at_val_text=c(at_val_text,l+seq_space[i])
		tmp_text=strsplit(analyses[id_param[i]],"_")
		val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
	}
}	
mtext(val_text,1,line=2,at=at_val_text,cex=0.8)
dev.off()

mean_tot_original=log10(apply(N_original_set[id,"coast",],2,mean))
amplitude_tot_original=apply(log10(apply(N_original_set[id,"coast",],2,range)),2,diff)
analyses=rownames(list_simulation)


####WARNING: this code is not flexible enough to take into account more than 2 values per parameter (min and max). If we want to try other values, we will need to increase the 2nd dimension of the two matrices
matrix_mean=array(NA,dim=c(nrow(free_param),2,length(mean_tot_original)),dimnames=list(rownames(free_param),c("min","max"),names(mean_tot_original)))
matrix_amplitude=array(NA,dim=c(nrow(free_param),2,length(amplitude_tot_original)),dimnames=list(rownames(free_param),c("min","max"),names(amplitude_tot_original)))

for(param_to_move in rownames(free_param)){
        id_param=grep(paste("^",param_to_move,sep=""),analyses) #This will help issue an error if we have more than 2 values per parameter
        for(i in 1:length(id_param)){
                mean_sens=log10(apply(N_sensitivity[id,"coast",,id_param[i]],2,mean))
                tmp_mean=100*(mean_tot_original-mean_sens)/mean_tot_original
		matrix_mean[param_to_move,i,]=tmp_mean
		amp_sens=apply(log10(apply(N_sensitivity[id,"coast",,id_param[[i]]],2,range)),2,diff)
		tmp_amplitude=100*(amplitude_tot_original-amp_sens)/amplitude_tot_original
		matrix_amplitude[param_to_move,i,]=tmp_amplitude

                #Now, if the species disappeared, we set their stats to 0
                matrix_mean[param_to_move,i,N_sensitivity[id[length(id)],"coast",,id_param[i]]==0]=NA
                matrix_amplitude[param_to_move,i,N_sensitivity[id[length(id)],"coast",,id_param[i]]==0]=NA

	}
}

pdf("mean_abundance_amplitude_sensitivity.pdf",width=17)
par(mfrow=c(1,2))
plot(0,0,t="n",xlim=c(0.5,nrow(free_param)+0.5),ylim=c(-7,4),xaxt="n",ylab="%Change in average abundance",xlab="")
axis(1,labels=rownames(free_param),at=1:nrow(free_param))
abline(h=0)
mtext(c(all_others),1,line=3,at=1:nrow(free_param),cex=0.9)
l=0
val_text=c()
at_val_text=c()
for(param_to_move in rownames(free_param)){
        l=l+1
        id_param=grep(paste("^",param_to_move,sep=""),analyses)
        space=0.5/(2*length(id_param))
        seq_space=seq(-space,space,length.out=length(id_param))
        for(i in 1:length(id_param)){
		per_change=mean(matrix_mean[param_to_move,i,]*is.finite(matrix_mean[param_to_move,i,]),na.rm=T)
               	rect(l+0.75*seq_space[i],min(c(0,per_change)),l+1.25*seq_space[i],max(c(per_change,0)),pch=16,col="grey")
                if(tab_summary[id_param[i],4]==0){
                        points(l+seq_space[i],per_change+(0.1*11)*sign(per_change),t="p",pch=16,col="red")
                }
                at_val_text=c(at_val_text,l+seq_space[i])
                tmp_text=strsplit(analyses[id_param[i]],"_")
                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
}
mtext(val_text,1,line=2,at=at_val_text,cex=0.9)

plot(0,0,t="n",xlim=c(0.5,nrow(free_param)+0.5),ylim=c(-40,25),xaxt="n",ylab="%Change in average amplitude",xlab="")
axis(1,labels=rownames(free_param),at=1:nrow(free_param))
abline(h=0)
mtext(c(all_others),1,line=3,at=1:nrow(free_param),cex=0.8)
l=0
val_text=c()
at_val_text=c()
for(param_to_move in rownames(free_param)){
        l=l+1
        id_param=grep(paste("^",param_to_move,sep=""),analyses)
        space=0.5/(2*length(id_param))
        seq_space=seq(-space,space,length.out=length(id_param))
        for(i in 1:length(id_param)){
                per_change=mean(matrix_amplitude[param_to_move,i,],na.rm=T)
                rect(l+0.75*seq_space[i],min(c(0,per_change)),l+1.25*seq_space[i],max(c(per_change,0)),pch=16,col="grey")
                if(tab_summary[id_param[i],4]==0){
                        points(l+seq_space[i],per_change+(0.1*65)*sign(per_change),t="p",pch=16,col="red")
                }
                at_val_text=c(at_val_text,l+seq_space[i])
                tmp_text=strsplit(analyses[id_param[i]],"_")
                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
}
mtext(val_text,1,line=2,at=at_val_text,cex=0.8)
dev.off()

#ydelim_mean=range(c(matrix_mean*is.finite(matrix_mean)),na.rm=T)
#ydelim_amp=range(c(matrix_amplitude*is.finite(matrix_amplitude)),na.rm=T)

pdf("species_abundance_amplitude_sensitivity.pdf",width=12,height=10)
dim=c(4,4,3)
s=0
for(d in 1:length(dim)){
par(mfrow=c(dim[d],2),mar=c(5,4,2,1))
sp=s+1
for(s in seq(sp,sp+dim[d]-1)){
	ydelim_mean=range(c(matrix_mean[,,s]*is.finite(matrix_mean[,,s])),na.rm=T)
	plot(0,0,t="n",xlim=c(0.5,nrow(free_param)+0.5),ylim=ydelim_mean,xaxt="n",ylab="%Change in abundance",xlab="",main=dimnames(matrix_mean)[[3]][s])
	if((s-sp+1)==dim[d]){
	axis(1,labels=rownames(free_param),at=1:nrow(free_param))
	mtext(c(all_others),1,line=3,at=1:nrow(free_param),cex=0.6)
	mtext(val_text,1,line=2,at=at_val_text,cex=0.6)
	}
	abline(h=0)
	l=0
	val_text=c()
	at_val_text=c()
	for(param_to_move in rownames(free_param)){
        	l=l+1
	        id_param=grep(paste("^",param_to_move,sep=""),analyses)
        	space=0.5/(2*length(id_param))
	        seq_space=seq(-space,space,length.out=length(id_param))
        	for(i in 1:length(id_param)){
                	per_change=matrix_mean[param_to_move,i,s]
	                rect(l+0.75*seq_space[i],min(c(0,per_change),na.rm=T),l+1.25*seq_space[i],max(c(per_change,0),na.rm=T),pch=16,col="grey")
	                if(is.infinite(per_change)|is.na(per_change)){
        	                points(l+seq_space[i],0,t="p",pch=16,col="red")
                	}
                at_val_text=c(at_val_text,l+seq_space[i])
                tmp_text=strsplit(analyses[id_param[i]],"_")
                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
	}
	ydelim_amp=range(c(matrix_amplitude[,,s]*is.finite(matrix_amplitude[,,s])),na.rm=T)
	plot(0,0,t="n",xlim=c(0.5,nrow(free_param)+0.5),ylim=ydelim_amp,xaxt="n",ylab="%Change in amplitude",xlab="",main=dimnames(matrix_mean)[[3]][s])
	if((s-sp+1)==dim[d]){
	axis(1,labels=rownames(free_param),at=1:nrow(free_param))
	mtext(c(all_others),1,line=3,at=1:nrow(free_param),cex=0.6)
	mtext(val_text,1,line=2,at=at_val_text,cex=0.6)
	}
	abline(h=0)
	l=0
	val_text=c()
	at_val_text=c()
	for(param_to_move in rownames(free_param)){
        	l=l+1
      	 	id_param=grep(paste("^",param_to_move,sep=""),analyses)
        	space=0.5/(2*length(id_param))
        	seq_space=seq(-space,space,length.out=length(id_param))
        	for(i in 1:length(id_param)){
                	per_change=matrix_amplitude[param_to_move,i,s]
                	rect(l+0.75*seq_space[i],min(c(0,per_change),na.rm=T),l+1.25*seq_space[i],max(c(per_change,0),na.rm=T),pch=16,col="grey")
                	if(is.infinite(per_change)|is.na(per_change)){
                        	points(l+seq_space[i],0,t="p",pch=16,col="red")
                	}
                	at_val_text=c(at_val_text,l+seq_space[i])
                	tmp_text=strsplit(analyses[id_param[i]],"_")
                	val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        	}
	}
}

}
dev.off()

pdf("mean_amplitude_sensitivity_v1.pdf",width=15)
analyses=rownames(list_simulation)
plot(0,0,t="n",xlim=c(1,nrow(free_param)),ylim=c(1,2),xaxt="n",ylab="Average amplitude of total community",xlab="")
axis(1,labels=rownames(free_param),at=1:nrow(free_param))
mtext(c(all_others),1,line=3,at=1:nrow(free_param),cex=0.8)
abline(h=mean_tot_original)
l=0
val_text=c()
at_val_text=c()
for(param_to_move in rownames(free_param)){
        l=l+1
        id_param=grep(paste("^",param_to_move,sep=""),analyses)
        space=0.5/(2*length(id_param))
        seq_space=seq(-space,space,length.out=length(id_param))
        for(i in 1:length(id_param)){
                mean_tot=diff(range(log10(apply(N_sensitivity[id,"coast",,id_param[i]],1,sum))))
                print(mean_tot)
                points(l+seq_space[i],mean_tot,t="p",pch=16,col="black")
                if(tab_summary[id_param[i],4]==0){
                        points(l+seq_space[i],mean_tot,t="p",pch=16,col="red")
                }
                at_val_text=c(at_val_text,l+seq_space[i])
                tmp_text=strsplit(analyses[id_param[i]],"_")
                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
}
mtext(val_text,1,line=2,at=at_val_text,cex=0.8)
dev.off()

pdf("summary_statistics_for_sensitivity.pdf",width=15)
par(mfrow=c(2,2))
analyses=rownames(list_simulation)
for(summary in 1:(ncol(tab_summary)-1)){
	plot(0,0,t="n",xlim=c(1,nrow(free_param)),ylim=range(tab_summary[,summary]),xaxt="n",ylab=colnames(tab_summary)[summary],xlab="")
	axis(1,labels=rownames(free_param),at=1:nrow(free_param))
	mtext(c(all_others),1,line=3,at=1:nrow(free_param),cex=0.8)
	abline(h=tab_summary[nrow(tab_summary),summary])
	l=0
	val_text=c()
	at_val_text=c()
	for(param_to_move in rownames(free_param)){
        	l=l+1
	        id_param=grep(paste("^",param_to_move,sep=""),analyses)
        	space=0.5/(2*length(id_param))
	        seq_space=seq(-space,space,length.out=length(id_param))
        	for(i in 1:length(id_param)){
			if(tab_summary[id_param[i],4]==0){
				text(l+seq_space[i],tab_summary[id_param[i],summary],labels=tab_summary[id_param[i],5],col="red")
			}else{
                		points(l+seq_space[i],tab_summary[id_param[i],summary],col="black",t="p",pch=16)
			}
	                at_val_text=c(at_val_text,l+seq_space[i])
        	        tmp_text=strsplit(analyses[id_param[i]],"_")
	                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
}
mtext(val_text,1,line=2,at=at_val_text,cex=0.8)
}
dev.off()

write.table(list_simulation,paste("list_simulation_sensitivity.csv",sep=""),sep=";",dec=".")
write.table(tab_summary,paste("list_statistics_sensitivity.csv",sep=""),sep=";",dec=".")
