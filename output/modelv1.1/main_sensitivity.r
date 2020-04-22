#################
###CP 18/02/2020: Clean version of the main function
###CP 21/04/2020: Version for main sensitivity
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
####Two options are possible for mean values
pop_table=read.table("../../param/mean_interpolated_abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
pop_table=pop_table[name_spp,]

#Interaction
tmp_inter=read.table("matrix_A_after_calibration.csv",sep=";",dec=".")
name_spp=colnames(tmp_inter)
nspp=length(name_spp)

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
tab=read.table(paste("species_specific_parameters.txt",sep=""),sep=";",dec=".",header=T)
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
#free_param=read.table("free_parameters.txt",sep=";",dec=".",header=T,row.names=1)
free_param=read.table("free_parameters_with_lowburial_testmortality.txt",sep=";",dec=".",header=T,row.names=1)
list_simulation=matrix(NA,nrow=(ncol(free_param)-1)*nrow(free_param),ncol=nrow(free_param))
colnames(list_simulation)=rownames(free_param)
rownames(list_simulation)=1:nrow(list_simulation)

###############Store statistics
tab_summary=matrix(NA,nrow=nrow(list_simulation)+1,ncol=5)
rownames(tab_summary)=1:(nrow(list_simulation)+1)
colnames(tab_summary)=c("Abundance","Amplitude","Phenology","Diff","Persistence")

tab_coast=N_original_set[,1,]
final_summary=summary_statistics(pop_table,tab_pheno,tab_coast,nb_year)
tab_summary[nrow(list_simulation)+1,1:3]=final_summary
tab_summary[nrow(list_simulation)+1,4]=sum(final_summary)
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
		tab_summary[nb_simu,1:3]=final_summary
		if(sum(tab_coast==0)>0){ #One species has died
			tab_summary[nb_simu,4]=0
		}else{
			tab_summary[nb_simu,4]=sum(final_summary)
		}
		tab_summary[nb_simu,5]=sum(tab_coast[nrow(tab_coast),]>0)
}
}

id=seq(n_iter-365,n_iter)
mean_tot_original=mean(log10(apply(N_original_set[id,"coast",],1,sum)))

pdf("mean_abundance_with_free_parameters_with_lossrate.pdf",width=15)
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

