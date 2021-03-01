#############
# 20/05/20 CP: attempt at a general sensitivity script
##############

graphics.off()
rm(list=ls())

source("step_functions.r")
source("../summary_statistics.r")

doyouload=TRUE #TRUE if taking previous simulations to draw new plots, FALSE if we want to re-run the simulations

if(!doyouload){

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
Gamma=resuspension*germination
S=S_max*tab$S

####Model I
#Interaction
tmp_inter_I=read.table("../../output/modelv1.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
tmp_inter_I=as.matrix(tmp_inter_I)
list_inter_I=list(tmp_inter_I,k_coast2ocean*tmp_inter_I) #Interaction list contains coastal and oceanic interactions

####ModelII
list_H_tmp=read.table("../../output/modelv2.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
list_H_II=list(list_H_tmp,list_H_tmp)
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
                
		var_tmp=step1_modelII(N_original_set[t,,,'modelII'],list_H_II,type_inter_II,temp_model[t],M,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
                N_original_set[t+1,,,'modelII']=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}

#Parameters to move
free_param=read.table("sensitivity_parameters_woutcystmortality.txt",sep=";",dec=".",header=T,row.names=1)

list_simulation=matrix(NA,nrow=(ncol(free_param)-1)*nrow(free_param),ncol=nrow(free_param))
colnames(list_simulation)=rownames(free_param)
rownames(list_simulation)=1:nrow(list_simulation)

###############Store statistics
tab_summary=array(NA,dim=c(nrow(list_simulation)+1,6,2),dimnames=list(1:(nrow(list_simulation)+1),c("Abundance","Amplitude","Phenology","Diff","Persistence_coast","Persistence_ocean"),c("modelI","modelII")))

for(m in 1:2){
tab_coast=N_original_set[,cpt,,m]
final_summary=summary_statistics(pop_table,tab_pheno,tab_coast,nb_year)
tab_summary[nrow(list_simulation)+1,1:3,m]=c(sqrt(sum(final_summary[[1]]^2)/nspp),sqrt(sum(final_summary[[2]]^2)/nspp),sqrt(sum(final_summary[[3]]^2)/nspp))
tab_summary[nrow(list_simulation)+1,4,m]=sum(tab_summary[nrow(list_simulation)+1,1:3,m])

tab_summary[nrow(list_simulation)+1,"Persistence_coast",m]=sum(apply(N_original_set[id_persistence,"coast",,m]==0,2,sum)<nb_persistence)
tab_summary[nrow(list_simulation)+1,"Persistence_ocean",m]=sum(apply(N_original_set[id_persistence,"ocean",,m]==0,2,sum)<nb_persistence)


}

#Initialize
N_sensitivity=array(NA,dim=c(length(id),3,nspp,nrow(list_simulation),2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,1:nrow(list_simulation),c("modelI","modelII")))

nb_simu=0
for(param_to_move in rownames(free_param)){
#for(param_to_move in rownames(free_param)[4]){
        for(value in free_param[param_to_move,c("Value_min","Value_max")]){
#        for(value in free_param[param_to_move,c("Value_max")]){
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
N_simu=array(NA,dim=c(n_iter,3,nspp,2),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,c("modelI","modelII")))
N_simu[1,,,]=10^3
for(t in 1:(n_iter-1)){

		var_tmp=step1_modelI(N_simu[t,,,1],list_inter_I,temp_model[t],M,morta,a_d,T_opt,B,threshold)
                Ntmp=var_tmp[[1]]
                N_simu[t+1,,,1]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)

                var_tmp=step1_modelII(N_simu[t,,,2],list_H_II,type_inter_II,temp_model[t],M,morta,a_d,T_opt,B)
                Ntmp=var_tmp[[1]]
                N_simu[t+1,,,2]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),e)
}
                N_sensitivity[,,,nb_simu,]=N_simu[id,,,]
		for(m in 1:2){
                tab_coast=N_simu[,cpt,,m]
                if(sum(tab_coast==0)>0){ #One species has died
                        tab_summary[nb_simu,1:3,m]=NA
                        tab_summary[nb_simu,4,m]=0
                }else{
                        final_summary=summary_statistics(pop_table,tab_pheno,tab_coast,nb_year)
                        tab_summary[nb_simu,1:3,m]=c(sqrt(sum(final_summary[[1]]^2)/nspp),sqrt(sum(final_summary[[2]]^2)/nspp),sqrt(sum(final_summary[[3]]^2)/nspp))
                        tab_summary[nb_simu,4,m]=sum(tab_summary[nb_simu,1:3,m])
                }
		tab_summary[nb_simu,"Persistence_coast",m]=sum(apply(N_simu[id_persistence,"coast",,m]==0,2,sum)<nb_persistence)
		tab_summary[nb_simu,"Persistence_ocean",m]=sum(apply(N_simu[id_persistence,"ocean",,m]==0,2,sum)<nb_persistence)

		}
#		if(nb_simu==8){stop()}
}
}
analyses=rownames(list_simulation)

save(list = ls(all.names = TRUE), file = "main_sensitivity_simu.RData", envir = .GlobalEnv)
}else{ #end if doyouload
load("main_sensitivity_simu.RData")
}

#Difference between abundance on the coast and on the ocean
tmp_var=N_original_set[id,c("coast"),,]-N_original_set[id,c("ocean"),,]
diff_tot_original=log10(abs(apply(tmp_var,c(2,3),mean)))

#Abundance total in the ocean
tmp_var=apply(N_original_set[id,"ocean",,],c(1,3),sum)
ab_tot_original=log10(apply(tmp_var,2,mean))

for(cpt in c("coast")){
mean_tot_original=log10(apply(N_original_set[id,cpt,,],c(2,3),mean))
amplitude_tot_original=apply(log10(apply(N_original_set[id,cpt,,],c(2,3),range)),2,diff)

####WARNING: this code is not flexible enough to take into account more than 2 values per parameter (min and max). If we want to try other values, we will need to increase the 2nd dimension of the two matrices
matrix_mean=array(NA,dim=c(nrow(free_param),2,nrow(mean_tot_original),2),dimnames=list(rownames(free_param),c("min","max"),rownames(mean_tot_original),c("modelI","modelII")))
matrix_amplitude=array(NA,dim=c(nrow(free_param),2,ncol(amplitude_tot_original),2),dimnames=list(rownames(free_param),c("min","max"),colnames(amplitude_tot_original),c("modelI","modelII")))
matrix_diff=array(NA,dim=c(nrow(free_param),2,ncol(amplitude_tot_original),2),dimnames=list(rownames(free_param),c("min","max"),colnames(amplitude_tot_original),c("modelI","modelII")))
matrix_abundance_tot=array(NA,dim=c(nrow(free_param),2,ncol(amplitude_tot_original),2),dimnames=list(rownames(free_param),c("min","max"),colnames(amplitude_tot_original),c("modelI","modelII")))

for(param_to_move in rownames(free_param)){
        id_param=grep(paste("^",param_to_move,sep=""),analyses) #This will help issue an error if we have more than 2 values per parameter
        for(i in 1:length(id_param)){
                mean_sens=log10(apply(N_sensitivity[,cpt,,id_param[i],],c(2,3),mean)) #mean_sens contains the log10 of average abundance per species, per model
                tmp_mean=100*(mean_tot_original-mean_sens)/mean_tot_original
                matrix_mean[param_to_move,i,,]=tmp_mean

                amp_sens=apply(log10(apply(N_sensitivity[,cpt,,id_param[i],],c(2,3),range)),2,diff) #amp_sens contains the log-ratio maximum/minimum per species, per model
                tmp_amplitude=100*(amplitude_tot_original-amp_sens)/amplitude_tot_original
                matrix_amplitude[param_to_move,i,,]=t(tmp_amplitude)

                ##Now, if the species disappeared, we set their stats to 0
		#for(m in 1:2){
                #matrix_mean[param_to_move,i,N_sensitivity[id[length(id)],cpt,,id_param[i],m]==0,m]=NA
                #matrix_amplitude[param_to_move,i,N_sensitivity[id[length(id)],cpt,,id_param[i],m]==0,m]=NA
		#}
		tmp_var=N_sensitivity[,c("coast"),,id_param[i],]-N_sensitivity[,c("ocean"),,id_param[i],]
		matrix_diff[param_to_move,i,,]=log10(abs(apply(tmp_var,c(2,3),mean,na.rm=T)))

		ab_tot=log10(apply(apply(N_sensitivity[,c("ocean"),,id_param[i],],c(1,3),sum,na.rm=T),2,mean))
		matrix_abundance_tot[param_to_move,i,,]=100*(ab_tot_original-ab_tot)/ab_tot_original

        }
}

let_panels=c("a","b")

pdf(paste("mean_abundance_amplitude_sensitivity.pdf",sep=""),width=12,height=10)
par(mfrow=c(2,1),mar=c(3,5,3,1))
plot(0,0,t="n",xlim=c(0.75,nrow(free_param)+0.25),ylim=c(-5,2),xaxt="n",ylab="%Diff log10(mean abundance)",xlab="",cex.lab=1.65,cex.axis=1.65)
#axis(1,labels=rownames(free_param),at=1:nrow(free_param),cex.axis=1.525)
abline(h=0)
#mtext(c(all_others),1,line=3.5,at=1:nrow(free_param),cex=1.3)
mtext("a",3,line=0.5,at=0.6,cex=1.75,font=2)
l=0
val_text=c()
at_val_text=c()
for(param_to_move in rownames(free_param)){
        l=l+1
        id_param=grep(paste("^",param_to_move,sep=""),analyses)
        space=0.5/(2*length(id_param))
        space_txt=0.8/(2*length(id_param))
        seq_space=seq(-space,space,length.out=length(id_param))
        seq_space_txt=seq(-space_txt,space_txt,length.out=length(id_param))
        for(i in 1:length(id_param)){

		for(mod in 1:2){
                per_change=mean(matrix_mean[param_to_move,i,,mod]*is.finite(matrix_mean[param_to_move,i,,mod]),na.rm=T)
		if(mod==1){
			colm="lightgrey"
			idplus=idplus_text=0
		}else{
			colm="darkgrey"
			idplus=0.08
			idplus_text=0.
		}
                rect(l+0.75*seq_space[i]+idplus,min(c(0,per_change)),l+1.25*seq_space[i]+idplus,max(c(per_change,0)),col=colm)
		print(analyses[i])
		print(tab_summary[id_param[i],,mod])
                if(tab_summary[id_param[i],4,mod]==0){
			if((tab_summary[id_param[i],"Persistence_coast",mod]<nspp)|(tab_summary[id_param[i],"Persistence_ocean",mod]<nspp)){
                        #text(l+seq_space[i]+idplus,per_change+(0.1*9)*sign(per_change),tab_summary[id_param[i],"Persistence",mod],col="red")
                        text(l+seq_space[i]+idplus,2,tab_summary[id_param[i],"Persistence_ocean",mod],col="red",cex=1.5)
#                        text(l+seq_space[i]+idplus_text,3-0.4,tab_summary[id_param[i],"Persistence_coast",mod],col="blue")
			}else{
                        #points(l+seq_space[i]+idplus,per_change+(0.1*11)*sign(per_change),t="p",pch=16,col="red",cex=0.5)
                        points(l+seq_space[i]+idplus,2,t="p",pch=16,col="red",cex=1)
			}
                }
		}
                at_val_text=c(at_val_text,l+seq_space_txt[i])
                tmp_text=strsplit(analyses[id_param[i]],"_")
                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
}
#mtext(val_text,1,line=2.25,at=at_val_text,cex=1.3)
legend("bottomleft",c("Model I","Model II"),col=c("lightgrey","darkgrey"),pch=16,bty="n",cex=1.75)

par(mar=c(4.5,5,1.5,1))
plot(0,0,t="n",xlim=c(0.75,nrow(free_param)+0.25),ylim=c(-40,30),xaxt="n",ylab="%Diff log.amplitude",xlab="",cex.axis=1.65,cex.lab=1.65)
axis(1,labels=rownames(free_param),at=1:nrow(free_param),cex.axis=1.525)
abline(h=0)
mtext(c(all_others),1,line=3.5,at=1:nrow(free_param),cex=1.4)
mtext("b",3,line=0.5,at=0.6,cex=1.75,font=2)
l=0
val_text=c()
at_val_text=c()
for(param_to_move in rownames(free_param)){
        l=l+1
        id_param=grep(paste("^",param_to_move,sep=""),analyses)
        space=0.5/(2*length(id_param))
        seq_space=seq(-space,space,length.out=length(id_param))
        for(i in 1:length(id_param)){
		for(mod in 1:2){
		print(paste(analyses[id_param],mod))
                per_change=mean(matrix_amplitude[param_to_move,i,,mod]*is.finite(matrix_mean[param_to_move,i,,mod]),na.rm=T)
		print(per_change)
                if(mod==1){
                        colm="lightgrey"
                        idplus=0
                }else{
                        colm="darkgrey"
                        idplus=0.08
                } 
                rect(l+0.75*seq_space[i]+idplus,min(c(0,per_change)),l+1.25*seq_space[i]+idplus,max(c(per_change,0)),col=colm)
                #if(tab_summary[id_param[i],4,mod]==0){
			#if(tab_summary[id_param[i],"Persistence",mod]<nspp){
                        #text(l+seq_space[i]+idplus,per_change+(0.1*9)*sign(per_change),tab_summary[id_param[i],"Persistence",mod],col="red")
                        #text(l+seq_space[i]+idplus,30,tab_summary[id_param[i],"Persistence",mod],col="red")
			#}else{
                        #points(l+seq_space[i]+idplus,per_change+(0.1*70)*sign(per_change),t="p",pch=16,col="red",cex=0.5)
                        #points(l+seq_space[i]+idplus,30,t="p",pch=16,col="red",cex=0.5)
			#}
                #}
		}
                at_val_text=c(at_val_text,l+seq_space_txt[i])
                tmp_text=strsplit(analyses[id_param[i]],"_")
                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
}
mtext(val_text,1,line=2.25,at=at_val_text,cex=1.3)
dev.off()
} #end of  cpt="coast

pdf("sensitivity_totabundance_ocean_diff_with_coast.pdf",width=13)
par(mfrow=c(1,2))
plot(0,0,t="n",xlim=c(0.75,nrow(free_param)+0.25),ylim=c(-4,5),xaxt="n",ylab="%Change in log10(mean(tot(abundance))) ocean",xlab="")
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

                for(mod in 1:2){
                per_change=mean(matrix_abundance_tot[param_to_move,i,,mod]*is.finite(matrix_abundance_tot[param_to_move,i,,mod]),na.rm=T)
                if(mod==1){
                        colm="lightgrey"
                        idplus=0
                }else{
                        colm="darkgrey"
                        idplus=0.08
                }
                rect(l+0.75*seq_space[i]+idplus,min(c(0,per_change)),l+1.25*seq_space[i]+idplus,max(c(per_change,0)),col=colm)
                if(tab_summary[id_param[i],4,mod]==0){
			if(tab_summary[id_param[i],"Persistence_ocean",mod]<nspp){
                        text(l+seq_space[i]+idplus,per_change+(0.1*9)*sign(per_change),tab_summary[id_param[i],"Persistence_ocean",mod],col="red")
			}else{
                        points(l+seq_space[i]+idplus,per_change+(0.1*9)*sign(per_change),t="p",pch=16,col="red",cex=0.5)
			}
                }
                }
                at_val_text=c(at_val_text,l+seq_space[i])
                tmp_text=strsplit(analyses[id_param[i]],"_")
                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
}
mtext(val_text,1,line=2,at=at_val_text,cex=0.9)

plot(0,0,t="n",xlim=c(0.75,nrow(free_param)+0.25),ylim=c(-3,4),xaxt="n",ylab="log10(|abundance coast-abundance_ocean|)",xlab="")
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

                for(mod in 1:2){
                per_change=mean(matrix_diff[param_to_move,i,,mod]*is.finite(matrix_diff[param_to_move,i,,mod]),na.rm=T)
                if(mod==1){
                        colm="lightgrey"
                        idplus=0
                }else{
                        colm="darkgrey"
                        idplus=0.08
                }
                rect(l+0.75*seq_space[i]+idplus,min(c(0,per_change)),l+1.25*seq_space[i]+idplus,max(c(per_change,0)),col=colm)
                if(tab_summary[id_param[i],4,mod]==0){
                        points(l+seq_space[i]+idplus,per_change+(0.1*1)*sign(per_change),t="p",pch=16,col="red",cex=0.5)
                }
                }
                at_val_text=c(at_val_text,l+seq_space[i])
                tmp_text=strsplit(analyses[id_param[i]],"_")
                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
}
mtext(val_text,1,line=2,at=at_val_text,cex=0.9)

dev.off()

#######################################################################
#Now, same thing for the seed bank
cpt="seed"
mean_tot_original=log10(apply(N_original_set[id,cpt,,],c(2,3),mean))
amplitude_tot_original=apply(log10(apply(N_original_set[id,cpt,,],c(2,3),range)),2,diff)

####WARNING: this code is not flexible enough to take into account more than 2 values per parameter (min and max). If we want to try other values, we will need to increase the 2nd dimension of the two matrices
matrix_mean=array(NA,dim=c(nrow(free_param),2,nrow(mean_tot_original),2),dimnames=list(rownames(free_param),c("min","max"),rownames(mean_tot_original),c("modelI","modelII")))
matrix_amplitude=array(NA,dim=c(nrow(free_param),2,ncol(amplitude_tot_original),2),dimnames=list(rownames(free_param),c("min","max"),colnames(amplitude_tot_original),c("modelI","modelII")))

for(param_to_move in rownames(free_param)){
        id_param=grep(paste("^",param_to_move,sep=""),analyses) #This will help issue an error if we have more than 2 values per parameter
        for(i in 1:length(id_param)){
                mean_sens=log10(apply(N_sensitivity[,cpt,,id_param[i],],c(2,3),mean)) #mean_sens contains the log10 of average abundance per species, per model
                tmp_mean=100*(mean_tot_original-mean_sens)/mean_tot_original
                matrix_mean[param_to_move,i,,]=tmp_mean

                amp_sens=apply(log10(apply(N_sensitivity[,cpt,,id_param[i],],c(2,3),range)),2,diff) #amp_sens contains the log-ratio maximum/minimum per species, per model
                tmp_amplitude=100*(amplitude_tot_original-amp_sens)/amplitude_tot_original
                matrix_amplitude[param_to_move,i,,]=t(tmp_amplitude)

        }
}

pdf(paste("mean_abundance_amplitude_sensitivity_seed.pdf",sep=""),width=12,height=10)
par(mfrow=c(2,1),mar=c(3,5,3,1))
plot(0,0,t="n",xlim=c(0.75,nrow(free_param)+0.25),ylim=c(-22,28),xaxt="n",ylab="%Diff log10(mean abundance)",xlab="",cex.lab=1.65,cex.axis=1.65)
abline(h=0)
mtext("a",3,line=0.5,at=0.6,cex=1.75,font=2)
l=0
val_text=c()
at_val_text=c()
for(param_to_move in rownames(free_param)){
        l=l+1
        id_param=grep(paste("^",param_to_move,sep=""),analyses)
        space=0.5/(2*length(id_param))
        space_txt=0.8/(2*length(id_param))
        seq_space=seq(-space,space,length.out=length(id_param))
        seq_space_txt=seq(-space_txt,space_txt,length.out=length(id_param))
        for(i in 1:length(id_param)){

                for(mod in 1:2){
                per_change=mean(matrix_mean[param_to_move,i,,mod]*is.finite(matrix_mean[param_to_move,i,,mod]),na.rm=T)
                if(mod==1){
                        colm="lightgrey"
                        idplus=idplus_text=0
                }else{
                        colm="darkgrey"
                        idplus=0.08
                        idplus_text=0.
                }
                rect(l+0.75*seq_space[i]+idplus,min(c(0,per_change)),l+1.25*seq_space[i]+idplus,max(c(per_change,0)),col=colm)
                print(analyses[i])
                print(tab_summary[id_param[i],,mod])
                if(tab_summary[id_param[i],4,mod]==0){
                        if((tab_summary[id_param[i],"Persistence_coast",mod]<nspp)|(tab_summary[id_param[i],"Persistence_ocean",mod]<nspp)){
                        text(l+seq_space[i]+idplus,27,tab_summary[id_param[i],"Persistence_ocean",mod],col="red",cex=1.5)
                        }else{
                        points(l+seq_space[i]+idplus,27,t="p",pch=16,col="red",cex=1)
                        }
                }
                }
                at_val_text=c(at_val_text,l+seq_space_txt[i])
                tmp_text=strsplit(analyses[id_param[i]],"_")
                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
}
legend("bottomleft",c("Model I","Model II"),col=c("lightgrey","darkgrey"),pch=16,bty="n",cex=1.75)
par(mar=c(4.5,5,1.5,1))
plot(0,0,t="n",xlim=c(0.75,nrow(free_param)+0.25),ylim=c(-85,85),xaxt="n",ylab="%Diff log.amplitude",xlab="",cex.axis=1.65,cex.lab=1.65)
axis(1,labels=rownames(free_param),at=1:nrow(free_param),cex.axis=1.525)
abline(h=0)
mtext(c(all_others),1,line=3.5,at=1:nrow(free_param),cex=1.4)
mtext("b",3,line=0.5,at=0.6,cex=1.75,font=2)
l=0
val_text=c()
at_val_text=c()
for(param_to_move in rownames(free_param)){
        l=l+1
        id_param=grep(paste("^",param_to_move,sep=""),analyses)
        space=0.5/(2*length(id_param))
        seq_space=seq(-space,space,length.out=length(id_param))
        for(i in 1:length(id_param)){
                for(mod in 1:2){
                print(paste(analyses[id_param],mod))
                per_change=mean(matrix_amplitude[param_to_move,i,,mod]*is.finite(matrix_mean[param_to_move,i,,mod]),na.rm=T)
                print(per_change)
                if(mod==1){
                        colm="lightgrey"
                        idplus=0
                }else{
                        colm="darkgrey"
                        idplus=0.08
                }
                rect(l+0.75*seq_space[i]+idplus,min(c(0,per_change)),l+1.25*seq_space[i]+idplus,max(c(per_change,0)),col=colm)
		if(abs(per_change)>85){
			text(l+0.75*seq_space[i]+idplus,sign(per_change)*79+idplus*100,paste(format(per_change,digits=2),"%"))
		}
                }
                at_val_text=c(at_val_text,l+seq_space_txt[i])
                tmp_text=strsplit(analyses[id_param[i]],"_")
                val_text=c(val_text,tmp_text[[1]][length(tmp_text[[1]])])
        }
}
mtext(val_text,1,line=2.25,at=at_val_text,cex=1.3)
dev.off()
