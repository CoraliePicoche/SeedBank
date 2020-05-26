########
# CP 19/05/2020: variation of the seed bank access (mediated through mortality M) and the exchange rate with the ocean (mediated through e) to see the effect on richness and abundance
########

rm(list=ls())
graphics.off()
set.seed(42)


source("step_functions.r")


#######################Param used in every simulations
#Fixed parameters: golden parameter set
tab=read.table("simu.csv",sep=";",dec=".",header=T)
cyst_mortality=as.numeric(as.character(tab[tab[,1]=="cyst_mortality",2]))
cyst_burial=as.numeric(as.character(tab[tab[,1]=="cyst_burial",2]))
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
O_y=as.numeric(as.character(tab[tab[,1]=="overyielding",2]))
ratio_pos=as.numeric(as.character(tab[tab[,1]=="ratio_pos",2]))
coeff_Hij=as.numeric(as.character(tab[tab[,1]=="coeff_Hij",2]))

#Interaction
list_H_tmp=read.table("interaction_matrix_after_calibration.csv",sep=";",dec=".")
list_H=list(list_H_tmp,list_H_tmp)
type_inter_tmp=read.table("facilitation_and_competition.csv",sep=";",dec=".")
type_inter=list(type_inter_tmp,type_inter_tmp)
name_spp=rownames(list_H_tmp)
nspp=length(name_spp)

#Data to use (Auger)
a=as.matrix(read.table("../../param/reconstructed_temperature_Auger.txt", row.names=1,header=T,sep=";",dec="."))
temp_model=a[1:n_iter]

#Sinking rates and T_opt
tab=read.table(paste("../../param/species_specific_parameters.txt",sep=""),sep=";",dec=".",header=T)
tab=tab[name_spp,]
S=S_max*tab$S
T_opt=tab$T_opt+273.15+5
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


#param to change
M_bis=seq(0.1,1,length.out=10)
exch_bis=seq(0,1,length.out=10)

N_array=array(NA,dim=c(n_iter,3,nspp,length(M_bis),length(exch_bis)),dimnames=list(NULL,c("coast","ocean","seed"),name_spp,format(M_bis,digits=2),format(exch_bis,digits=2)))
N_array[1,,,,]=10^3

for(m_bis in 1:length(M_bis)){
	print(paste("M",M_bis[m_bis]))
	for(e_bis in 1:length(exch_bis)){
	        for(t in 1:(n_iter-1)){
                	var_tmp=step1(N_array[t,,,m_bis,e_bis],list_H,type_inter,temp_model[t],M_bis[m_bis],morta,a_d,T_opt,B)
               		Ntmp=var_tmp[[1]]
               		N_array[t+1,,,m_bis,e_bis]=step2(Ntmp,S,Gamma*(temp_model[t]>=temp_germin),exch_bis[e_bis])
		}
	}
}

##Look for richness
richness=matrix(NA,nrow=length(M_bis),ncol=length(exch_bis),dimnames=list(format(M_bis,digits=2),format(exch_bis,digits=2)))
average_abs=array(NA,dim=c(3,length(M_bis),length(exch_bis)),dimnames=list(c("coast","ocean","seed"),format(M_bis,digits=2),format(exch_bis,digits=2)))

pdf("abundance_f_e_M.pdf",width=10,height=10)
for(s in 1:length(name_spp)){
for(m_bis in 1:length(M_bis)){
	for(e_bis in 1:length(exch_bis)){
	#	richness[m_bis,e_bis]=sum(N_array[n_iter,"ocean",,m_bis,e_bis]>0)
		if(N_array[n_iter,'ocean',s,m_bis,e_bis]>0){
		#tmp=apply(N_array[id,'ocean',s,m_bis,e_bis],1,sum)
		tmp=N_array[id,'ocean',s,m_bis,e_bis]
		average_abs['ocean',m_bis,e_bis]=log10(mean(tmp))
		tmp=N_array[id,'coast',s,m_bis,e_bis]
		#tmp=apply(N_array[id,'coast',s,m_bis,e_bis],1,sum)
		average_abs['coast',m_bis,e_bis]=log10(mean(tmp))
		aplot=T
		}else{
		aplot=F
		}	
	}
}
	if(aplot){
#	heatmap(average_abs['ocean',,],main=name_spp[s])
    layout(matrix(c(1,2),nrow=1,ncol=2),widths=c(10,1),heights=1)
        par(mar=c(2, 2, 2, 2),oma=c(2,3,0,0))

        image(x=M_bis,y=exch_bis,average_abs["ocean",,],axes=FALSE,xlab="Morta",ylab="Exchange",col=heat.colors(10))
	title(name_spp[s])
        par(mar=c(2,.5,2,3))
	        seqi=seq(min(c(average_abs["ocean",,]),na.rm=T),max(c(average_abs["ocean",,]),na.rm=T),length.out=100)

        image(x=1,y=seqi,matrix(seqi,nrow=1,ncol=100),col=heat.colors(100),axes=FALSE,xlab="",ylab="")
axis(4,cex.axis=0.5,las=1)



	}
}
dev.off()

			


