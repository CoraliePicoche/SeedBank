###10/03/2020 CP: just drawing random, but definitive values for sinking rates of all species, and niche width
###20/04/20 CP: Also computing a temperature time series once and for all so that we remove any source of variability other than the parameters and formulation of the different model
###22/04/CP: Correction to harmonize files used in different models

rm(list=ls())
graphics.off()
source("script/matrix_MAR_clean.r")
set.seed(42)

pop_table=read.table("param/mean_interpolated_abundances_Auger.txt",sep=",",header=T)
name_spp=pop_table$sp
nspp=length(name_spp)

##Sinking rate
S=rbeta(nspp,0.55,1.25)
names(S)=name_spp
#Manip so that Chaetoceros and Thalassa spp have the highest sinking rate
tmp_S=S
or=order(S,decreasing=T)
tmp_S[c("CHA","THP")]=S[c(or[1],or[2])]
tmp_S[or[1]]=S[1]
tmp_S[or[2]]=S[2]
S=tmp_S
names(S)=name_spp


##Phenology
table_phenology=read.table("param/generalist_specialist_spp_added_amplitude_season.csv",sep=";",dec=".",header=T)

tab_G=subset(table_phenology,Type=="G")
generalist_range=runif(nrow(tab_G),180,1500)
g_order=order(generalist_range)
length_order=order(tab_G$Mean_length)
Val_b=rep(NA,nrow(tab_G))
tab_bisG=cbind(tab_G,Val_b)
tab_bisG[length_order,6]=generalist_range[g_order]

tab_S=subset(table_phenology,Type=="S")
specialist_range=runif(nrow(tab_S),7,55)
s_order=order(specialist_range)
length_order=order(tab_S$Mean_length)
Val_b=rep(NA,nrow(tab_S))
tab_bisS=cbind(tab_S,Val_b)
tab_bisS[length_order,6]=specialist_range[s_order]

tmp=rbind(tab_bisG,tab_bisS)
tmp=tmp[as.character(name_spp),]
tmp=cbind(tmp,S)

write.table(tmp,"param/species_specific_parameters.txt",sep=";",dec=".")

#Data to use (Auger)
n_iter=10000
evt_tab=read.table(paste("param/Augerhydro.txt",sep=""),sep=";",header=T)
temp=evt_tab[,"TEMP"]
#We need to build a simulation for temperatures as we don't have enough data in the real dataset
min_temp=min(temp,na.rm=T)
max_temp=max(temp,na.rm=T)
sd_temp=sd(temp,na.rm=T)
mean_temp=mean(temp,na.rm=T)
theta=1.3
temp_model=273.15+mean_temp+theta*sd_temp*sin(2*pi*1:n_iter/365)+rnorm(n_iter,0,sd_temp*sqrt(1-theta^2/2))

write.table(temp_model,"param/reconstructed_temperature_Auger.txt",sep=";",dec=".")

#Interaction matrix: instead of using the .RData everytime, we write an interaction matrix once and for all
load(paste("param/Auger_pencen_null_regular_common_MO.RData",sep=""))
name_spp=colnames(cis$call$model$B)
B_matrix=clean_matrix(cis,signif=F)
rownames(B_matrix)=colnames(B_matrix)=name_spp
nspp=length(name_spp)
write.table(B_matrix,"param/community_matrix_B_Auger.csv",dec=".",sep=";")

