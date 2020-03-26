###10/03/2020 CP: just drawing random, but definitive values for sinking rates of all species

rm(list=ls())
graphics.off()
set.seed(42)

pop_table=read.table("param/abundances_Auger.txt",sep=",",header=T)
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
table_phenology=read.table("param/generalist_specialist_spp.csv",sep=";",dec=".",header=T)

tab_G=subset(table_phenology,Type=="G")
generalist_range=runif(nrow(tab_G),180,1500)
g_order=order(generalist_range)
length_order=order(tab_G$Mean_length)
Val_b=rep(NA,nrow(tab_G))
tab_bisG=cbind(tab_G,Val_b)
tab_bisG[length_order,4]=generalist_range[g_order]

tab_S=subset(table_phenology,Type=="S")
specialist_range=runif(nrow(tab_S),7,55)
s_order=order(specialist_range)
length_order=order(tab_S$Mean_length)
Val_b=rep(NA,nrow(tab_S))
tab_bisS=cbind(tab_S,Val_b)
tab_bisS[length_order,4]=specialist_range[s_order]

tmp=rbind(tab_bisG,tab_bisS)
tmp=tmp[as.character(name_spp),]
tmp=cbind(tmp,S)

write.table(tmp,"param/species_specific_parameters_simulated.txt",sep=";",dec=".")
