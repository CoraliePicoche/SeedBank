###10/03/2020 CP: just drawing random, but definitive values for sinking rates of all species

rm(list=ls())
graphics.off()

pop_table=read.table("param/abundances_Auger.txt")
name_spp=rownames(pop_table)[!is.na(pop_table[,1])]
name_spp=name_spp[-1]
nspp=length(name_spp)

S=rbeta(nspp,0.55,1.25)*30/100
names(S)=name_spp
#Manip so that Chaetoceros spp have the highest sinking rate
tmp_S=S
or=order(S,decreasing=T)
tmp_S[c("CHA","THP")]=S[c(or[1],or[2])]
tmp_S[or[1]]=S[1]
tmp_S[or[2]]=S[2]
S=tmp_S
names(S)=name_spp

write.table(S,"param/sinking_rates_simulated.txt",sep=";",dec=".")
