########################
#20/02/2020 CP: "Sensitivity" analysis for certain parameters
########################

rm(list=ls())
graphics.off()
library(RColorBrewer)

simu=1:9 #Germination and resuspension are the ones we want to look at
diag="CHA"

list_germination=rep(NA,length(simu))
list_resuspension=rep(NA,length(simu))

tab_coast=list()
tab_ocean=list()
tab_seed=list()

for(n_simulation in simu){
	tab=read.table(paste("param/simu",n_simulation,".csv",sep=""),sep=";",dec=".",header=T)
	list_germination[n_simulation]=as.numeric(as.character(tab[tab[,1]=="germination",2]))
	list_resuspension[n_simulation]=as.numeric(as.character(tab[tab[,1]=="resuspension",2]))

	tab_coast[[n_simulation]]=read.table(paste("output/out_coast",n_simulation,".csv",sep=""),sep=";",dec=".")
	tab_ocean[[n_simulation]]=read.table(paste("output/out_ocean",n_simulation,".csv",sep=""),sep=";",dec=".")
	tab_seed[[n_simulation]]=read.table(paste("output/out_seed",n_simulation,".csv",sep=""),sep=";",dec=".")
}

name_spp=colnames(tab_coast[[1]])

n_iter=nrow(tab_coast[[1]])

id=(n_iter-365):n_iter

pdf("output/comparison_abundances_CHA_germ_res.pdf")
par(mfrow=c(3,3))

for(j in 1:length(simu)){
	transfo_N_coast=log10(tab_coast[[j]][id,]+10^(-5))
#transfo_N_ocean=log10(tab_ocean[id,]+10^(-5))
	transfo_N_seed=log10(tab_seed[[j]][id,]+10^(-5))

	plot(id,transfo_N_coast[,diag],col="lightblue",t="o",main=paste("Germ",list_germination[j],"Res",list_resuspension[j]),pch=16,ylim=c(4,5))
#	points(id,transfo_N_seed[,diag],col="red",t="o",pch=16)
}
dev.off()

pdf("output/comparison_abundances_one_panel_coast.pdf")
acol=brewer.pal(9,"Blues")
leg=c()
par(mfrow=c(1,1))
plot(id,id,t="o",ylim=c(4,5),xlab="",ylab="")
for(j in 1:length(simu)){
	leg=c(leg,paste("Germ",list_germination[j],"Res",list_resuspension[j]))
       transfo_N_coast=log10(tab_coast[[j]][id,]+10^(-5))

       points(id,transfo_N_coast[,diag],col=acol[j],t="o",pch=16,cex=1)
}
legend("topleft",leg,pch=16,col=acol)
dev.off()

pdf("output/comparison_abundances_one_panel_seed.pdf")
acol=brewer.pal(9,"Oranges")
par(mfrow=c(1,1))
plot(id,id,t="o",ylim=c(4,7),xlab="",ylab="")
for(j in 1:length(simu)){
       transfo_N_seed=log10(tab_seed[[j]][id,]+10^(-5))

       points(id,transfo_N_seed[,diag],col=acol[j],t="o",pch=16,cex=1)
}
legend("bottomleft",leg,pch=16,col=acol,ncol=2)

dev.off()
