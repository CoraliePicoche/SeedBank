rm(list=ls())
graphics.off()

####Model I
#Interaction
tmp_inter_I=read.table("../../output/modelv1.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
tmp_inter=as.matrix(tmp_inter_I)
diag(tmp_inter)=NA
####ModelII
list_H_tmp=read.table("../../output/modelv2.1/interaction_matrix_after_calibration.csv",sep=";",dec=".")
list_H=list(list_H_tmp,list_H_tmp)
type_inter_tmp=read.table("../../output/modelv2.1/facilitation_and_competition.csv",sep=";",dec=".")


list_H_tmp=as.matrix(type_inter_tmp/list_H_tmp) #We compute at low densities
diag(list_H_tmp)=NA
name_spp=rownames(tmp_inter_I)

quantif_inter=array(NA,dim=c(length(name_spp),10,2),dimnames=list(name_spp,c("AbsImpact","AbsVul","MeanImpact","MeanVul","PercPosImpact","PercPosVul","MeanImpactComp","MeanImpactFacil","MeanVulComp","MeanVulFacil"),c("ModelI","ModelII")))


for(s in 1:length(name_spp)){
	
	quantif_inter[s,"AbsImpact","ModelI"]=mean(abs(tmp_inter[,s]),na.rm=T)
	quantif_inter[s,"AbsImpact","ModelII"]=mean(abs(list_H_tmp[,s]),na.rm=T)
	quantif_inter[s,"AbsVul","ModelI"]=mean(abs(tmp_inter[s,]),na.rm=T)
	quantif_inter[s,"AbsVul","ModelII"]=mean(abs(list_H_tmp[s,]),na.rm=T)
	quantif_inter[s,"MeanImpact","ModelI"]=mean(tmp_inter[,s],na.rm=T)
	quantif_inter[s,"MeanImpact","ModelII"]=mean(list_H_tmp[,s],na.rm=T)
	quantif_inter[s,"MeanVul","ModelI"]=mean(tmp_inter[s,],na.rm=T)
	quantif_inter[s,"MeanVul","ModelII"]=mean(list_H_tmp[s,],na.rm=T)

	
	quantif_inter[s,"MeanImpactComp","ModelI"]=mean(tmp_inter[which(tmp_inter[,s]>0),s],na.rm=T)
	quantif_inter[s,"MeanImpactComp","ModelII"]=mean(list_H_tmp[which(list_H_tmp[,s]>0),s],na.rm=T)
	quantif_inter[s,"MeanImpactFacil","ModelI"]=mean(tmp_inter[which(tmp_inter[,s]<0),s],na.rm=T)
	quantif_inter[s,"MeanImpactFacil","ModelII"]=mean(list_H_tmp[which(list_H_tmp[,s]<0),s],na.rm=T)
	
	quantif_inter[s,"MeanVulComp","ModelI"]=mean(tmp_inter[s,which(tmp_inter[s,]>0)],na.rm=T)
	quantif_inter[s,"MeanVulComp","ModelII"]=mean(list_H_tmp[s,which(list_H_tmp[s,]>0)],na.rm=T)
	quantif_inter[s,"MeanVulFacil","ModelI"]=mean(tmp_inter[s,which(tmp_inter[s,]<0)],na.rm=T)
	quantif_inter[s,"MeanVulFacil","ModelII"]=mean(list_H_tmp[s,which(list_H_tmp[s,]<0)],na.rm=T)


	quantif_inter[s,"PercPosImpact","ModelI"]=sum(tmp_inter[,s]>0,na.rm=T)/(10)*100
	quantif_inter[s,"PercPosImpact","ModelII"]=sum(list_H_tmp[,s]>0,na.rm=T)/(10)*100
	quantif_inter[s,"PercPosVul","ModelI"]=sum(tmp_inter[s,]>0,na.rm=T)/(10)*100
	quantif_inter[s,"PercPosVul","ModelII"]=sum(list_H_tmp[s,]>0,na.rm=T)/(10)*100


}

pdf("inter_per_species.pdf")
par(mfrow=c(2,3))
for(model in c("ModelI","ModelII")){
	plot(0,0,t="n",xlab=name_spp,xlim=c(0,12),ylab=c("AllInter"),ylim=range(c(quantif_inter[,1:4,model]),na.rm=T))
	for(s in 1:length(name_spp)){
		if(name_spp[s] %in% c("DIT","GUI","LEP","SKE","PRO","PRP")){
			abline(v=s,lty=2,lwd=0.5)
		}
		points(s,quantif_inter[s,"AbsImpact",model],col="red",pch=15)
		points(s,quantif_inter[s,"AbsVul",model],col="blue",pch=15)
		points(s,quantif_inter[s,"MeanImpact",model],col="red",pch=16)
		points(s,quantif_inter[s,"MeanVul",model],col="blue",pch=16)
	}
	axis(1,at=1:11,labels=name_spp)
	legend("bottomleft",c("Impact","Vulnerability","Absolute","Mean"),col=c("red","blue","black","black"),pch=c(15,15,16,16),bty="n")
	
	plot(0,0,t="n",xlab=name_spp,xlim=c(0,12),ylab=c("% Pos"),ylim=range(c(quantif_inter[,5:6,model]),na.rm=T))
	for(s in 1:length(name_spp)){
		if(name_spp[s] %in% c("DIT","GUI","LEP","SKE","PRO","PRP")){
			abline(v=s,lty=2,lwd=0.5)
		}
		points(s,quantif_inter[s,"PercPosImpact",model],col="red",pch=16)
		points(s,quantif_inter[s,"PercPosVul",model],col="blue",pch=16)
	}
	axis(1,at=1:11,labels=name_spp)
	
	plot(0,0,t="n",xlab="",xlim=c(0,12),ylab=c("Inter per sign"),ylim=range(c(quantif_inter[,7:10,model]),na.rm=T))
	for(s in 1:length(name_spp)){
		if(name_spp[s] %in% c("DIT","GUI","LEP","SKE","PRO","PRP")){
			abline(v=s,lty=2,lwd=0.5)
		}
		points(s,quantif_inter[s,"MeanImpactComp",model],col="red",pch=17)
		points(s,quantif_inter[s,"MeanImpactFacil",model],col="red",pch=16)
		points(s,quantif_inter[s,"MeanVulComp",model],col="blue",pch=17)
		points(s,quantif_inter[s,"MeanVulFacil",model],col="blue",pch=16)
	}
	axis(1,at=1:11,labels=name_spp)
	legend("bottomleft",c("Impact","Vulnerability","Competitive","Facilitative"),col=c("red","blue","black","black"),pch=c(16,16,17,16),bty="n")
}
dev.off()
