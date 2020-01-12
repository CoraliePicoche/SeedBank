graphics.off()
rm(list=ls())

list_lookalike=c("Chaetoceros diadema","Chaetoceros socialis","Thalassiosira eccentrica","Skeletonema","Gyrosigma attenuatum","Pseudo-nitzschia","Nitzschia","Navicula transitans","Gymnodinium","Protoperidinium")
#For Gyrosigma and Navicula, I am specific because I want a certain biovolume

f_size=read.csv("PEG_BVOL2016.csv",sep=";",header=TRUE)
biovolume=c()
sur=c()
for (sp in list_lookalike){
	id=grep(sp,f_size$Species)
	bio=as.numeric(as.character(f_size$Calculated_volume_µm3[id]))
	l1=as.numeric(as.character(f_size$Length.l1.µm[id]))
	l2=as.numeric(as.character(f_size$Length.l2.µm[id]))
	
	d1=as.numeric(as.character(f_size$Diameter.d1.µm[id]))
	d2=as.numeric(as.character(f_size$Diameter.d2.µm[id]))

	surface1=l1*l2
	surface2=pi*(d1/2)^2
	surface3=pi*(d2/2)^2

	if((length(surface1)>0)){
		if(!is.na(surface1)){surface=surface1}
	}
	if((length(surface2)>0)){
		if(!is.na(surface2)){surface=surface2}
	}
	if((length(surface3)>0)){
		if(!is.na(surface3)){surface=surface3}
	}

	if(is.na(median(bio,na.rm=T))){stop()}

	biovolume=c(biovolume,median(bio,na.rm=T))
	sur=c(sur,median(surface,na.rm=T))
}

data=cbind(list_lookalike,biovolume,sur)

write.csv(data,"param_sinking.csv")
