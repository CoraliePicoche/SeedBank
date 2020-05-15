#####
#14/05/2020 CP: Illustrate the thermal preference of each species
#####

rm(list=ls())
graphics.off()

growth_rate_noMTE_Bissinger=function(temp,T_opt,B,d){ 
        ftmp=rep(NA,length(T_opt))
        rtmp=rep(NA,length(T_opt))
        metabolism=0.81*exp(0.0631*(temp-273))*d
        for(i in 1:length(T_opt)){
                if(temp<=T_opt[i]){
                        ftmp[i]=exp(-(abs(temp-T_opt[i]))^3/B[i])
                }else{
                        ftmp[i]=exp(-5*(abs(temp-T_opt[i]))^3/B[i])
                }
                rtmp[i]=ftmp[i]*metabolism
        }
        return(rtmp)
}

temp=seq(0,30,length.out=1000)+273

tab=read.table(paste("param/species_specific_parameters.txt",sep=""),sep=";",dec=".",header=T)
T_opt=tab$T_opt+273+5
B=tab$Val_b
sp=rownames(tab)

mat_val=matrix(NA,ncol=length(sp),nrow=length(temp))

for(t in 1:length(temp)){
	mat_val[t,]=growth_rate_noMTE_Bissinger(temp[t],T_opt,B,0.5)
}


aylim=c(0,max(c(mat_val)))

pdf("output/thermal_niche.pdf",width=10,height=8)
par(mfrow=c(3,4),mar=c(4,4.5,2,1))
for(s in 1:length(sp)){
	if(s%%4==1){
		yl=expression(paste("Growth rate (d"^"-1",")",sep=""))
	}else{
		yl=""
	}
	if(s>8){
		xl="Temperature"
	}else{
		xl=""
	}
	plot(temp-273,mat_val[,s],t="l",main=sp[s],xlab=xl,ylab=yl,ylim=aylim)
}
dev.off()


