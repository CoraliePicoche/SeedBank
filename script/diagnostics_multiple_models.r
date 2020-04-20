#####
#CP 08/04/2020  Compares the output of different models
####

rm(list=ls())
graphics.off()

m_list=c("v1.1","v2.1")
model_list=paste("model",m_list,sep="")
name_model=c("Classical BH","Saturating interactions")

collist=c("orchid","cyan","red")

#Just looking for the first result to be able to initialize the table
#Can also be seen as a check that dimensions match between models
tab_coast_init=read.table(paste("output/",model_list[1],"/out_coast.csv",sep=""),sep=";",dec=".")
name_spp=colnames(tab_coast_init)
n_iter=nrow(tab_coast_init)

important_id=365 #takes the last 365 days of the simulation
id=(n_iter-important_id+1):n_iter
seq_id=1:important_id

transfo_N=array(NA,dim=c(length(model_list),3,important_id,length(name_spp)),dimnames=list(name_model,c("coast","ocean","seed"),1:important_id,name_spp))

#Register data
for (m in 1:length(model_list)){
	tab_coast=as.matrix(read.table(paste("output/",model_list[m],"/out_coast.csv",sep=""),sep=";",dec="."))
	tab_ocean=as.matrix(read.table(paste("output/",model_list[m],"/out_ocean.csv",sep=""),sep=";",dec="."))
	tab_seed=as.matrix(read.table(paste("output/",model_list[m],"/out_seed.csv",sep=""),sep=";",dec="."))
	
	transfo_N[m,"coast",,name_spp]=log10(tab_coast[id,name_spp]+10^(-5))
	transfo_N[m,"ocean",,name_spp]=log10(tab_ocean[id,name_spp]+10^(-5))
	transfo_N[m,"seed",,name_spp]=log10(tab_seed[id,name_spp]+10^(-5))
}

#Plot data

####Two options are possible for mean values
#If we use raw values of corres_hernandez, we avoid the artefacts created by the interpolation and the random value used when gaps are over 2 points in the time series, but we increase the mean value artificially as cells are counted only when they are numerous. The inverse is true when using interpolated data. This is only a matter of choice.
#abundances_tab=read.table(paste("param/","corres_hernandez_Auger.txt",sep=""),sep=";",header=T)
#dates=as.Date(abundances_tab$Date)
#abundances_tab=abundances_tab[year(dates)>=1996,name_spp]#Using data from 1996
#dates=dates[year(dates)>=1996]
#x_obs=apply(abundances_tab,mean,na.rm=T) 

pop_table=read.table("param/mean_interpolated_abundances_Auger.txt",sep=",",dec=".",header=T)
rownames(pop_table)=pop_table$sp
x_obs=pop_table[name_spp,"Mean_abundances"]
names(x_obs)=name_spp
tab_mean=read.table("param/mean_monthly_abundance.txt",sep=";",dec=".",header=T)


pdf(paste("output/comparison_model/timeseries_one_by_one",paste(m_list,collapse="-"),".pdf",sep=""),width=10)
par(mfrow=c(1,1))

id_observed=seq(13,12*30.5,length.out=12)

	for(i in 1:length(name_spp)){
        	plot(1:important_id,transfo_N[1,"coast",,name_spp[i]],main=name_spp[i],col=collist[1],t="o",pch=16,lty=1,xlab="time",ylab="abundance",ylim=c(range(c(transfo_N[,"coast",,name_spp[i]]))),xaxt="n")
	        axis(1,at=seq(seq_id[1],seq_id[length(seq_id)],by=30),labels=seq(1,365,by=30))
        	abline(h=log10(x_obs[name_spp[i]]))
	        points(id_observed,log10(tab_mean[,name_spp[i]]+10^(-5)),pch=17,col="red",cex=2)
		for(m in 2:length(model_list)){
			lines(1:important_id,transfo_N[m,"coast",,name_spp[i]],col=collist[m],t="o",pch=16,lty=1)
		}
                legend("bottomright",name_model,col=collist[1:length(model_list)],pch=c(16),bty="n")
	}
dev.off()

