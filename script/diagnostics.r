######
##06/01/2020: CP. Just graphs and diag for outputs
#####


colo=c(rep(c("red","orange","green","blue"),2),"red","orange")
apch=c(rep(16,4),rep(17,4),rep(18,2))
alty=c(rep(1,4),rep(2,4),rep(3,2))

N_coast=read.table("output/out_coast.csv",sep=";",dec=".")
N_ocean=read.table("output/out_ocean.csv",sep=";",dec=".")
N_seed=read.table("output/out_seed.csv",sep=";",dec=".")

sp=colnames(N_coast)
n_iter=nrow(N_coast)

transfo_N_coast=log10(N_coast+10^(-5))
transfo_N_ocean=log10(N_ocean+10^(-5))
transfo_N_seed=log10(N_seed+10^(-5))

id=n_iter:(n_iter-365)

pdf("exploratory/fig/all_in_out.pdf",width=16,height=16)
par(mfrow=c(3,1))

plot(id,transfo_N_coast[id,1],col=colo[1],t="o",pch=apch[1],ylim=range(transfo_N_coast[id,]),xaxt="n",ylab="Coast",xlab="",lty=alty[1])
for(i in 2:10){
points(id,transfo_N_coast[id,i],col=colo[i],t="o",pch=apch[i],lty=alty[i])
}

plot(id,transfo_N_ocean[id,1],col=colo[1],t="o",pch=apch[1],ylim=range(transfo_N_ocean[id,]),xaxt="n",ylab="Ocean",xlab="",lty=alty[1])
for(i in 2:10){
points(id,transfo_N_ocean[id,i],col=colo[i],pch=apch[i],t="o",lty=alty[i])
}

plot(id,transfo_N_seed[id,1],col=colo[1],t="p",pch=apch[1],ylim=range(transfo_N_seed[id,]),ylab="Seed",xlab="time",lty=alty[1])
for(i in 2:10){
points(id,transfo_N_seed[id,i],col=colo[i],pch=apch[i],t="o",lty=alty[i])
}

legend("bottomright",sp,col=colo,pch=apch,lty=alty)

dev.off()

pdf("exploratory/fig/one_by_one.pdf",width=10)
par(mfrow=c(1,1))

for(i in 1:length(sp)){
	aylim=range(c(transfo_N_coast[id,i]),c(transfo_N_ocean[id,i]),c(transfo_N_seed[id,i]))
	plot(id,transfo_N_coast[id,i],main=sp[i],col="lightblue",t="o",pch=16,lty=1,xlab="time",ylab="abundance",ylim=aylim)
	lines(id,transfo_N_ocean[id,i],col="darkblue",t="o",pch=16,lty=1)
	lines(id,transfo_N_seed[id,i],col="brown",t="o",pch=16,lty=1)
}
dev.off()

print(max(transfo_N_coast[id,"CHD"]))
print(mean(transfo_N_coast[id,"CHD"]))
print(mean(transfo_N_coast[id,"CHS"]))

pdf("exploratory/fig/growth_rate.pdf",width=10)
par(mfrow=c(4,3))
id_hot=(temp>293)[id[1:(length(id)-1)]]
for(i in 1:10){
	growth_rate_coast=diff(log(N_coast[id,i]))
	plot(log(N_coast[id[1:(length(id)-1)],i]),growth_rate_coast)
	points(log(N_coast[id[1:(length(id)-1)],i][id_hot]),growth_rate_coast[id_hot],pch=16,col="red")
	points(log(N_coast[id[1:(length(id)-1)],i][!id_hot]),growth_rate_coast[!id_hot],pch=16,col="blue")
}
dev.off()
