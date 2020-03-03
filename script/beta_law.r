rm(list=ls())
graphics.off()

val=c(9.5,5.1,9.6,8.1,6.0,8.3,15.5,7.8,19.7,9.6,6.6,3.7)
val_m=mean(val)

dist_sinking=rbeta(10000,0.55,1.25)*30

#We have val_m prop to mean(dist_sinking)

pdf("output/beta_law.pdf")
hist(dist_sinking,freq=F,xlab="Sinking rate (%)",ylab="Freq",main="")
dev.off()
