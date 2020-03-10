#CP 10/03/2020: Fake figures to have an idea of what we're aiming at

rm(list=ls())
graphics.off()

set.seed(42)

high_compet_with_seed=rnorm(100,10,2)
low_compet_with_seed=rnorm(100,11,1)

high_compet_without_seed=rnorm(100,6,2) 
low_compet_without_seed=rnorm(100,8,1) 

high_temp_with_seed=runif(100,4,6)
low_temp_with_seed=runif(100,5,7)

high_temp_without_seed=runif(100,2,5)
low_temp_without_seed=runif(100,3,5)

pdf("output/expected_change_in_tolerance.pf",width=10,height=6)
par(mfrow=c(1,2),mar=c(4,4,2,2))
boxplot(list(high_compet_with_seed,low_compet_with_seed,high_compet_without_seed,low_compet_without_seed),at=c(1,2,4,5),ylab="Final richness",xaxt="n",col=c("red","blue"))
axis(1,at=c(1.5,4.5),labels=c("W/ Seed","W/o Seed"))
legend("topright",c("High compet","Low compet"),fill=c("red","blue"),bty="n")

boxplot(list(high_temp_with_seed,low_temp_with_seed,high_temp_without_seed,low_temp_without_seed),at=c(1,2,4,5),ylab="log10(Primary production)",xaxt="n",col=c("red","blue"))
axis(1,at=c(1.5,4.5),labels=c("W/ Seed","W/o Seed"))
legend("bottomleft",c("High temp","Low temp"),fill=c("red","blue"),bty="n")
dev.off()

seed_morta=c(0.1,0.5,0.9)
exchange=c(0.1,0.5,0.9)

possible_richness=matrix(c(5,7,8,9,10,11,10,11,11),3,3,byrow=T)

rbPal <- colorRampPalette(c('blue','red'))

#This adds a column of color values
# based on the y values
val_col <- (matrix(rbPal(9)[as.numeric(cut(c(possible_richness),breaks = 9))],3,3,byrow=F))

pdf("output/expected_change_in_richness.pf")
plot(0,0,xlim=c(0,1),ylim=c(0,1),ylab="seed_morta",xlab="exchange",t="n")
for(i in 1:3){
	for(j in 1:3){
		points(exchange[j],seed_morta[i],col=val_col[j,4-i],pch=16,cex=3)
	}
}
dev.off()
