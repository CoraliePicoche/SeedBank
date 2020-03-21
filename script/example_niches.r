#CP 10/03/2020: Just a small script to give a conceptual examples of the different niches that can be seen
#21/03/2020: Just changed the limits to have a nicer plot

rm(list=ls())
graphics.off()
source("script/infer_interaction_matrix_growth_rate.r")

T_min=0+273
T_max=30+273
T_opt=seq(T_min+5,T_max-5,length.out=4)
T_seq=seq(T_min-5,T_max+5,length.out=1000)
nspp=length(T_opt)

#SV model, keeping the same area
A=15
B=rep(NA,nspp)
for(i in 1:nspp){
	B[i]=optimize(f_to_optimize_B,T_min,T_max,T_opt[i],A,interval=c(0,100000))$minimum
}

mat_sv=matrix(NA,length(T_seq),nspp)
for(i in 1:length(T_seq)){
	mat_sv[i,]=growth_rate_SV(T_seq[i],T_opt,B)
}

#Generalist vs specialist
poly=function(T_opt,h,w,x){
		a=h
		T1=T_opt-w/2
		
		den=(T1-T_opt/2-T1*T1/(2*T_opt))
		b=-h/den

		a=h-b*T_opt/2
		c=-b/(T_opt*2)

		tmp=a+b*x+c*x*x
		return(tmp)
	}

#2 generalist
#other 
h=c(1.5,1.2,0.75,1.5)
w=c(10,45,10,45)

mat_genspe=matrix(NA,length(T_seq),nspp)
for(i in 1:nspp){
	mat_genspe[,i]=poly(T_opt[i],h[i],w[i],T_seq)
}

mat_regular=matrix(NA,length(T_seq),nspp)
for(i in 1:nspp){
	mat_regular[,i]=poly(T_opt[i],1.,15,T_seq)
}

pdf("output/examples_niches.pdf",width=12,height=5)
acol=c("blue","purple","orange","red")
par(mfrow=c(1,3))
plot(0,0,xlim=range(T_seq),ylim=range(mat_sv),t="o",xlab="Temperature",ylab="growth_rate",main="Fixed niche area")
for(i in 1:nspp){
	lines(T_seq,mat_sv[,i],col=acol[i],lty=1,lwd=2)
}
plot(0,0,xlim=range(T_seq),ylim=range(mat_sv),t="o",xlab="Temperature",ylab="",main="Generalist vs Specialist")
for(i in 1:nspp){
	lines(T_seq,mat_genspe[,i],col=acol[i],lwd=2)
}
plot(0,0,xlim=range(T_seq),ylim=range(mat_sv),t="o",xlab="Temperature",ylab="",main="Same growth and tolerance")
for(i in 1:nspp){
	lines(T_seq,mat_regular[,i],col=acol[i],lwd=2)
}
dev.off()
