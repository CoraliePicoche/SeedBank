rm(list=ls())
graphics.off()

#Pop would be array [date,sp]
BC_index=function(pop_A,pop_B){
	id_den=0
	id_num=0
	sp=colnames(pop_A)
	for (s in sp){
		mean_A=mean(pop_A[,sp])	
		mean_B=mean(pop_B[,sp])
		id_num=id_num+min(mean_A,mean_B)
		id_den=id_den+mean_A+mean_B
	}
	id_bc=1-(2*id_num)/id_den
	return(id_bc)
}
