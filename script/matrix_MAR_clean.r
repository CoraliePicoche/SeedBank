#2018/09/07 CP 
#I need a function that turns the fit_log$B format to a readily usable interaction matrix

clean_matrix=function(cis, diag_bool=TRUE,signif=FALSE){
	require('stringr')
                        sp=dimnames(cis$model$data)[[1]] #Studied species
                        var=dimnames(cis$call$model$c)[[1]] #Studied covariates
                        nom=dimnames(cis$par$B)[[1]] #Names of the interactions
                        #get the value from the MARSS object, according to the names of the coefficients
			if(diag_bool){
                        	B=diag(-1,length(sp),length(sp))
			}else{
                        	B=matrix(0,nrow=length(sp),ncol=length(sp))
			}
			for (n in 1:length(nom)){
                                if(grepl("^\\(",nom[n])){ #default when using unconstrained or diagonal, names are of the form (1,2) for (i,j)
                                        i=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][2])
                                        j=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][3])
                                        if(is.na(i)){ #i and j might not be numeric
                                                i=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][2]))))
                                                j=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][3]))))
                                        }
                                  }else{
                                        a=str_split(nom[n],sp) #The coefficient is XY where X eats Y, X and Y being in the species list
                                        for (abis in 1:length(a)){
                                                if(length(a[[abis]])==2){
                                                        if((a[[abis]][1]=="")&&(a[[abis]][2]!="")){
                                                                j=abis
                                                        }
                                                        if((a[[abis]][2]=="")&&(a[[abis]][1]!="")){
                                                                i=abis
                                                        }
                                                }else if(length(a[[abis]])==3){
                                                        if((a[[abis]][2]=="")&&(a[[abis]][1]=="")&&(a[[abis]][3]=="")){
                                                                j=abis
                                                                i=abis
                                                        }
                                                }
                                        }
                                }
				if(signif){
					if(cis$par.upCI$B[n]*cis$par.lowCI$B[n]>0){
	                                	B[i,j]=B[i,j]+cis$par$B[n]
					}
				}else{
	                                B[i,j]=B[i,j]+cis$par$B[n]
				}
                        }
	return(B)
}
