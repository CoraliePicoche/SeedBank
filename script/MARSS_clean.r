#Analyse
analyse_MARSS=function(tab_sp,tab_cov,B1,filename,C1="unconstrained",amethod="kem",pi1="zero",V1="identity",aalpha=0.05,min_iter=15,max_iter=500,iter_boot=1000,aconv.test.slop.tol=0.001,aabstol=0.001,boot=TRUE,Q1="diagonal and unequal",aicb=FALSE){
#Matrix        
	U1="zero" #Offset for processes
        Z1=diag(1,dim(tab_sp)[1],dim(tab_sp)[1]) #Idenity between Y and X
        A1="zero" #Offset for observations
        R1="zero" #Covariance matrix for the observation matrix, zero is wrong but like Griffiths, seem that we're not able to do before
	c1=tab_cov

        if(amethod=="kem"){
		cntl.list=list(conv.test.slope.tol=aconv.test.slop.tol,minit=min_iter,maxit=max_iter,abstol=aabstol)
	}else{
		cntl.list=list()
	}
	model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1,C=C1,c=c1)

	print(Sys.time())
	print(paste("Start estimation for",filename,sep=" "))
        fit_log=MARSS(tab_sp, method=amethod,model=model.list,control=cntl.list)
	print("Stop estimation")

	if(boot){
	        print("Start bootstrap")
        	cis=MARSSparamCIs(fit_log,method="parametric",nboot=iter_boot,alpha=aalpha)
        	save(fit_log,cis,file=filename)
		}
	else{
        	save(fit_log,file=filename)
	}
	if(aicb){
	        print("Start AIC")
        	aic=MARSSaic(fit_log,output=c("AIC","AICc", "AICbp", "boot.params"),Options=list(nboot=iter_boot))
        	save(fit_log,aic,file=paste("aic_",filename,sep=""))

	}
        print(paste("Sauvegarde de",filename,sep=" "))
}


