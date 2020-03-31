############
## CP 18/01/2020: Cleaner version of quadratic programming and growth rate inference
############

library(limSolve)
library(Matrix)
source("step_functions.r")
#source("../../script/step_functions.r")

#compute b parameter in Scranton & Vasseur model
f_to_optimize_B=function(b,T_min,T_max,T_opt,A){
        f1=integrate(Vectorize(growth_rate_SV),lower=T_min-5,upper=T_max+5,T_opt,b)$val
        tmp=abs(f1-A)
        return(tmp)
}

#Turn the community matrix B from a MAR model to the interaction matrix A in a Beverton-Holt model. Based on matrix B and abundances at equilibrium N
MAR2BH=function(B,N){
        A=matrix(NA,nrow=nrow(B),ncol=ncol(B))
	rownames(A)=colnames(A)=rownames(B)
        for(i in 1:dim(B)[1]){
                sumBJ=sum(B[i,]*N)
                sumB=sum(B[i,])
                for(j in 1:dim(B)[2]){
                        A[i,j]=-1/N[j]*B[i,j]/(1+sumB)
                }
        }

        return(A)
}

#Uses quadratic programming (Maynard et al. 2019) to improve the estimates for the interaction matrix A, and, optionnally, for the growth rate r. It requires the abundances at equilibrium N. tol is the tolerance for the relative change in parameters (10% by default)
quadprog=function(A,N,r_mean,r_calibrate=F,tol=0.1){
	nspp=nrow(A)
	name_spp=colnames(A)
	Avec=as.numeric(t(A))
	Alow <- Avec-abs(Avec*tol)
	Aupp <- Avec+abs(Avec*tol)
        # constrain the diagonals to be positive (density-dependance
        Alow[as.logical(as.numeric(diag(nspp)))] = rep(0,nspp)

        #Also keep the sign from A: species that compete should keep a positive coefficient. Same for species that are mutualist (negative coefficients).
        Alow[Avec>0]=0
        Aupp[Avec<0]=0

	parvec <- c(Avec)
	r_2=1-exp(r_mean)
	parvec <- c(parvec,r_2)

	if(r_calibrate){
		rupp=rep(0,length(r_2)) #Negative value only for r2 so that actual r>1 -> species have a positive maximum growth rate
		rlow <- r_2-abs(r_2*tol) #We have a fixed value of r
	}else{
		no_tol=10e-2
		rlow <- r_2-abs(r_2*no_tol) #We have a fixed value of r
		rupp <- r_2+abs(r_2*no_tol) #We have a fixed value of r
	
	}
	h <- c(Alow,rlow,-Aupp,-rupp)
	npar <- length(parvec)

	# inequality constraints
	G <- rbind(diag(npar),-diag(npar))
	# equality constraints, assuming at equilibirum
	E <- cbind(as.matrix(bdiag(replicate(nspp,matrix(N,nrow=1),simplify = F))),diag(nspp))
	f <- rep(0,nrow(E))
	# fit the qp model, returning a silent warning if it doesn't converge
	fit=limSolve::lsei(A = diag(npar), B = parvec, E = E, F = f, G = G, H = h,type=1,verbose=TRUE,fulloutput=TRUE)
	
	A_after=t(matrix(fit$X[1:nspp^2],nspp,nspp))
	colnames(A_after)=rownames(A_after)=name_spp

	r_after=log(1-fit$X[(nspp^2+1):length(fit$X)])
	
	return(list(A_after,r_after))
}
	
