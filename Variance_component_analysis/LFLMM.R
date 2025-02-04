LFLMM <-
function(Y, X, a=rep(0,nrow(Y)), beta=NULL, psi=NULL, delta=NULL, Psi=NULL, theta=NULL, omega=NULL, rho=1, subset=NULL, nLatent=0, forced=F, ITRMAX=30, PLOT=T, nexpcells=NULL, xgamma, batch=rep(1,ncol(Y)), s0=NULL){

	if(nrow(X)!=ncol(Y)){warning("Input TPM matrix Y is not compatible with covariate matrix X."); return()}
	if(!is.data.frame(X)){warning("Covariate matrix X is not a data frame."); return()}
	
	if(is.null(subset)){subset=rep(T, nrow(X))}
	nh = 1
	Z = rep(1,sum(subset))
	isnum = 0
	for(i in seq(ncol(X))){
print(class(X[,i]))
		if(is.character(X[,i]) || is.factor(X[,i])){
			isnum = c(isnum, 0)
			Z1 = array(0,c(sum(subset),length(table(X[subset,i]))))
			Z1[!is.na(X[subset,i]),] = model.matrix(~0+X[subset,i])
		}else{
			isnum = c(isnum, 1)
			Z1 = matrix(scale(as.numeric(X[subset,i])),sum(subset))
			Z1[is.na(Z1)]=0
		}
		Z  = cbind(Z, Z1)
		nh = c(nh, ncol(Z1))
	}
	names(nh)=c("Intercept", names(X))
	
	isnum = cumsum(1-isnum)*isnum
	print(isnum)
	Collapse=function(x,flag){
		x1=x[cumsum(flag)==0]
		x2=x[rev(cumsum(rev(flag)))==0]
		res = sum(x[flag])
		names(res) = paste(names(x[flag]),collapse="/")
		return(c(x1, res, x2))
	}
	if(sum(isnum)>0){for(i in seq(max(isnum))){
		if(sum(isnum==i)>2){
			if(!forced){
				mflag = readline(paste("Do you want to merge ", paste(names(nh)[isnum==i],collapse=","), " factors? [N/y]: ",sep=""))
			}else{
				mflag="y"
			}
			if(sum(mflag%in%c("y","Y"))){
				nh    = Collapse(nh, isnum==i)
				isnum = Collapse(isnum, isnum==i)
			}
		}
	}}
	if(nLatent>0){nh=c(nh, LatentFactor=nLatent);}
print(nh)	
	K = length(nh)
	H = diag(K)[rep(seq(K),nh),]
	X = Z

	# Latent factors
	#Z = cbind(Z, array(0,c(nrow(Z),nLatent)))
	if(!is.null(Psi)){
		Z = cbind(Z, Psi)
	}else{
		Z = cbind(Z, array(0,c(nrow(Z),nLatent)))
	}

	# Likelihood
	logDet = function(x){eval = eigen(x,T)[[1]]; sum(log(eval))}
	Kurt=function(V){diag(V)%*%t(diag(V))+2*V^2}
	Solve = function(X){if(length(X)==1){c(1/X)}else{X=eigen(X,T); r=X[[1]]; r[r<1e-7]=1; X[[2]]%*%diag(1/r)%*%t(X[[2]])}}
	lkhd = function(beta, psi, delta, omega, theta, s0, rho, batch, Y, X, Z, H, Yt, YYt, YX, YZ){
		J = nrow(Y)
		N = ncol(Y)
		M = ncol(Z)
		L = length(unique(batch))
		z = c(apply(Z, 2, sum))
		R = rho[batch]
		ZtZ  = t(Z)%*%Z # for W
		
		# A = (ZtZ+U^-2) MAT
		u = c(H%*%psi)
		A = ZtZ+diag(1/u^2)

		print("sigma2")
		# signa^2 | beta delta Psi(Z) psi(A)
		ZtRYt   = t(as.matrix(Y%*%(Z*R)))
		ZtRX    = t(Z*R)%*%X
		ZtRYtmb = (ZtRYt - c(ZtRX%*%beta))
		AinvZtRYtmb = solve(A, ZtRYtmb)
		yR2y = colSums((Yt*R)^2)
		Xb   = c(X%*%beta)
		YR2Xb = as.matrix(Y%*%(Xb*R^2))
		if(is.null(s0)){
			s2 = c( ( yR2y - 2*YR2Xb + sum((Xb*R)^2) ) - colSums(AinvZtRYtmb*ZtRYtmb) )
			if(!is.null(nexpcells)){
				tau    = (N/2+theta) / (s2/2+theta*c(nexpcells%*%omega))
				logtau = log(tau/(N/2+theta)) + digamma((N/2+theta))
				resglm = GammaGlm(tau[xgamma>0.05],logtau[xgamma>0.05],nexpcells[xgamma>0.05,,drop=F],PLOT=1,x=xgamma[xgamma>0.05]); omega=resglm$omega; theta=resglm$theta
				cat("omega theta=");print(c(omega,theta))
				s2 = (s2/2+theta*c(nexpcells%*%omega)) / (N/2+theta)
			}else{
				s2 = (s2+1)/(N+1)
			}
		}else{s2=s0}

		print("beta")
		# beta | delta sigma2 Psi(Z) psi(A)
		XtR2X = t(X*R)%*%(X*R)
		y = colSums(Y/s2) / sum(1/s2)
		beta0=0.01
		beta = c(Solve(diag(ncol(X))*beta0 + XtR2X-t(ZtRX)%*%solve(A)%*%ZtRX) %*% (t(X*R^2)%*%y-t(ZtRX)%*%solve(A,t(Z*R)%*%y)))
		Xb   = c(X%*%beta)

		print("Psi")
		# Psi | beta delta psi(A) sigma2
		PL=NULL
		if(nLatent>0){
			cat("Latent factor")
			M2 = M-nLatent
			A2 = t(Z[,1:M2])%*%Z[,1:M2] + diag(1/u[1:M2]^2)
			if(length(rho)>1){
				RXb   = c((X*R)%*%beta)
				YR2Xb = c(as.matrix(Y%*%c(R*RXb)))
				YRZ   = as.matrix(Y%*%(Z[,1:M2]*R))
				ZtRXb = c(t(Z[,1:M2])%*%RXb)
				CAind = c(YRZ%*%solve(A2, ZtRXb))
				YVYt = as.matrix(Y%*%(Yt*R^2))
			}else{
				RXb   = c(X%*%beta)
				YR2Xb = c(YX%*%beta)
				YRZ   = YZ[,1:M2]
				ZtRXb = c(t(Z[,1:M2])%*%RXb)
				CAind = c(YRZ%*%solve(A2, ZtRXb))
				YVYt  = YYt
			}

			YVYt = t(YVYt-YR2Xb)-YR2Xb + sum(RXb^2) - YRZ%*%solve(A2)%*%t(YRZ)
			YVYt = t(YVYt+CAind)+CAind - sum(solve(A2,ZtRXb)*ZtRXb)
			YVYt = t(YVYt/sqrt(s2))/sqrt(s2)/J
			cat("Eigen decomp")
			PL = Eigen(YVYt,nLatent+5)
			
			#YtmbR = t( t( (as.matrix(t(Y)) - Xb)*R )/sqrt(s2) )/sqrt(J)
			#YtmbRZ2 = t(YtmbR)%*%Z[,1:M2]
			#A2 = t(Z[,1:M2])%*%Z[,1:M2] + diag(1/u[1:M2]^2)
			#cat("Eigen decomp")
			#PL = Eigen(t(YtmbR)%*%YtmbR-YtmbRZ2%*%solve(A2)%*%t(YtmbRZ2),nLatent+5)
			#PL = Eigen(t(YtmbR)%*%YtmbR-YtmbRZ2%*%solve(A2)%*%t(YtmbRZ2),nLatent+5)

			flage=(apply(PL$evec^2,2,max)/apply(PL$evec^2,2,sum))<0.3
			cat("good pcs=");print(sum(flage[1:nLatent]))
			lambda = PL$eval[1:nLatent]; lambda[lambda<1]=1
			#Z[,(M2+1):M] = YtmbR%*%PL$evec[,1:nLatent]%*%diag(sqrt(1-1/lambda))
			Z[,(M2+1):M] = ( as.matrix((Yt*R)%*%(PL$evec[,1:nLatent]/sqrt(s2))) - RXb%*%t(colSums(PL$evec[,1:nLatent]/sqrt(s2))) )     %*%diag(sqrt(1-1/lambda))
			ZtZ = t(Z)%*%Z
			A = ZtZ+diag(1/u^2)
			ZtRX  = t(Z)%*%(X*R) # for W
			ZtRYt = t(as.matrix(Y%*%(Z*R))) # for W
		}
		
		# Var(xi)
		Vx = ZtZ - ZtZ %*% solve(A, ZtZ)
		hess = diag(psi)%*%t(H)%*%Kurt(Vx)%*%H%*%diag(psi)*2/N + t(H)%*%(Vx * t(solve(A)/u)/u)%*%H*2 - diag(c(t(H)%*%diag(Vx)))*2

		# rho
		if(length(rho)>1){print("rho")
			B=diag(L)
			ZtYtmbAll = NULL
			for(l in seq(L)){
				Xb     = c(X[batch==l,]%*%beta)
				Yl=Y[,batch==l]
				B[l,l] = sum(rowSums(Yl^2)/s2) - 2*sum((Yl%*%Xb)/s2) + sum(Xb^2)*sum(1/s2) # tauj * ytt %*% yt
				ZtYt   = t(as.matrix(Yl%*%Z[batch==l,]))
				ZtXb   = c( t(Z[batch==l,])%*%Xb )
				ZtYtmbAll = c(ZtYtmbAll, list(ZtYt - ZtXb))
			}
			for(l in seq(L)){
				for(ll in c(l:L)){
					B[l,ll] = B[l,ll] - sum( (solve(A, ZtYtmbAll[[l]])*ZtYtmbAll[[ll]])%*%(1/s2) )
					B[ll,l] = B[l,ll]
				}
			}
			rho = MaxRho0(B, as.numeric(table(batch)), J, rho)
			rho = rho/sum(rho)*length(rho)
			cat("rho=");print(rho)
		}

		# likelihood
		res = sum(log(s2))*N/2 + J * logDet(diag(M) + t(ZtZ*u)*u)/2
		
		print("psi")
		# psi | beta delta Psi(Z) sigma2
		ZtRYtmb = (ZtRYt - c(ZtRX%*%beta)) 
		AinvZtRYtmb = solve(A, ZtRYtmb)
		attr(res, "gradient") = - (c( t(H)%*%apply(t(t(AinvZtRYtmb^2)/s2)/u^3,1,sum) ) - J * c( t(H)%*%(diag(solve(A,ZtZ))/u) ))
		attr(res, "hessian") = - (hess+t(hess))/2 * J
		attr(res, "beta") = beta
		attr(res, "omega") = omega
		attr(res, "theta") = theta
		attr(res, "psi") = psi
		attr(res, "rho") = rho
		attr(res, "sigma2") = s2
		attr(res, "Z") = Z
		attr(res, "X") = X
		attr(res, "H") = H
		attr(res, "delta") = delta
		attr(res, "PL") = PL
		res
	}
	
	# Matrix prep
	lkhd.all=NULL
	res.min=NULL
	if(is.null(psi)){psi=rep(0.01,length(nh)); if(nLatent>0){psi[length(psi)]=1}}
	if(is.null(beta)){beta=rep(0.1,ncol(X))}
	if(is.null(delta)){delta=rep(0,nrow(Z))}
	if(is.null(omega) && !is.null(nexpcells)){
		s2 = apply(Y,1,var)
		omega = coef(lm(s2~0+nexpcells))
		lambda = 1/c(nexpcells%*%omega)
		theta = mean(lambda^2/(1/s2-lambda)^2)
		print(c(omega,theta))
	}
	Y  = Y[,subset]
	Yt = t(Y)
	YYt=NULL; YX=NULL; YZ=NULL;
	if(length(rho)==1 && nLatent>0){print("matrix prep")
		YYt=as.matrix(Y%*%Yt)
		YX=as.matrix(Y%*%X)
		YZ=as.matrix(Y%*%Z)
	}

	# Iteration start
	for(i in 1:ITRMAX){
		cat(paste("[",i,"]"))
		tmp = lkhd(beta=beta, psi=psi, delta=delta, omega=omega, theta=theta, s0=s0, Y=Y, Z=Z, X=X, H=H, rho=rho, batch=batch, Yt=Yt, YYt=YYt, YZ=YZ, YX=YX)
		#return(tmp)
		if(is.null(lkhd.all) || tmp<min(lkhd.all)){res.min=tmp}
		lkhd.all=c(lkhd.all,as.numeric(tmp))
		if(PLOT)plot(lkhd.all)
		if(nLatent>0){
			K2 = ncol(H)-1
			psi[1:K2] = psi[1:K2] - solve(attr(tmp,"hessian")[1:K2,1:K2]+2*diag(0.00001/psi[1:K2]^2-1),attr(tmp,"gradient")[1:K2]+2*(-0.00001/psi[1:K2]+psi[1:K2]))
		}else{
			psi = psi - solve(attr(tmp,"hessian")+2*diag(0.00001/psi^2-1),attr(tmp,"gradient")+2*(-0.00001/psi+psi))
		}
		cat("psi=");print(psi)
		beta =attr(tmp, "beta")
		delta =attr(tmp, "delta")
		omega =attr(tmp, "omega")
		theta =attr(tmp, "theta")
		rho =attr(tmp, "rho")
		cat("omega=");print(c(omega,theta))	
		Z = attr(tmp, "Z") # Psi
		if(length(lkhd.all)>50){lkhd.all=rev(rev(lkhd.all)[1:50])}
		if(i>1 & abs(diff(rev(lkhd.all)[1:2])/rev(lkhd.all)[1])<1e-7){cat("Converged\n");break}
	}
	c(list(lkhd=-as.numeric(res.min)), attributes(res.min), list(nh=nh))
}
