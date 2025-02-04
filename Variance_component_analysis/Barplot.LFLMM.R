Barplot <-
function(x,col=1){
	
	getUpper95=function(psi, V){
		denom = sum(psi^2)+1
		psivar=NULL
		for(k in 1:length(psi)){
			grad    =         - (2*psi)*psi[k]^2/denom^2
			grad[k] = grad[k] + (2*psi[k])/denom
			psivar = c(psivar, sum(t(V*grad)*grad))
		}
		
		psi^2/denom + 1.96 * sqrt(psivar)
	}

	remids=c(1)
	V = solve(x$hessian)[-remids,-remids]
	labs = names(x$nh)[-remids]
	psi  = x$psi[-remids]
	ord  = rev(order(psi^2))
	par(mgp=c(2,0.5,0))
	
	varexp = psi[ord]^2/(1+sum(psi^2))*100
	varse = getUpper95(psi,V)[ord]*100
	xcoords = barplot(varexp, name=labs[ord],las=2, ylim=c(0,max(varse)*1.01), col=col,border=col,ylab="Variance explained %",xaxt="n",cex.lab=1.5,cex.axis=1.2)
	text(cex=1, x=xcoords , y=par("usr")[3] - 1, labs[ord], xpd=NA, srt=45, adj = 0.965,cex.axis=1.2)
		
	segments(xcoords, varexp, xcoords, varse)
	segments(xcoords-0.3, varse, xcoords+0.3, varse)
	
	res = cbind(varexp, varse)
	rownames(res)=labs[ord]
	res
	
	#ord
}
