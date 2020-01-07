##### code modified from metaMA

#' effect size from moderate t test
#' @param tstat vector; t value of limma's moderate t test
#' @param ntilde doule; .
#' @param m doule; . 
#' @details calculate the effect size
#' @return a matrix 
#' @export
effectsize <- function(tstat,ntilde,m)
{
	# cm=gamma(m/2)/(sqrt(m/2)*gamma((m-1)/2))
	ln.cm<-lgamma(m/2)-log(sqrt(m/2))-lgamma((m-1)/2)  
	cm<-exp(ln.cm)
	d=tstat/sqrt(ntilde)
	dprime=cm*d
	terme1=m/((m-2)*ntilde)
	vard=terme1+d^2*(terme1*ntilde-1/cm^2)
	vardprime=cm^2*(terme1+dprime^2*(terme1*ntilde-1/cm^2))
	result=cbind(d,vard,dprime,vardprime)
	colnames(result)=c("d","vard","dprime","vardprime")
	result
}

#' directional Effect Size combination
#' @param ES data frame or matrix; row for genes, column for studies
#' @param varES data frame or matrix; row for genes, column for studies
#' @param useREM logical; . [default TRUE]
#' @importFrom plyr ldply
#' @details calculate the weight averaged effect size from multiple studies
#' @return a data.frame
#' @export
directEScombi <- function(ES,varES,useREM=TRUE)
{
	listres <- vector("list",2)
	f.Q.NA <- function (dadj, varadj)
	{
		#w <- 1/varadj
		#tmp1 <- w * dadj
		#mu <- rowSums(tmp1,na.rm=TRUE)/rowSums(w,na.rm=TRUE)
		#Q <- rowSums(w * (dadj - mu)^2,na.rm=TRUE)
		Q <- sapply(seq_len(nrow(dadj)),function(i){
					   f.use <- varadj[i,]!=0
					   w <- 1/varadj[i,f.use]
					   tmp1 <- w * dadj[i,f.use]
					   mu <- sum(tmp1,na.rm=TRUE)/sum(w,na.rm=TRUE)
					   sum(w * (dadj[i,f.use] - mu)^2,na.rm=TRUE)
					})
	}
	tau2.NA <- function(Q, num.studies, my.weights)
	{
			vwts <- rowSums(my.weights,na.rm=TRUE)
			tmp2 <- rowSums(my.weights^2,na.rm=TRUE)
			tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
			return(tau2)
	}
	tau2.NA.oneRow <- function(Q, num.studies, my.weights)
	{
		f.use <- is.finite(my.weights)
		vwts <- sum(my.weights[f.use],na.rm=TRUE)
		tmp2 <- sum(my.weights[f.use]^2,na.rm=TRUE)
		tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
		return(tau2)
	}
	Qvals <- f.Q.NA(ES, varES)

#	num.studies <- dim(ES)[2]
#	if (useREM) { varES <- varES + tau2.NA(Qvals, num.studies, my.weights = 1/varES) }
#	wt <- 1/varES
#	MUvals <- rowSums(ES * wt,na.rm=TRUE)/rowSums(wt,na.rm=TRUE)
#	MUsES <- sqrt(1/rowSums(wt,na.rm=TRUE))
#	zSco <- MUvals/MUsES
#	###rpvalESc <- 2*(1-pnorm(abs(zSco)))
#	rpvalESc <- 2*(pnorm(-abs(zSco)))
#	rpadjEsc <- p.adjust(rpvalESc,method="BH")
#	###res=which(p.adjust(rpvalESc,method="BH")<=BHth)
#	result=cbind(MUvals,MUsES,zSco,rpvalESc,rpadjEsc)

	result <- ldply(seq_len(nrow(varES)),function(i){
						ES.i <- ES[i,]
						varES.i <- varES[i,]
						f.use <- varES.i!=0
						if(useREM){
							varES.i <- varES.i + tau2.NA.oneRow(Qvals[i],sum(f.use),1/varES.i[f.use])
						}
						wt.i <- 1/varES.i[f.use]
						MUvals <- sum(ES.i[f.use] * wt.i,na.rm=TRUE)/sum(wt.i,na.rm=TRUE)
						MUsES <- sqrt(1/sum(wt.i,na.rm=TRUE))
						zSco <- MUvals/MUsES
						rpvalESc <- 2*(pnorm(-abs(zSco)))
						cbind(MUvals,MUsES,zSco,rpvalESc)
					})
	colnames(result) <- c("comb.ES","comb.ES.sd","comb.Z","comb.p")
	rownames(result) <- rownames(varES)
	result[,"comb.padj"] <- p.adjust(result[,"comb.p"],method="BH")

	result
}



