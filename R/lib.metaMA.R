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
#' @details calculate the weight averaged effect size from multiple studies
#' @return a matrix 
#' @export
directEScombi <- function(ES,varES,useREM=TRUE)
{
	listres=vector("list",2)
	f.Q.NA=function (dadj, varadj) 
	{
		w <- 1/varadj
		tmp1 <- w * dadj
		mu <- rowSums(tmp1,na.rm=TRUE)/rowSums(w,na.rm=TRUE)
		Q <- rowSums(w * (dadj - mu)^2,na.rm=TRUE)
	}
	tau2.NA <- function(Q, num.studies, my.weights)
	{
			vwts <- rowSums(my.weights,na.rm=TRUE)
			tmp2 <- rowSums(my.weights^2,na.rm=TRUE)
			tau2 <- pmax(0, (Q - (num.studies - 1))/(vwts - tmp2/vwts))
			return(tau2)
	}
	num.studies <- dim(ES)[2]
	Qvals <- f.Q.NA(ES, varES)
	if (useREM) { varES <- varES + tau2.NA(Qvals, num.studies, my.weights = 1/varES) }
	wt <- 1/varES
	MUvals <- rowSums(ES * wt,na.rm=TRUE)/rowSums(wt,na.rm=TRUE)
	MUsES <- sqrt(1/rowSums(wt,na.rm=TRUE))
	zSco <- MUvals/MUsES
	###rpvalESc <- 2*(1-pnorm(abs(zSco)))
	rpvalESc <- 2*(pnorm(-abs(zSco)))
	rpadjEsc <- p.adjust(rpvalESc,method="BH")
	###res=which(p.adjust(rpvalESc,method="BH")<=BHth)
	result=cbind(MUvals,MUsES,zSco,rpvalESc,rpadjEsc)
	colnames(result)=c("comb.ES","comb.ES.sd","comb.Z","comb.p","comb.padj")
	result
}



