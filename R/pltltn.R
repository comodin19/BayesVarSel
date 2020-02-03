#' Correction for p>>n for an object of class \code{Bvs}
#'
#' In cases where p>>n and the true model is expected to be sparse, it is very unlikely that the Gibbs sampling
#' will sample models in the singular subset of the model space (models with k>n). Nevertheless, depending on
#' how large is p/n and the strenght of the signal, this part of the model space could be very influential in the
#' final response. 
#' 
#' From an object created with GibbsBvs and prior probabilities specified as Scott-Berger, 
#' this function provides an estimation of the posterior probability of models with k>n which is a measure of the 
#' importance of these models. In summary, when this probability is large,  the sample size is not large enough to beat
#' such large p.
#' Additionally, \code{pltltn} gives corrections of the posterior inclusion probabilities and posterior probabilities 
#' of dimension of the true model.
#'
#' @export
#' @param object An object of class \code{Bvs} obtained with \code{GibbsBvs}
#' @return \code{pltltn} returns a list with the following elements: 
#' \item{pS }{An estimation of the probability that the true model is irregular (k>n)} 
#' \item{postprobdim }{A corrected estimation of the posterior probabilities over the dimensions}
#' \item{inclprob }{A corrected estimation of the posterior inclusion probabilities}
#' @author Gonzalo Garcia-Donato
#'
#' Maintainer: <gonzalo.garciadonato@uclm.es>
#' @seealso See 
#'   \code{\link[BayesVarSel]{GibbsBvs}} for creating objects of the class
#'   \code{Bvs}.
#' @examples
#'
#' \dontrun{
#' #Analysis of
#' }
#'
pltltn<- function(object){
	#Corrected posterior inclusion probabilities and probabilities of dimension
	#for the case where p>>n.
	#Ms=model space with singular models; Mr=model space with regular models
  if (!inherits(object, "Bvs"))
    stop("calling summary.Bvs(<fake-Bvs-x>) ...")
	
	if (object$method != "gibbs") 
		stop("This corrected estimates are for Gibbs sampling objects.")
	
	if (object$priorprobs!="ScottBerger"){
		stop("This function was conceived to be used in combination with Scott-Berger prior\n")
	}
	
	cat("Info: Use this function only for problems with p>>n (not just p>n)\n")
	
	#The method weitghs results on Mr (provided by Gibbs sampling) with those in Ms (theoretical)
	
	#first obtain the estimates conditionall on Mr
	kgamma<- rowSums(object$modelslogBF[,1:object$p])
	isMr<- kgamma < object$n-length(object$lmnull$coefficients)
	inclprobMr<- colMeans(object$modelslogBF[isMr,1:object$p])
	postprobdimMr<- table(kgamma[isMr])/sum(isMr)
	names(postprobdimMr)<- as.numeric(names(postprobdimMr))+length(object$lmnull$coefficients)
		
	#ratio of prior probabilities of Ms to Mr
	qSR<- (object$p-object$n)/(object$n+1)
	#The posterior probabililty of Ms:
	pS<- qSR/(qSR+object$C)
	cat(paste("Estimate of the posterior probability of the\n model space with singular models is:", round(pS,3),"\n"))

	#Corrected posterior probability over the dimension:
	postprobdim<- postprobdimMr*(1-pS)
	#Corrected inclusion probabilities:
	inclprob<- inclprobMr*(1-pS) + 0.5*pS
	result<- list()
	result$pS<- pS
	result$postprobdim<- postprobdim
	result$inclprob<- inclprob
	return(result)
}

