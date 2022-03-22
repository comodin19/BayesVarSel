#' Bayesian Variable Selection for linear regression models
#'
#' Exact computation of summaries of the posterior distribution using
#' sequential computation.
#'
#' The model space is the set of all models, Mi, that contain the intercept and
#' are nested in that specified by \code{formula}. The simplest of such models,
#' M0, contains only the intercept. Then \code{Bvs} provides exact summaries of
#' the posterior distribution over this model space, that is, summaries of the
#' discrete distribution which assigns to each model Mi its probability given
#' the data:
#'
#' Pr(Mi | \code{data})=Pr(Mi)*Bi/C,
#'
#' where Bi is the Bayes factor of Mi to M0, Pr(Mi) is the prior probability of
#' Mi and C is the normalizing constant.
#'
#' The Bayes factor B_i depends on the prior assigned for the parameters in the regression
#' models Mi and \code{Bvs} implements a number of popular choices. The "Robust"
#' prior by Bayarri, Berger, Forte and Garcia-Donato (2012) is the recommended (and default) choice.
#' This prior prior can be implemented in a more stable way using the derivations in Greenaway (2019) 
#' and that are available in BayesVarSel since version 2.2.x setting the argument to "Robust.G". 
#'
#' Additional options are "gZellner" a prior which 
#' corresponds to the
#' prior in Zellner (1986) with g=n. Also "Liangetal" prior is the
#' hyper-g/n with a=3 (see the original paper Liang et al 2008, for details).
#' "ZellnerSiow" is the multivariate Cauchy prior proposed by Zellner and Siow
#' (1980, 1984), further studied by Bayarri and Garcia-Donato (2007).
#' "FLS" is the (benchmark) prior recommended by Fernandez, Ley and Steel (2001) which is
#' the prior in Zellner (1986) with g=max(n, p*p) p being the number of
#' covariates to choose from (the most complex model has p+number of fixed
#' covariates).
#' "intrinsic.MGC" is the intrinsic prior derived by Moreno, Giron, Casella (2015)
#' and "IHG" corresponds to the intrinsic hyper-g prior derived in Berger, Garcia-Donato, 
#' Moreno and Pericchi (2022).
#'
#' With respect to the prior over the model space Pr(Mi) three possibilities
#' are implemented: "Constant", under which every model has the same prior
#' probability, "ScottBerger" under which Pr(Mi) is inversely proportional to
#' the number of models of that dimension, and "User" (see below). The
#' "ScottBerger" prior was studied by Scott and Berger (2010) and controls for
#' multiplicity (default choice since version 1.7.0).
#'
#' When the parameter \code{prior.models}="User", the prior probabilities are
#' defined through the p+1 dimensional parameter vector \code{priorprobs}. Let
#' k be the number of explanatory variables in the simplest model (the one
#' defined by \code{fixed.cov}) then except for the normalizing constant, the
#' first component of \code{priorprobs} must contain the probability of each
#' model with k covariates (there is only one); the second component of
#' \code{priorprobs} should contain the probability of each model with k+1
#' covariates and so on. Finally, the p+1 component in \code{priorprobs}
#' defined the probability of the most complex model (that defined by
#' \code{formula}. That is
#'
#' \code{priorprobs}[j]=Cprior*Pr(M_i such that M_i has j-1+k explanatory variables)
#'
#' where Cprior is the normalizing constant for the prior, i.e
#' \code{Cprior=1/sum(priorprobs*choose(p,0:p)}.
#'
#' Note that \code{prior.models}="Constant" is equivalent to the combination
#' \code{prior.models}="User" and \code{priorprobs=rep(1,(p+1))} but the
#' internal functions are not the same and you can obtain small variations in
#' results due to these differences in the implementation.
#'
#' Similarly, \code{prior.models} = "ScottBerger" is equivalent to the
#' combination \code{prior.models}= "User" and \code{priorprobs} =
#' \code{1/choose(p,0:p)}.
#'
#' The case where n<p is handled assigning to the Bayes factors of models with k regressors
#' with n<k a value of 1. This should be interpreted as a generalization of the
#' null predictive matching in Bayarri et al (2012). Use \code{\link[BayesVarSel]{GibbsBvs}} for cases where
#' p>>.
#' 
#' Limitations: about the error "A Bayes factor is infinite.". Bayes factors can be
#' extremely big numbers if i) the sample size is large or if
#' ii) a competing model is much better (in terms of fit) than the model taken as the
#' null model. If you see this error, try to use the more stable version of the
#' robust prior "Robust.g" and/or reconisder using more accurate and realistic definitions
#' of the simplest model (via the \code{null.model} argument).
#'
#'
#' @export
#' @param formula Formula defining the most complex (full) regression model in the
#' analysis. See details.
#' @param data data frame containing the data.
#' @param null.model A formula defining which is the simplest (null) model.
#' It should be nested in the full model. By default, the null model is defined
#' to be the one with just the intercept.
#' @param prior.betas Prior distribution for regression parameters within each
#' model (to be literally specified). Possible choices include "Robust", "Robust.G", "Liangetal", "gZellner",
#' "ZellnerSiow", "FLS", "intrinsic.MGC" and "IHG" (see details).
#' @param prior.models Prior distribution over the model space (to be literally specified). Possible
#' choices are "Constant", "ScottBerger" and "User" (see details).
#' @param n.keep How many of the most probable models are to be kept? By
#' default is set to 10, which is automatically adjusted if 10 is greater than
#' the total number of models.
#' @param time.test If TRUE and the number of variables is moderately large
#' (>=18) a preliminary test to estimate computational time is performed.
#' @param priorprobs A p+1 (p is the number of non-fixed covariates)
#' dimensional vector defining the prior probabilities Pr(M_i) (should be used
#' in the case where \code{prior.models}= "User"; see details.)
#' @param parallel A logical parameter specifying whether parallel computation
#' must be used (if set to TRUE)
#' @param n.nodes The number of cores to be used if parallel computation is used.
#' @return \code{Bvs} returns an object of class \code{Bvs} with the following
#' elements: \item{time }{The internal time consumed in solving the problem}
#' \item{lmfull }{The \code{lm} class object that results when the model
#' defined by \code{formula} is fitted by \code{lm}}
#' \item{lmnull }{The
#' \code{lm} class object that results when the model defined by
#' \code{null.model} is fitted by \code{lm}}
#' \item{variables }{The name of all
#' the potential explanatory variables (the set of variables to select from).}
#' \item{n }{Number of observations} \item{p }{Number of explanatory variables
#' to select from} \item{k }{Number of fixed variables} \item{HPMbin }{The
#' binary expression of the Highest Posterior Probability model}
#' \item{modelsprob }{A \code{data.frame} which summaries the \code{n.keep}
#' most probable, a posteriori models, and their associated probability.}
#' \item{inclprob }{A named vector with the inclusion probabilities of all
#' the variables.} \item{jointinclprob }{A \code{data.frame} with the joint
#' inclusion probabilities of all the variables.} \item{postprobdim }{Posterior
#' probabilities of the dimension of the true model} \item{call }{The
#' \code{call} to the function} 
#' \item{C}{The value of the normalizing constant (C=sum BiPr(Mi), for Mi in the model space)}
#' \item{method }{\code{full} or \code{parallel} in case of
#' parallel computation}
#' \item{prior.betas}{prior.betas}
#' \item{prior.models}{prior.models}
#' \item{priorprobs}{priorprobs}
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#' Maintainer: <anabel.forte@@uv.es>
#' @seealso Use \code{\link[BayesVarSel]{print.Bvs}} for the best visited models and an 
#' estimation of their posterior probabilities and  \code{\link[BayesVarSel]{summary.Bvs}} for
#' summaries of the posterior distribution.
#' 
#' \code{\link[BayesVarSel]{plot.Bvs}} for several plots of the result,
#' \code{\link[BayesVarSel]{BMAcoeff}} for obtaining model averaged simulations
#' of regression coefficients and \code{\link[BayesVarSel]{predict.Bvs}} for
#' predictions.
#'
#' \code{\link[BayesVarSel]{GibbsBvs}} for a heuristic approximation based on
#' Gibbs sampling (recommended when p>20, no other possibilities when p>31).
#' 
#' See \code{\link[BayesVarSel]{GibbsBvsF}} if there are factors among the explanatory variables
#' @references Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1557.
#'
#' Berger, J., Garcıa-Donato, G., Moreno, E., and Pericchi, L. (2022). 
#' The intrinsic hyper-g prior for normal linear models. in preparation.
#' 
#' Bayarri, M.J. and Garcia-Donato, G. (2007)<DOI:10.1093/biomet/asm014>
#' Extending conventional priors for testing general hypotheses in linear
#' models. Biometrika, 94:135-152.
#'
#' Barbieri, M and Berger, J (2004)<DOI:10.1214/009053604000000238> Optimal
#' Predictive Model Selection. The Annals of Statistics, 32, 870-897.
#'
#' Fernandez, C., Ley, E. and Steel, M.F.J.
#' (2001)<DOI:10.1016/s0304-4076(00)00076-2> Benchmark priors for Bayesian
#' model averaging. Journal of Econometrics, 100, 381-427.
#'
#' Greenaway, M. (2019) Numerically stable approximate Bayesian methods for
#' generalized linear mixed models and linear model selection. Thesis (Department of Statistics,
#' University of Sydney).
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger,J.O.
#' (2008)<DOI:10.1198/016214507000001337> Mixtures of g-priors for Bayesian
#' Variable Selection. Journal of the American Statistical Association.
#' 103:410-423
#'
#' Moreno, E., Giron, J. and Casella, G. (2015) Posterior model consistency
#' in variable selection as the model dimension grows. Statistical Science. 30: 228-241.
#'
#' Zellner, A. and Siow, A. (1980)<DOI:10.1007/bf02888369> Posterior Odds Ratio
#' for Selected Regression Hypotheses. In Bayesian Statistics 1 (J.M. Bernardo,
#' M. H. DeGroot, D. V. Lindley and A. F. M. Smith, eds.) 585-603. Valencia:
#' University Press.
#'
#' Zellner, A. and Siow, A. (1984). Basic Issues in Econometrics. Chicago:
#' University of Chicago Press.
#'
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#' @keywords package
#' @examples
#'
#' \dontrun{
#' #Analysis of Crime Data
#' #load data
#' data(UScrime)
#'
#' #Default arguments are Robust prior for the regression parameters
#' #and constant prior over the model space
#' #Here we keep the 1000 most probable models a posteriori:
#' crime.Bvs<- Bvs(formula= y ~ ., data=UScrime, n.keep=1000)
#'
#' #A look at the results:
#' crime.Bvs
#'
#' summary(crime.Bvs)
#'
#' #A plot with the posterior probabilities of the dimension of the
#' #true model:
#' plot(crime.Bvs, option="dimension")
#'
#' #Two image plots of the conditional inclusion probabilities:
#' plot(crime.Bvs, option="conditional")
#' plot(crime.Bvs, option="not")
#' }
#'
Bvs <-
  function(formula,
           data,
           null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
           prior.betas = "Robust",
           prior.models = "ScottBerger",
           n.keep = 10,
           time.test = TRUE,
           priorprobs = NULL,
           parallel = FALSE,
           n.nodes = detectCores()) {

    formula <- as.formula(formula)
    null.model<- as.formula(null.model)

    #The response in the null model and in the full model must coincide
    if (formula[[2]] != null.model[[2]]){
      stop("The response in the full and null model does not coincide.\n")
    }

    #Let's define the result
    result <- list()

    #at least two nodes are needed in the parallel
    if (parallel) {
      if (n.nodes < 2) {
        stop("At least 2 nodes are needed\n")
      }
    }

    #Get a tempdir as working directory
    wd <- tempdir()
    #remove all the previous documents in the working directory
    unlink(paste(wd, "*", sep = "/"))

    #evaluate the null model:
    lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)
    fixed.cov <- dimnames(lmnull$x)[[2]]

    #eval the full model
    #Set the design matrix if fixed covariates present:
    if (!is.null(fixed.cov)) {
      lmfull <- lm(formula, data, y = TRUE, x = TRUE)
      X.full <- lmfull$x
      namesx <- dimnames(X.full)[[2]]

      #check if null model is contained in the full one:

      namesnull <- dimnames(lmnull$x)[[2]]
      "%notin%" <- function(x, table) match(x, table, nomatch = 0) == 0
      for (i in 1:length(namesnull)){
        if (namesnull[i] %notin% namesx) {
          cat("Error in var: ", namesnull[i], "\n")
          stop("null model is not nested in full model\n")
        }
      }

      #Is there any variable to select from?
      if (length(namesx) == length(namesnull)) {
        stop(
          "The number of fixed covariates is equal to the number of covariates in the full model. No model selection can be done\n"
        )
      }

      #position for fixed variables in the full model
      fixed.pos <- which(namesx %in% namesnull)

      n <- dim(data)[1]
			
			#warning about n<p situation
			if (n < dim(X.full)[2]) cat("In this dataset n<p and unitary Bayes factors are used for models with k>n.\n")

      #the response variable for the C code
      Y <- lmnull$residuals

      #Design matrix of the null model
      X0 <- lmnull$x
      P0 <-
        X0 %*% (solve(t(X0) %*% X0)) %*% t(X0)#Intentar mejorar aprovechando lmnull
      knull <- dim(X0)[2]

      #matrix containing the covariates from which we want to select
      X1 <- lmfull$x[, -fixed.pos]
      if (dim(X1)[1] < n) {
        stop("NA values found for some of the competing variables")
      }

      #Design matrix for the C-code
      X <- (diag(n) - P0) %*% X1 #equivalent to X<- (I-P0)X
      namesx <- dimnames(X)[[2]]
      if (namesx[1] == "(Intercept)") {
        namesx[1] <-
          "Intercept" #namesx contains the name of variables including the intercept
      }
			
	    #Check if, among the competing variables, there are factors
			if (sum(attr(lmfull$terms, "dataClasses")=="factor") & sum(attr(lmfull$terms, "dataClasses")=="factor")-sum(attr(lmnull$terms, "dataClasses")=="factor")>0){
				cat("--------------\n")
				cat("The competing variables contain factors, a situation for which we recommend using\n")
				cat("GibbsBvsF()\n")
	      ANSWER <-
	        readline("Do you want to continue with Bvs()?(y/n) then press enter.\n")
	      while (substr(ANSWER, 1, 1) != "n" &
	             substr(ANSWER, 1, 1) != "y") {
	        ANSWER <- readline("")
	      }

	      if (substr(ANSWER, 1, 1) == "n")
	      {
	        return(NULL)
	      }        
			}

      p <- dim(X)[2] #Number of covariates to select from

      #check if the number of models to save is correct
      if (!parallel & n.keep > 2 ^ (p)) {
        warning(
          paste(
            "The number of models to keep (",
            n.keep,
            ") is larger than the total number of models (",
            2 ^ (p),
            ") and it has been set to ",
            2 ^ (p) ,
            sep = ""
          )
        )
        n.keep <- 2 ^ p
      }
      if (parallel & n.keep > 2 ^ (p) / n.nodes) {
        warning(
          paste(
            "The number of models to keep (",
            n.keep,
            ") must be smaller than the total number of models per node(",
            2 ^ (p) / n.nodes,
            ") and it has been set to ",
            2 ^ (p) / n.nodes ,
            sep = ""
          )
        )
        n.keep <- 2 ^ p  / n.nodes
      }


    }

    #If no fixed covariates considered
    if (is.null(fixed.cov)) {
      #Check that all the fixed covariables are included in the full model
      lmfull = lm(formula, data, y = TRUE, x = TRUE)
      X.full <- lmfull$x
      namesx <- dimnames(X.full)[[2]]

      #remove the brackets in "(Intercept)" if present.
      if (namesx[1] == "(Intercept)") {
        namesx[1] <-
          "Intercept" #namesx contains the name of variables including the intercept
      }


      X <- lmfull$x
      knull <- 0
      Y <- lmfull$y
      p <- dim(X)[2]
      n <- dim(X)[1]
      #check if the number of models to save is correct
      if (!parallel & n.keep > 2 ^ (p)) {
        warning(
          paste(
            "The number of models to keep (",
            n.keep,
            ") is larger than the total number of models (",
            2 ^ (p),
            ") and it has been set to ",
            2 ^ (p) ,
            sep = ""
          )
        )
        n.keep <- 2 ^ p
      }
      if (parallel & n.keep > 2 ^ (p) / n.nodes) {
        warning(
          paste(
            "The number of models to keep (",
            n.keep,
            ") must be smaller than the total number of models per node(",
            2 ^ (p) / n.nodes,
            ") and it has been set to ",
            2 ^ (p) / n.nodes ,
            sep = ""
          )
        )
        n.keep <- 2 ^ p  / n.nodes
      }

    }

    #write the data files in the working directory
    write(Y,
          ncolumns = 1,
          file = paste(wd, "/Dependent.txt", sep = ""))
    write(t(X),
          ncolumns = p,
          file = paste(wd, "/Design.txt", sep = ""))

    #prior for betas: (#pfb will contain the internal representation of the prior.betas argument)
	pfb <- check.prior.betas(prior.betas) 
    #prior for model space: (#pfms will contain the internal representation of the prior.models argument)
	#this function also writes the vector of prior probabilities for usage of C functions
	pfms <- check.prior.modelspace(prior.models, priorprobs, p, wd)
 
    method <- paste(pfb, pfms, sep = "")

    #Info:
    cat("Info. . . .\n")
    cat("Most complex model has", p + knull, "covariates\n")
    if (!is.null(fixed.cov)) {
      if (knull > 1) {
        cat("From those",
            knull,
            "are fixed and we should select from the remaining",
            p,
            "\n")
      }
      if (knull == 1) {
        cat("From those",
            knull,
            "is fixed and we should select from the remaining",
            p,
            "\n")
      }
      cat(paste(paste(
        namesx, collapse = ", ", sep = ""
      ), "\n", sep = ""))
    }
    else {
      cat("None is fixed and we should select from the remaining",
          p,
          "\n")

      cat(paste(paste(
        namesx, collapse = ", ", sep = ""
      ), "\n", sep = ""))
    }

    cat("The problem has a total of", 2 ^ (p), "competing models\n")
    cat("Of these, the ",
        n.keep,
        "most probable (a posteriori) are kept\n")


    #check if the number of covariates is too big.
    if (p > 30) {
      stop("Number of covariates too big. . . consider using GibbsBvs\n")
    }


    myfun <- function(name.start.end, method) {
      #if name.start.end has length 3 it comes from the parallel (which needs a suffix to differentiate)
      #each of the processes and is a number in the first position in this vector. Otherwise
      #it has dimension two. To reconciliate both in a common function do the following

      #non-parallel case:
      if (length(name.start.end) == 2){
        name.start.end<- c("", name.start.end)
      }

      Cresult <- switch(method,
                        "gc" = .C(
                          "gConst",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "gs" = .C(
                          "gSB",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "gu" = .C(
                          "gUser",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "rc" = .C(
                          "RobustConst",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "rs" = .C(
                          "RobustSB",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "ru" = .C(
                          "RobustUser",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "lc" = .C(
                          "LiangConst",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "ls" = .C(
                          "LiangSB",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "lu" = .C(
                          "LiangUser",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "zc" = .C(
                          "ZSConst",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "zs" = .C(
                          "ZSSB",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "zu" = .C(
                          "ZSUser",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "fc" = .C(
                          "flsConst",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "fs" = .C(
                          "flsSB",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "fu" = .C(
                          "flsUser",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "iu" = .C(
                          "intrinsicUser",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "is" = .C(
                          "intrinsicSB",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "ic" = .C(
                          "intrinsicConst",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),						
                        "xu" = .C(
                          "geointrinsicUser",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "xs" = .C(
                          "geointrinsicSB",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "xc" = .C(
                          "geointrinsicConst",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),												
                        "yu" = .C(
                          "geointrinsic2User",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "ys" = .C(
                          "geointrinsic2SB",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "yc" = .C(
                          "geointrinsic2Const",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "wu" = .C(
                          "Robust2User",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "ws" = .C(
                          "Robust2SB",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        ),
                        "wc" = .C(
                          "Robust2Const",
                          as.character(name.start.end[1]),
                          as.integer(n),
                          as.integer(p),
                          as.integer(n.keep),
                          as.integer(name.start.end[2]),
                          as.integer(name.start.end[3]),
                          as.character(wd),
                          as.double(estim.time),
                          as.integer(knull)
                        )
      )
      return(Cresult)
    }

    #The previous test (for time)
    estim.time <- 0

    if (time.test && p >= 18) {
      cat("Time test. . . .\n")
      smallset <- c(2 ^ (p - 1) - 1999, (2 ^ (p - 1) + 2000))
      Cresult <- myfun(smallset, method = method)

      if (parallel){
        estim.time <- Cresult[[8]] * 2 ^ (p) / (60 * 4000 * n.nodes)
      }
      else estim.time <- Cresult[[8]] * 2 ^ (p) / (60 * 4000)

      cat("The problem would take ",
          estim.time,
          "minutes (approx.) to run\n")
      ANSWER <-
        readline("Do you want to continue?(y/n) then press enter.\n")
      while (substr(ANSWER, 1, 1) != "n" &
             substr(ANSWER, 1, 1) != "y") {
        ANSWER <- readline("")
      }

      if (substr(ANSWER, 1, 1) == "n")
      {
        return(NULL)
      }

    }

    #if the answer is yes work on the problem
    cat("Working on the problem...please wait.\n")

    if (parallel){

      #Calculate how to distribute the model space through the nodes:

      iterperproc <- round((2 ^ (p) - 1) / n.nodes)
      if (n.keep > iterperproc)
        stop("Number of kept models should be smaller than the number of models per node\n")
      distrib <- list()
      for (i in 1:(n.nodes - 1)) {
        distrib[[i]] <- c(i, (i - 1) * iterperproc + 1, i * iterperproc)
      }
      distrib[[n.nodes]] <-
        c(n.nodes, (n.nodes - 1) * iterperproc + 1, 2 ^ (p) - 1)

      cl <- makeCluster(n.nodes)
      #Load the library in the different nodes
      clusterEvalQ(cl, library(BayesVarSel))

      clusterApply(cl, distrib, myfun, method = method)
      stopCluster(cl)

      ##############Put together the results

      #next is the prior probability for the null model Pr(M_0)=p_0/sum(p_j)
      if (pfms == "c") {
        PrM0 <- 1 / 2 ^ (p - 1)
        #the unnormalized prior prob for M0:
        p0 <- 1
      }

      if (pfms == "s") {
        PrM0 <- 1 / (p + 1)
        #the unnormalized prior prob for M0:
        p0 <- 1
      }

      if (pfms == "u") {
        PrM0 <- priorprobs[1] / sum(choose(p, 0:p) * priorprobs)
        #the unnormalized prior prob for M0:
        p0 <- priorprobs[1]
      }

      fPostProb <- paste(wd, "PostProb", sep = "/")
      fInclusionProb <- paste(wd, "InclusionProb", sep = "/")
      fMostProbModels <- paste(wd, "MostProbModels", sep = "/")
      fNormConstant <- paste(wd, "NormConstant", sep = "/")
      fNormConstantPrior <- paste(wd, "NormConstantPrior", sep = "/")
      fProbDimension <- paste(wd, "ProbDimension", sep = "/")
      fJointInclusionProb <- paste(wd, "JointInclusionProb", sep = "/")
      fBetahat <- paste(wd, "betahat", sep = "/")

      #Obtain the normalizing constant (say E) for the prior probabilities:
      #Pr(Ml)=p_l/E
      E <- 0
      for (i in 1:n.nodes) {
        E <- E + scan(
          file = paste(fNormConstantPrior, i, sep = ""),
          n = 1,
          quiet = T
        )
      }
      E <- E - (n.nodes - 1) * p0

      #Obtain the normalizing constant (say D) for the posterior probabilities:
      #Pr(Ml|data)=B_{l0}*Pr(M_l)/D, where B_{l0}=m_l(data)/m_0(data)
      D <- 0
      for (i in 1:n.nodes) {
        D <-
          D + scan(
            file = paste(fNormConstant, i, sep = ""),
            n = 1,
            quiet = T
          ) * scan(
            file = paste(fNormConstantPrior, i, sep = ""),
            n = 1,
            quiet = T
          )
      }
      D <- (D - (n.nodes - 1) * PrM0) / E


      #Now obtain the n.keep most probable models
      i <- 1
      thisNormConstant <-
        scan(
          file = paste(fNormConstant, i, sep = ""),
          n = 1,
          quiet = T
        )
      thisNormConstantPrior <-
        scan(file = paste(fNormConstantPrior, i, sep = ""),
             quiet = T)
      #next is the Bayes factor times the (unnormalized) prior for this model
      #(see the main.c code to see how is the unnormalized prior). So, if
      #the unnormalized prior is=1, then next is the Bayes factor*1 and so on
      thisUnnorPostProb <-
        read.table(file = paste(fPostProb, i, sep = ""),
                   colClasses = "numeric")[[1]] * thisNormConstant * thisNormConstantPrior

      thisMostProbModels <-
        read.table(file = paste(fMostProbModels, 1, sep = ""),
                   colClasses = "numeric")[[1]]

      for (i in 2:n.nodes) {
        readNormConstant <-
          scan(
            file = paste(fNormConstant, i, sep = ""),
            n = 1,
            quiet = T
          )
        readNormConstantPrior <-
          scan(file = paste(fNormConstantPrior, i, sep = ""),
               quiet = T)
        readUnnorPostProb <-
          read.table(file = paste(fPostProb, i, sep = ""),
                     colClasses = "numeric")[[1]] * readNormConstant * readNormConstantPrior
        readMostProbModels <-
          read.table(file = paste(fMostProbModels, i, sep = ""),
                     colClasses = "numeric")[[1]]

        jointUnnorPostProb <- c(readUnnorPostProb, thisUnnorPostProb)
        jointModels <- c(readMostProbModels, thisMostProbModels)
        reorder <- order(jointUnnorPostProb, decreasing = T)
        thisUnnorPostProb <- jointUnnorPostProb[reorder[1:n.keep]]
        thisMostProbModels <- jointModels[reorder[1:n.keep]]
      }


      #The inclusion probabilities
      accum.InclusionProb <-
        read.table(file = paste(fInclusionProb, i, sep = ""))[[1]] * 0
      for (i in 1:n.nodes) {
        readNormConstant <-
          scan(
            file = paste(fNormConstant, i, sep = ""),
            n = 1,
            quiet = T
          )
        readNormConstantPrior <-
          scan(file = paste(fNormConstantPrior, i, sep = ""),
               quiet = T)
        accum.InclusionProb <-
          accum.InclusionProb + read.table(file = paste(fInclusionProb, i, sep =
                                                          ""),
                                           colClasses = "numeric")[[1]] * readNormConstant * readNormConstantPrior
      }

      accum.InclusionProb <- accum.InclusionProb / (D * E)

      #The joint inclusion probs:
      accum.JointInclusionProb <-
        as.matrix(read.table(
          file = paste(fJointInclusionProb, i, sep = ""),
          colClasses = "numeric"
        )) * 0
      for (i in 1:n.nodes) {
        readNormConstant <-
          scan(
            file = paste(fNormConstant, i, sep = ""),
            n = 1,
            quiet = T
          )
        readNormConstantPrior <-
          scan(file = paste(fNormConstantPrior, i, sep = ""),
               quiet = T)
        accum.JointInclusionProb <- accum.JointInclusionProb +
          as.matrix(read.table(
            file = paste(fJointInclusionProb, i, sep = ""),
            colClasses = "numeric"
          )) * readNormConstant * readNormConstantPrior
      }

      accum.JointInclusionProb <- accum.JointInclusionProb / (D * E)
      #-----

      #The dimension probabilities
      accum.ProbDimension <-
        read.table(file = paste(fProbDimension, i, sep = ""),
                   colClasses = "numeric")[[1]] * 0
      for (i in 1:n.nodes) {
        readNormConstant <-
          scan(
            file = paste(fNormConstant, i, sep = ""),
            n = 1,
            quiet = T
          )
        readNormConstantPrior <-
          scan(file = paste(fNormConstantPrior, i, sep = ""),
               quiet = T)
        accum.ProbDimension <-
          accum.ProbDimension + read.table(file = paste(fProbDimension, i, sep =
                                                          ""),
                                           colClasses = "numeric")[[1]] * readNormConstant * readNormConstantPrior
      }

      accum.ProbDimension <- accum.ProbDimension / (D * E)

      betahat <-
        read.table(file = paste(fBetahat, i, sep = ""), colClasses = "numeric")[[1]] *
        0
      ac <- 0
      for (i in 1:n.nodes) {
        readNormConstant <-
          scan(
            file = paste(fNormConstant, i, sep = ""),
            n = 1,
            quiet = T
          )
        readNormConstantPrior <-
          scan(file = paste(fNormConstantPrior, i, sep = ""),
               quiet = T)
        betahat <-
          betahat + read.table(file = paste(fBetahat, i, sep = ""),
                               colClasses = "numeric")[[1]] * readNormConstant * readNormConstantPrior
        ac <- ac + readNormConstant
      }
      betahat <- betahat / (D * E)

      write.table(
        file = fMostProbModels,
        thisMostProbModels,
        row.names = F,
        col.names = F
      )
      write.table(
        file = fPostProb,
        thisUnnorPostProb / (D * E),
        row.names = F,
        col.names = F
      )
      write.table(
        file = fInclusionProb,
        accum.InclusionProb,
        row.names = F,
        col.names = F
      )
      write.table(
        file = fProbDimension,
        accum.ProbDimension,
        row.names = F,
        col.names = F
      )
      write.table(file = fNormConstant,
                  D,
                  row.names = F,
                  col.names = F)
      write.table(file = fNormConstant,
                  D,
                  row.names = F,
                  col.names = F)
      write.table(file = fNormConstantPrior,
                  E,
                  row.names = F,
                  col.names = F)
      write.table(file = fBetahat,
                  betahat,
                  row.names = F,
                  col.names = F)
      write.table(
        file = fJointInclusionProb,
        accum.JointInclusionProb,
        row.names = F,
        col.names = F
      )

      ##############End of put together the results


    }
    else {

      Cresult <- myfun(c(1, (2^p - 1)), method = method)
      time <- Cresult[[8]]

    }


    #a function to transform the number of the model into a binary number.
    integer.base.b_C <- function(x, k) {
      #x is the number we want to express in binary
      #k is the number positions we need
      if (x == 0)
        return(rep(0, k))
      else{
        ndigits <- (floor(logb(x, base = 2)) + 1)
        res <- rep(0, ndigits)
        for (i in 1:ndigits) {
          #i <- 1
          res[i] <- (x %% 2)
          x <- (x %/% 2)
        }
        return(c(res, rep(0, k - ndigits)))
      }
    }

    #read the files given by C
    models <-
      as.vector(t(read.table(
        paste(wd, "/MostProbModels", sep = ""), colClasses = "numeric"
      )))
    prob <-
      as.vector(t(read.table(
        paste(wd, "/PostProb", sep = ""), colClasses = "numeric"
      )))
    incl <-
      as.vector(t(read.table(
        paste(wd, "/InclusionProb", sep = ""), colClasses = "numeric"
      )))
    joint <-
      as.matrix(read.table(paste(wd, "/JointInclusionProb", sep = ""), colClasses =
                             "numeric"))
    dimen <-
      as.vector(t(read.table(
        paste(wd, "/ProbDimension", sep = ""), colClasses = "numeric"
      )))
    betahat <-
      as.vector(t(read.table(
        paste(wd, "/betahat", sep = ""), colClasses = "numeric"
      )))


    #data.frame with Most probable models
    mod.mat <- as.data.frame(cbind(t(rep(0, (
      p + 1
    )))))

    names(mod.mat) <- c(namesx, "prob")

    N <- n.keep

    for (i in 1:N) {
      mod.mat[i, 1:p] <- integer.base.b_C(models[i], p)
      varnames.aux <- rep("", p)
      varnames.aux[mod.mat[i, 1:p] == 1] <- "*"
      mod.mat[i, 1:p] <- varnames.aux
    }

    mod.mat[, (p + 1)] <- prob[]

    inclusion <- incl #inclusion probabilities except for the intercept

    #the final result

    result <- list()

    result$time <- time #The time it took the programm to finish
    result$lmfull <- lmfull # The lm object for the full model
    if (!is.null(fixed.cov)) {
      result$lmnull <- lmnull # The lm object for the null model
    }

    result$variables <- namesx #The name of the competing variables
    result$n <- n #number of observations
    result$p <- p #number of competing variables
    result$k <- knull
    result$HPMbin <-
      integer.base.b_C(models[1], (p)) #The binary code for the HPM model
    names(result$HPMbin) <- namesx
    result$modelsprob <-
      mod.mat #A table with the n.keep most probable models ands its probability
    result$inclprob <-
      inclusion #inclusion probability for each variable
    names(result$inclprob) <- namesx

    result$jointinclprob <-
      data.frame(joint[1:p, 1:p], row.names = namesx)#data.frame for the joint inclusion probabilities
    names(result$jointinclprob) <- namesx
    #
    result$postprobdim <-
      dimen #vector with the dimension probabilities.
    names(result$postprobdim) <-
      (0:p) + knull #dimension of the true model
    #
		
		result$C<- scan(file=paste(wd,"/NormConstant", sep=""), quiet = T)
		
    #result$betahat <- betahat
    #rownames(result$betahat)<-namesx
    #names(result$betahat) <- "BetaHat"
    result$call <- match.call()
    if (!parallel){
      result$method <- "full"
    }
    else result$method <- "parallel"
		
	result$prior.betas<- prior.betas	
	result$prior.models<- prior.models
	result$priorprobs<- priorprobs
		

    class(result) <- "Bvs"
    result

  }
#' Print an object of class \code{Bvs}
#'
#' Print an object of class \code{Bvs}. The ten most probable models (among the visited ones if the object was created with
#' GibbsBvs) are shown jointly with their Bayes factors and an estimation of their posterior probability based on the estimation
#' of the normalizing constant.
#'
#' @export
#' @param x An object of class \code{Bvs}
#' @param ... Additional parameters to be passed
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#'   Maintainer: <anabel.forte@@uv.es>
#' @seealso See \code{\link[BayesVarSel]{Bvs}},
#'   \code{\link[BayesVarSel]{GibbsBvs}} for creating objects of the class
#'   \code{Bvs}.
#' @examples
#'
#' \dontrun{
#' #Analysis of Crime Data
#' #load data
#' data(UScrime)
#'
#' #Default arguments are Robust prior for the regression parameters
#' #and constant prior over the model space
#' #Here we keep the 1000 most probable models a posteriori:
#' crime.Bvs<- Bvs(formula= y ~ ., data=UScrime, n.keep=1000)
#'
#' #A look at the results:
#' print(crime.Bvs)
#' }
#'
print.Bvs <-
  function(x, ...){
    cat("\n")
    cat("Call:\n")
    print(x$call)

    if(x$method=="gibbs"){
       p <- x$p
       n.iter <- dim(x$modelslogBF)[1]
       dimmodels<- rowSums(x$modelslogBF[,1:p])+1
       logpostprob<- x$modelslogBF[, (p+1)] + log(x$priorprobs[dimmodels])-log(x$C)-log(sum(x$priorprobs*choose(p,0:p)))
 
       ordenado <- x$modelslogBF[order(logpostprob,decreasing = T),]
       ordenado<- cbind(x$modelslogBF, logpostprob)
       ordenado <- ordenado[order(logpostprob,decreasing = T),]
       models<- ordenado[!duplicated(ordenado),]
 
			ten<- min(10, dim(models)[1])
      #mod.mat <- as.data.frame(models[1:ten,1:p])
      mod.mat <- as.data.frame(models[1:ten,])

      for (i in 1:ten) {
        varnames.aux <- rep("", p)
        varnames.aux[models[i,1:p] == 1] <- "*"
        mod.mat[i,1:p] <- varnames.aux
      }
      cat("\nThe ", ten, " most probable models among the visited ones are:\n")
      print(mod.mat)
			cat("---\n")
			cat("Code: Column logBFi0 is the log of Bayes factor and\n")
			cat("column logpostprob is an estimation of posterior probabilities (log-scale)\n")
			cat("based on the normalizing constant.")

      }

	  if(x$method=="gibbsWithFactors"){
			summary(x)
		}

    if(x$method=="full" | x$method=="parallel"){
	    n.keep<-dim(x$modelsprob)[1]
      if(n.keep<=10){
        cat(paste("\nThe",n.keep,"most probable models and their probabilities are:\n",sep=" "))
        print(x$modelsprob)
      }else{
        cat("\nThe 10 most probable models and their probabilities are:\n")
        print(x$modelsprob[1:10, ])
        cat("\n(The remanining", n.keep - 10, "models are kept but omitted in this print)")
      }
    }
    cat("\n\n")

  }



#' Summary of an object of class \code{Bvs}
#'
#' Summary of an object of class \code{Bvs}, providing inclusion probabilities and a representation of
#' the Median Probability Model and the Highest Posterior probability Model.
#'
#' @export
#' @param object An object of class \code{Bvs}
#' @param ... Additional parameters to be passed
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#'   Maintainer: <anabel.forte@@uv.es>
#' @seealso See \code{\link[BayesVarSel]{Bvs}},
#'   \code{\link[BayesVarSel]{GibbsBvs}} for creating objects of the class
#'   \code{Bvs}.
#' @examples
#'
#' \dontrun{
#' #Analysis of Crime Data
#' #load data
#' data(UScrime)
#'
#' #Default arguments are Robust prior for the regression parameters
#' #and constant prior over the model space
#' #Here we keep the 1000 most probable models a posteriori:
#' crime.Bvs<- Bvs(formula= y ~ ., data=UScrime, n.keep=1000)
#'
#' #A look at the results:
#' summary(crime.Bvs)
#' }
#'
summary.Bvs <-
  function(object,...){

    #we use object because it is requiered by S3 methods
    z <- object
    p <- z$p
    if (!inherits(object, "Bvs"))
      warning("calling summary.Bvs(<fake-Bvs-x>) ...")
    ans <- list()
    #ans$coefficients <- z$betahat
    #dimnames(ans$coefficients) <- list(names(z$lm$coefficients),"Estimate")

    HPM <- z$HPMbin
    MPM <- as.numeric(z$inclprob >= 0.5)
    astHPM <- matrix(" ", ncol = 1, nrow = p)
    astMPM <- matrix(" ", ncol = 1, nrow = p)
    astHPM[HPM == 1] <- "*"
    astMPM[MPM == 1] <- "*"

    incl.prob <- z$inclprob

    summ.Bvs <- as.data.frame(cbind(round(incl.prob ,digits = 4), astHPM, astMPM))
    dimnames(summ.Bvs) <- list(z$variables, c("Incl.prob.", "HPM", "MPM"))

    ans$summary <- summ.Bvs
    ans$method <- z$method
    ans$call <- z$call

    cat("\n")
    #cat("Call:\n")
    #print(ans$call)
    #cat("\n")

    cat("Inclusion Probabilities:\n")
    print(ans$summary)
    cat("---\n")
    cat("Code: HPM stands for Highest posterior Probability Model and\n MPM for Median Probability Model.\n ")
    if (ans$method == "gibbs") {
      cat("Results are estimates based on the visited models.\n")
    }
    class(ans) <- "summary.Bvs"
    return(invisible(ans))
  }

#
check.prior.betas <- function(prior.betas){
	#checking that the prior for betas is implemented (and is given the right label)
	#and returning the label for its internal representation
	success <- 0
	if (prior.betas == "Robust") {success <- 1; prior.betas.internal <- "r"}
	if (prior.betas == "gZellner") {success <- 1; prior.betas.internal <- "g"}
	if (prior.betas == "Liangetat") {success <- 1; prior.betas.internal <- "l"}
	if (prior.betas == "ZellnerSiow") {success <- 1; prior.betas.internal <- "z"}
	if (prior.betas == "FLS") {success <- 1; prior.betas.internal <- "f"}
	if (prior.betas == "intrinsic.MGC") {success <- 1; prior.betas.internal <- "i"}
	if (prior.betas == "Robust.G") {success <- 1; prior.betas.internal <- "w"}
	if (prior.betas == "IHG") {success <- 1; prior.betas.internal <- "x"}
		
	if (success != 1) stop("Option for argument prior.betas not recognized (use literally the label as defined in the help)")	

		return (prior.betas.internal)

}

check.prior.modelspace <- function(prior.models, priorprobs, p, wd){
	#checking that the prior for model space is implemented (and is given the right label)
	#and returning a label for its internal representation
	success <- 0
	if (prior.models == "Constant") {success <- 1; prior.models.internal <- "c"}
	if (prior.models == "ScottBerger") {success <- 1; prior.models.internal <- "s"}
	if (prior.models == "User") {
		success <- 1; prior.models.internal <- "u"
		if (is.null(priorprobs)) stop("A valid vector of prior probabilities must be provided\n")
		if (length(priorprobs) != (p + 1)) stop("Vector of prior probabilities with incorrect length\n")
		if (sum(priorprobs < 0) > 0) stop("Prior probabilities must be positive\n")
		if (priorprobs[1] == 0) stop("Vector of prior probabilities not valid: All the theory here implemented works with the implicit assumption that the null model could be the true model\n")
	}
	
	if (success != 1) stop("Option for argument prior.models not recognized (use literally the label as defined in the help)")	

    #The zero here added always is for C compatibility
      write(
        priorprobs,
        ncolumns = 1,
        file = paste(wd, "/priorprobs.txt", sep = "")
      )

	return (prior.models.internal)
	
	
}