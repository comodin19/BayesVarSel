#' Bayesian Variable Selection with Factors for linear regression models using Gibbs
#' sampling.
#'
#' Numerical and factor variable selection from a Bayesian perspective. The posterior distribution is approximated
#' with Gibbs sampling
#'
#' In practical terms, \code{GibbsBvsF} can be understood as a version of \code{\link[BayesVarSel]{GibbsBvs}} in the presence of factors.
#' The methodology implemented in \code{GibbsBvsF} to handle variable selection problems with factors
#' has been proposed in Garcia-Donato and Paulo (2018) leading to a method
#' for which results  do not depend on how the factors are
#' coded (eg. via \code{\link[stats]{contrast}}).
#'
#' Internally, a rank defficient representation
#' of factors using dummies is used and the number of competing models considered is
#'
#' 2^(pnum+sum_j L_j),
#'
#' where pnum is the number of numerical variables and L_j is the number of levels in factor j.
#'
#' A main difference with \code{Bvs} and \code{GibbsBvs} (due to the presence of factors) concerns the prior
#' probabilities on the model space:
#'
#' The options \code{prior.models="SBSB"}, \code{prior.models="ConstConst"} and \code{prior.models="SBConst"}
#' acknowledge the "grouped" nature of the dummy variables representing
#' factors through the use of two stage
#' priors described in Garcia-Donato and Paulo (2018). In the first stage probabilities over factors and numerical
#' variables are specified and (conditional on these) within the second stage
#' the probablities are apportioned over the different submodels defined
#' by the dummies. The default option is "SBSB" which uses in both stages an assignment
#' of the type Scott-Berger so inversely proportional to the number of models of the same dimension. The
#' option "ConstConst" implements a uniform prior for both stages while "SBConst" uses a Scott-Berger prior
#' in the first stage and it is uniform in the second stage. Within all these priors, the prior inclusion probabilities
#' of factors and numerical variables are 1/2.
#'
#' The options \code{prior.models="Const"} and
#' \code{prior.models="SB"} do not have a staged structure and "Const" apportions the prior probabilities
#' uniformly over all possible models (2^(pnum+sum_j L_j)) and in "SB" the probability
#' is inversely proportional to the number of any model of the same dimension. In these cases, prior inclusion probabilities
#' of factors and numerical variables depend on the number of levels of factors and, in general, are not 1/2.
#' @export
#' @param formula Formula defining the most complex linear model in the
#' analysis. See details.
#' @param data data frame containing the data.
#' @param null.model A formula defining which is the simplest (null) model.
#' It should be nested in the full model. It is compulsory that the null model
#' contains the intercept and by default, the null model is defined
#' to be the one with just the intercept
#' @param prior.betas Prior distribution for regression parameters within each
#' model. Possible choices include "Robust", "Liangetal", "gZellner",
#' and "ZellnerSiow" (see details in \code{\link[BayesVarSel]{Bvs}}).
#' @param prior.models Prior distribution over the model space. Possible
#' choices (see details) are "Const", "SB", "ConstConst", "SBConst" and "SBSB" (the default).
#' @param n.iter The total number of iterations performed after the burn in
#' process.
#' @param init.model The model at which the simulation process starts. Options
#' include "Null" (the model only with the covariates specified in
#' \code{fixed.cov}), "Full" (the model defined by \code{formula}), "Random" (a
#' randomly selected model) and a vector with (pnum+sum_j L_j) zeros and ones defining a model.
#' @param n.burnin Length of burn in, i.e. number of iterations to discard at
#' the beginning.
#' @param n.thin Thinning rate. Must be a positive integer.  Set 'n.thin' > 1
#' to save memory and computation time if 'n.iter' is large. Default is 1. This
#' parameter jointly with \code{n.iter} sets the number of simulations kept and
#' used to construct the estimates so is important to keep in mind that a large
#' value for 'n.thin' can reduce the precision of the results
#' @param time.test If TRUE and the number of variables is large (>=21) a
#' preliminary test to estimate computational time is performed.
#' @param seed A seed to initialize the random number generator
#' @return \code{GibbsBvsF} returns an object of class \code{Bvs} with the
#' following elements: \item{time }{The internal time consumed in solving the
#' problem} \item{lmfull }{The \code{lm} class object that results when the
#' model defined by \code{formula} is fitted by \code{lm}}
#' \item{lmnull }{The
#' \code{lm } class object that results when the model defined by
#' \code{fixed.cov} is fitted by \code{lm}}
#' \item{variables }{The name of all the potential explanatory variables (numerical or factors)}
#' \item{n }{Number of observations}
#' \item{p }{Number of explanatory variables (both numerical and factors) to select from}
#' \item{k }{Number of fixed variables}
#' \item{HPMbin }{The binary expression of the most
#' probable model found.}
#' \item{inclprob }{A named vector with the
#' estimates of the inclusion probabilities of all the variables.}
#' \item{jointinclprob }{A \code{data.frame} with the estimates of the joint
#' inclusion probabilities of all the variables.}
#' \item{postprobdim }{Estimates
#' of posterior probabilities of the number of active variables in the true model (hence ranking from
#' \code{k } to \code{k+p}).}
#' \item{modelslogBF }{A matrix with both the binary representation of the
#' active variables in the MCMC after the burning period and the Bayes factor (log scale) of
#' that model to the null model.}
#' \item{modelswllogBF }{A matrix with both the binary representation of the
#' active variables (at the level of the levels in the factors) in the MCMC after the burning period and the Bayes factor (log scale) of
#' that model to the null model.}
#' \item{call }{The \code{call} to the
#' function.}
#' \item{C }{An estimation of the normalizing constant (C=sum Bi Pr(Mi), for Mi in the model space) using the method in George and McCulloch (1997).}
#' \item{positions }{A binary matrix with \code{p} rows and (pnum+sum_j L_j) columns. The 1's identify, for each variable (row) the position (column)
#' of dummies (in case of factor) or of the numerical variable grouped on that variable. (Its use is conceived for internal purposes).}
#' \item{positionsx }{A \code{p} dimensional binary vector, stating which of the competing variables is a numerical variable. (Its use is conceived for internal purposes).}
#' \item{method }{\code{gibbsWithFactors}}
#' @author Gonzalo Garcia-Donato and Anabel Forte
#' @seealso \code{\link[BayesVarSel]{plot.Bvs}} for several plots of the result.
#'
#' Under construction: \code{\link[BayesVarSel]{BMAcoeff}} for obtaining model averaged simulations
#' of regression coefficients and \code{\link[BayesVarSel]{predict.Bvs}} for
#' predictions.
#'
#' See \code{\link[BayesVarSel]{GibbsBvs}} and \code{\link[BayesVarSel]{Bvs}} when no factors are involved.
#' @references Garcia-Donato, G. and Martinez-Beneito, M.A.
#' (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian
#' variable selection problems with large model spaces. Journal of the American
#' Statistical Association, 108: 340-352.
#'
#' Garcia-Donato, G. and Paulo, R. (2018) Including factors in Bayesian variable selection
#' problems. arXiv:1709.07238.
#'
#' George E. and McCulloch R. (1997) Approaches for Bayesian variable
#' selection. Statistica Sinica, 7, 339:372.
#' @keywords package
#' @examples
#'
#' \dontrun{
#' data(diabetes, package="faraway")
#'
#' #remove NA's and the column with the id of samples:
#' diabetes2<- na.omit(diabetes)[,-1]
#'
#' #For reproducibility:
#' set.seed(16091956)
#' #Now run the main instruction
#' diabetesVS<- GibbsBvsF(formula= glyhb ~ ., data=diabetes2, n.iter=100000, n.burnin=5000)
#'
#' summary(diabetesVS)
#'
#' #A plot of the dimension of the true model,
#' plot(diabetesVS, option="dimension")
#'
#' #A joint inclusion plot
#' plot(diabetesVS, option="joint")
#'
#' #Now a similar exercise but with fixed variables:
#' diabetesVS2<- GibbsBvsF(formula= glyhb ~ ., null.model= glyhb ~ chol+stab.glu,
#' 		                   data=diabetes2, n.iter=100000, n.burnin=5000)
#'
#'
#' #and with fixed factors:
#' diabetesVS3<- GibbsBvsF(formula= glyhb ~ ., null.model= glyhb ~ chol+stab.glu+location,
#' 		                   data=diabetes2, n.iter=100000, n.burnin=5000)
#'
#'
#' }
#'
GibbsBvsF <-
  function(formula,
           data,
           null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
           prior.betas = "Robust",
           prior.models = "SBSB",
           n.iter = 10000,
           init.model = "Full",
           n.burnin = 500,
           n.thin = 1,
           time.test = TRUE,
					 seed = runif(1, 0, 16091956)) {

    formula <- as.formula(formula)

    null.model<- as.formula(null.model)

    #The response in the null model and in the full model must coincide
    if (formula[[2]] != null.model[[2]]){
      stop("The response in the full and null model does not coincide.\n")
    }

    #Let's define the result
    result <- list()

    #Get a tempdir as working directory
    wd <- tempdir()
    #remove all the previous documents in the working directory
    unlink(paste(wd, "*", sep = "/"))

    #evaluate the null model:
    lmnull <- lm(formula = null.model, data, y = TRUE, x = TRUE)
    fixed.cov <- dimnames(lmnull$x)[[2]]

    #Set the design matrix if fixed covariates present:
      #Eval the full model
      lmfull = lm(formula,
                  data = data,
                  y = TRUE,
                  x = TRUE)

			#Variables or factors that are a linear function of the others:
			if (lmfull$rank!=dim(lmfull$x)[2])
				stop("Some of the explanatory variables are a linear function of others\n")
			#Factors:
			#before:X.full <- lmfull$x
    	X.full<- get_rdX(lmfull) #rank defficient paramet

      namesx <- dimnames(X.full)[[2]]

      #check if null model is contained in the full one:
      namesnull <- dimnames(lmnull$x)[[2]]
      "%notin%" <- function(x, table) match(x, table, nomatch = 0) == 0

			#check that the null model contains the intercept:
			if (is.null(namesnull))
				stop("The null model should contain the intercept\n")

			if (sum(namesnull=="Intercept")==0 & sum(namesnull=="(Intercept)")==0)
				stop("The null model should contain the intercept\n")

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

      #the response variable for the C code
      Y <- lmnull$residuals

      #Design matrix of the null model
      X0 <- lmnull$x
      P0 <-
        X0 %*% (solve(t(X0) %*% X0)) %*% t(X0)#Intentar mejorar aprovechando lmnull
      knull <- dim(X0)[2]

      #matrix containing the covariates from which we want to select
			#Factors:
      X1<- X.full[, -fixed.pos] #before:X1 <- lmfull$x[, -fixed.pos]

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

      p <- dim(X)[2]#Number of covariates to select from


		#Factors:
		#positions is a matrix with number of rows equal to the number of regressors
		#(either factor or numeric) and number of columns the number of columns of X
		#Each row describes the position (0-1) in X of a regressor (several positions in case
		#this regressor is a factor)
		#Works with only intercept, but not with other null models: depvars<- attr(lmfull$terms, "term.labels")
		depvars<- setdiff(attr(lmfull$terms, "term.labels"), attr(lmnull$terms, "term.labels"))

    #Check if, among the competing variables, there are factors
		if (sum(attr(lmfull$terms, "dataClasses")[depvars]=="factor")==0){
			stop("No Factors found among the competing variables. Use GibbsBvs() instead\n")
		}

		positions<- matrix(0, ncol=p, nrow=length(depvars))
		for (i in 1:length(depvars)){positions[i,]<- grepl(depvars[i], colnames(X), fixed=T)}
		#positionsX is a vector of the same length as columns has X

		#with 1 in the position with a numeric variable:
		positionsx<- as.numeric(colSums(positions%*%t(positions))==1)

		write(positionsx, ncolumns=1, file = paste(wd, "/positionsx.txt", sep = ""))
    write(t(positions),
          ncolumns = p,
          file = paste(wd, "/positions.txt", sep = ""))
	  #both files are used to obtain prior probabilities and rank of matrices
		rownames(positions)<- depvars

    #write the data files in the working directory
    write(Y,
          ncolumns = 1,
          file = paste(wd, "/Dependent.txt", sep = ""))
    write(t(X),
          ncolumns = p,
          file = paste(wd, "/Design.txt", sep = ""))

    #The initial model:
    if (is.character(init.model) == TRUE) {
      im <- substr(tolower(init.model), 1, 1)
      if (im != "n" &&
          im != "f" && im != "r") {
        stop("Initial model not valid\n")
      }
      if (im == "n") {
        init.model <- rep(0, p)
      }
      if (im == "f") {
        init.model <- rep(1, p)
      }
      if (im == "r") {
        init.model <- rbinom(n = p,
                             size = 1,
                             prob = .5)
      }
    }
    else{
      init.model <- as.numeric(init.model > 0)
      if (length(init.model) != p) {
        stop("Initial model with incorrect length\n")
      }
    }

    write(
      init.model,
      ncolumns = 1,
      file = paste(wd, "/initialmodel.txt", sep = "")
    )

    #Info:
    cat("Info. . . .\n")
    cat("Most complex model has", dim(positions)[1] + knull, "numerical covariates and factors\n")
    if (!is.null(fixed.cov)) {
      if (knull > 1) {
        cat("From those",
            knull,
            "are fixed and we should select from the remaining",
            dim(positions)[1],
            "\n")
      }
      if (knull == 1) {
        cat("From those",
            knull,
            "is fixed (the intercept) and we should select from the remaining",
            dim(positions)[1],
            "\n")
      }
			cat(" Covariates (numerical variables):\n", depvars[positionsx==1], "\n",
			    "Factors:\n", depvars[positionsx==0], "\n\n")

	    }
    cat("The problem has a total of", 2 ^ (p), "competing models\n")
    iter <- n.iter
    cat("Of these,", n.burnin + n.iter, "are sampled with replacement\n")

    cat("Then,",
        floor(iter / n.thin),
        "are kept and used to construct the summaries\n")


    #Note: priorprobs.txt is a file that is needed only by the "User" routine. Nevertheless, in order
    #to mantain a common unified version the source files of other routines also reads this file
    #although they do not use. Because of this we create this file anyway.
      priorprobs <- rep(0, p + 1)
      write(
        priorprobs,
        ncolumns = 1,
        file = paste(wd, "/priorprobs.txt", sep = "")
      )


		#Factors:
		#here the added index "2" makes reference of the hierarchical corresponding prior but only keeping
		#a model of the same class (copies are removed and only the full within each class is kept)
		if (prior.models!="SBSB" & prior.models!="ConstConst" & prior.models!="SB" & prior.models!="Const" & prior.models!="SBConst")
			{stop("Prior over the model space not supported\n")}


    if (prior.betas == "Unitary"){write(0, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}
    if (prior.betas == "Robust"){write(1, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}
    if (prior.betas == "Liangetal"){write(4, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}
    if (prior.betas == "gZellner"){write(2, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}
    if (prior.betas == "ZellnerSiow"){write(5, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}
    if (prior.betas == "Robust2"){write(6, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}

	  if (prior.betas == "FLS"){stop("Prior FLS not yet supported\n")}

		if (prior.betas != "Unitary" & prior.betas != "Robust" & prior.betas != "Liangetal" &
			    prior.betas != "gZellner" & prior.betas != "ZellnerSiow" & prior.betas != "FLS" &
					prior.betas != "Robust2") {stop("Dont recognize the prior for betas\n")}

		if (prior.models=="SBSB"){method<- "rSBSB"}
		if (prior.models=="ConstConst"){method<- "rConstConst"}
		if (prior.models=="SBConst"){method<- "rSBConst"}
		if (prior.models=="SB"){method<- "rSB"}
		if (prior.models=="Const"){method<- "rConst"}

    estim.time <- 0

    #Call the corresponding function:
    result <- switch(
      method,
      "rSBSB" = .C(
        "GibbsFSBSB",
        as.character(""),
        as.integer(n),
        as.integer(p),
        as.integer(floor(n.iter / n.thin)),
        as.character(wd),
        as.integer(n.burnin),
        as.double(estim.time),
        as.integer(knull),
        as.integer(n.thin),
        as.integer(seed)
      ),
      "rConstConst" = .C(
        "GibbsFConstConst",
        as.character(""),
        as.integer(n),
        as.integer(p),
        as.integer(floor(n.iter / n.thin)),
        as.character(wd),
        as.integer(n.burnin),
        as.double(estim.time),
        as.integer(knull),
        as.integer(n.thin),
        as.integer(seed)
      ),
      "rSBConst" = .C(
        "GibbsFSBConst",
        as.character(""),
        as.integer(n),
        as.integer(p),
        as.integer(floor(n.iter / n.thin)),
        as.character(wd),
        as.integer(n.burnin),
        as.double(estim.time),
        as.integer(knull),
        as.integer(n.thin),
        as.integer(seed)
      ),
      "rSB" = .C(
        "GibbsFSB",
        as.character(""),
        as.integer(n),
        as.integer(p),
        as.integer(floor(n.iter / n.thin)),
        as.character(wd),
        as.integer(n.burnin),
        as.double(estim.time),
        as.integer(knull),
        as.integer(n.thin),
        as.integer(seed)
      ),
      "rConst" = .C(
        "GibbsFConst",
        as.character(""),
        as.integer(n),
        as.integer(p),
        as.integer(floor(n.iter / n.thin)),
        as.character(wd),
        as.integer(n.burnin),
        as.double(estim.time),
        as.integer(knull),
        as.integer(n.thin),
        as.integer(seed)
      ))

    time <- result[[7]]


    #read the files given by C
    models <- as.vector(t(read.table(paste(wd,"/MostProbModels",sep=""),colClasses="numeric")))
    incl <- as.vector(t(read.table(paste(wd,"/InclusionProb",sep=""),colClasses="numeric")))
    joint <- as.matrix(read.table(paste(wd,"/JointInclusionProb",sep=""),colClasses="numeric"))
    dimen <- as.vector(t(read.table(paste(wd,"/ProbDimension",sep=""),colClasses="numeric")))
    betahat<- as.vector(t(read.table(paste(wd,"/betahat",sep=""),colClasses="numeric")))
    allmodels<- as.matrix(read.table(paste(wd,"/AllModels",sep=""),colClasses="numeric"))
    allBF<- as.vector(t(read.table(paste(wd,"/AllBF",sep=""),colClasses="numeric")))

    #Log(BF) for every model
    modelslBF<- cbind(allmodels, log(allBF))
    colnames(modelslBF)<- c(namesx, "logBFi0")

		############
		#Specific to factors:
		#
		#Recall now that the number of vars (either numerical or factors) is dim(positions)[1]

		#Resampling removing the saturated models (keeping the oversaturated)
		if (prior.models == "SBSB"){modelslBFwR<- resamplingSBSB(modelslBF, positions)}
		if (prior.models == "ConstConst"){modelslBFwR<- resamplingConstConst(modelslBF, positions)}
		if (prior.models == "SB"){modelslBFwR<- resamplingSB(modelslBF, positions)}
		if (prior.models == "Const"){modelslBFwR<- resamplingConst(modelslBF, positions)}
		if (prior.models == "SBConst"){modelslBFwR<- resamplingSBConst(modelslBF, positions)}

		#Now we convert the sampled models to vars and factors:
		cat(dim(positions),"\n")
		if (dim(positions)[1] > 1)
			modelslBFF<- t(apply(modelslBFwR[,-(p+1)], MARGIN=1, FUN=function(x,M){as.numeric(x%*%M>0)}, M=t(positions)))
		else
			modelslBFF<- matrix(as.numeric(modelslBFwR[,-(p+1)]%*%t(positions)>0), ncol=1)

        #Inclusion probabilities:

		inclusion <- colMeans(modelslBFF)
		names(inclusion)<- depvars

		#HPM with factors:
		HPMFbin<- models%*%t(positions)

		#joint inclusion probs:
		jointinclprob<- matrix(0, ncol=dim(positions)[1], nrow=dim(positions)[1])

		for (i in 1:dim(modelslBFF)[1]){
			jointinclprob<- jointinclprob + matrix(modelslBFF[i,], ncol=1)%*%matrix(modelslBFF[i,], nrow=1)
		}
			jointinclprob<- jointinclprob/dim(modelslBFF)[1]
		colnames(jointinclprob)<- depvars; rownames(jointinclprob)<- depvars

		#Dimension of the true model:
		dimenF<- c(rowSums(modelslBFF), 0:dim(positions)[1])
		dimenF<- (table(dimenF)-1)/dim(modelslBFF)[1]
		names(dimenF)<- (0:dim(positions)[1])+knull

		#Attach the column with the log(BF)
		modelslBFF<- cbind(modelslBFF, modelslBFwR[,"logBFi0"])
		colnames(modelslBFF)<- c(depvars, "logBFi0")

    result <- list()
    #
    result$time <- time #The time it took the programm to finish
    result$lmfull <- lmfull # The lm object for the full model
    if(!is.null(fixed.cov)){
      result$lmnull <- lmnull # The lm object for the null model
    }

    result$variables <- depvars #The name of the competing variables
    result$n <- n #number of observations
    result$p <- length(depvars) #number of competing variables
    result$k <- knull#number of fixed covariates
    result$HPMbin <- (HPMFbin)#The binary code for the HPM model

    result$modelslogBF <- modelslBFF#The binary code for all the visited models (after n.thin is applied) and the correspondent log(BF)

		#Keep the visited models at the level of levels for posterior analyses
		result$modelswllogBF<- modelslBFwR

    result$inclprob <- inclusion #inclusion probability for each variable

    result$jointinclprob <- data.frame(jointinclprob) #data.frame for the joint inclusion probabilities
    #
    result$postprobdim <- dimenF #vector with the dimension probabilities.

		result$positions<- positions
		result$positionsx<- positionsx

		#
		#
		#####################

    result$call <- match.call()

    result$method <- "gibbsWithFactors"
    class(result)<- "Bvs"
    result


  }


