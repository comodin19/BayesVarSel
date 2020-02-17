#' Bayesian Model Averaged estimations of regression coefficients
#'
#' Samples of the model averaged objective posterior distribution of regression
#' coefficients
#'
#'
#' The distribution that is sampled from is the discrete mixture of the
#' (objective) posterior distributions of the regression coefficients with
#' weights proportional to the posterior probabilities of each model. That is,
#' from
#'
#' \eqn{latex}{ sum_M f(\beta | data, M) Pr(M | data)}
#'
#' The models used in the mixture above are the retained best models (see the
#' argument \code{n.keep} in \link[BayesVarSel]{Bvs}) if \code{x} was generated
#' with \code{Bvs} and the sampled models with the associated frequencies if
#' \code{x} was generated with \code{GibbsBvs}. The formula for the objective
#' posterior distribution within each model \eqn{latex}{f(\beta | data, M)} is
#' taken from Bernardo and Smith (1994) page 442.
#'
#' Note: The above mixture is potentially highly multimodal and this command
#' ends with a multiple plot with the densities of the different regression
#' coefficients to show the user this peculiarity. Hence which summaries should
#' be used to describe this distribution is a delicate issue and standard
#' functions like the mean and variance are not recommendable.
#'
#' @param x An object of class \code{Bvs}
#' @param n.sim Number of simulations to be produced
#' @param method Text specifying the matrix decomposition used to determine the
#' matrix root of 'sigma' when simulating from the multivariate t distribution.
#' Possible methods are eigenvalue decomposition ('"eigen"', default), singular
#' value decomposition ('"svd"'), and Cholesky decomposition ('"chol"'). See
#' the help of command \code{rmvnorm} in package \code{mvtnorm} for more
#' details
#' @return \code{BMAcoeff} returns an object of class \code{bma.coeffs} which
#' is a matrix with \code{n.sim} rows with the simulations. Each column of the
#' matrix corresponds to a regression coefficient in the full model.
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#' Maintainer: <anabel.forte@@uv.es>
#' @export
#' @seealso See \code{\link[BayesVarSel]{histBMA}} for a histogram-like
#' representation of the columns in the object. See
#' \code{\link[BayesVarSel]{Bvs}} and
#' \code{\link[BayesVarSel]{GibbsBvs}} for creating objects of the class
#' \code{Bvs}. See \code{\link[mvtnorm]{rmvnorm} for details about argument
#' method.}
#' @examples
#'
#' \dontrun{
#'
#' #Analysis of Crime Data
#' #load data
#' data(UScrime)
#'
#' crime.Bvs<- Bvs(formula= y ~ ., data=UScrime, n.keep=1000)
#' crime.Bvs.BMA<- BMAcoeff(crime.Bvs, n.sim=10000)
#' #the best 1000 models are used in the mixture
#'
#' #We could force all  possible models to be included in the mixture
#' crime.Bvs.all<- Bvs(formula= y ~ ., data=UScrime, n.keep=2^15)
#' crime.Bvs.BMA<- BMAcoeff(crime.Bvs.all, n.sim=10000)
#' #(much slower as this implies ordering many more models...)
#'
#' #With the Gibbs algorithms:
#' data(Ozone35)
#'
#' Oz35.GibbsBvs<- GibbsBvs(formula= y ~ ., data=Ozone35, prior.betas="gZellner",
#' prior.models="Constant", n.iter=10000, init.model="Full", n.burnin=100,
#' time.test = FALSE)
#' Oz35.GibbsBvs.BMA<- BMAcoeff(Oz35.GibbsBvs, n.sim=10000)
#'
#'
#' }
#'
BMAcoeff <- function(x, n.sim = 10000, method = "svd") {
  #x is an object of class Bvs

  #simulations of the model averaging mixture of the objective posterior distributions of the regression
  #coefficients. Weights are taken from the previous Bvs or GibbsBvs run.

  #(next is for compatibility between notations for Intercept in lm and Bvs)
  if (!is.null(x$lmnull)) {
    if (colnames(x$lmnull$x)[1] == "(Intercept)") {
      colnames(x$lmnull$x)[1] <- "Intercept"
    }
  }
  if (colnames(x$lmfull$x)[1] == "(Intercept)") {
    colnames(x$lmfull$x)[1] <- "Intercept"
  }
  #name of dep var
  name.y <- colnames(x$lmfull$model[1])

  #results are given as a matrix bma.coeffs
  bma.coeffs <-
    matrix(0, nrow = n.sim, ncol = length(x$lmfull$coefficients))
  colnames(bma.coeffs) <- colnames(x$lmfull$x)

  #differentiate if method="full" (enumeration) or method="Gibbs" or method="gibbsWithFactors"
  if (x$method == "full" | x$method == "parallel") {
    cat("\n")
    cat("Simulations obtained using the best",
        dim(x$modelsprob)[1],
        "models\n")
    cat("that accumulate",
        round(sum(x$modelsprob[, "prob"]), 2),
        "of the total posterior probability\n")

    #draw n.sim models with replacement and with weights proportional
    #to their posterior prob:
    models <-
      sample(
        x = dim(x$modelsprob)[1],
        size = n.sim,
        replace = TRUE,
        prob = x$modelsprob[, "prob"]
      )
    #table of these models
    t.models <- table(models)
    cs.tmodels <- cumsum(t.models)

    X <- x$lmfull$x

    for (iter in 1:length(t.models)) {
      #rMD is model drawn (a number between 1 and n.keep)
      rMD <-
        as.numeric(names(t.models)[iter])
      howmany <- t.models[iter]

      #covs in that model rMD (apart from the fixed ones)
      covsrMD <- names(x$modelsprob[rMD, ])[x$modelsprob[rMD, ] == "*"]
      #the data in model drawn (dependent variable in first column)
      datarMD <-
        as.data.frame(cbind(x$lmfull$model[, name.y], X[, covsrMD]))

      colnames(datarMD) <- c(name.y, covsrMD)
      #now add the fixed.cov if any
      if (!is.null(x$lmnull)) {
        datarMD <- cbind(datarMD, x$lmnull$x)
      }

      #remove rare characters because a double application of lm command
      #(happens for instance in covariates with ":")
      colnames(datarMD) <- gsub("`", "", colnames(datarMD))

      #formula for that model
      formMD <- as.formula(paste(name.y, "~.-1", sep = ""))

      #fit
      fitrMD <-
        lm(formula = formMD,
           data = as.data.frame(datarMD),
           qr = TRUE)

      #simulated value for the beta:
      #simple version
      #rcoeff<- rmvnorm(n=howmany, mean=fitrMD$coefficients, sigma=vcov(fitrMD))
      #exact version
      Rinv <- qr.solve(qr.R(fitrMD$qr))
      iXtX <- Rinv %*% t(Rinv)
      Sigma <- sum(fitrMD$residuals * datarMD[, name.y]) * iXtX / fitrMD$df
      rcoeff <-
        rmvt(
          n = howmany,
          sigma = Sigma,
          df = fitrMD$df,
          delta = fitrMD$coefficients,
          type = "shifted",
          method = method
        )
      #      if(sum(names(fitrMD$coefficients)%in%colnames(bma.coeffs))!= length(names(fitrMD$coefficients))){
      #        stop("The names of your covariates may contain a non-standar charater ($,',:,ect.). Please check.")
      #      }

      bma.coeffs[(max(cs.tmodels[iter - 1], 0) + 1):cs.tmodels[iter], names(fitrMD$coefficients)] <-
        rcoeff

    }

  }
  if (x$method == "gibbs") {
    cat("\n")
    cat("Simulations obtained using the ",
        dim(x$modelslogBF)[1],
        " sampled models.\n")
    cat("Their frequencies are taken as the true posterior probabilities\n")

    #draw n.sim models with replacement (weights are implicitly proportional
    #to their posterior prob since repetitions are included in the sample):
    models <-
      sample(x = dim(x$modelslogBF)[1],
             size = n.sim,
             replace = TRUE)
    #table of these models
    t.models <- table(models)
    cs.tmodels <- cumsum(t.models)

    X <- x$lmfull$x

    for (iter in 1:length(t.models)) {
      #rMD is model drawn (a number between 1 and n.keep)
      rMD <-
        as.numeric(names(t.models)[iter])
      howmany <- t.models[iter]

      #covs in that model rMD (apart from the fixed ones)
      covsrMD <- names(x$modelslogBF[rMD, ])[x$modelslogBF[rMD, ] == "1"]
      #the data in model drawn (dependent variable in first column)
      datarMD <-
        as.data.frame(cbind(x$lmfull$model[, name.y], X[, covsrMD]))

      colnames(datarMD) <- c(name.y, covsrMD)
      #now add the fixed.cov if any
      if (!is.null(x$lmnull)) {
        datarMD <- cbind(datarMD, x$lmnull$x)
      }

      #remove rare characters because a double application of lm command
      #(happens for instance in covariates with ":")
      colnames(datarMD) <- gsub("`", "", colnames(datarMD))
      #formula for that model
      formMD <- as.formula(paste(name.y, "~.-1", sep = ""))

      #fit
      fitrMD <-
        lm(formula = formMD,
           data = as.data.frame(datarMD),
           qr = TRUE)

      #simulated value for the beta:
      #simple version
      #rcoeff<- rmvnorm(n=howmany, mean=fitrMD$coefficients, sigma=vcov(fitrMD))
      #exact version
      Rinv <- qr.solve(qr.R(fitrMD$qr))
      iXtX <- Rinv %*% t(Rinv)
      Sigma <- sum(fitrMD$residuals * datarMD[, name.y]) * iXtX / fitrMD$df
      rcoeff <-
        rmvt(
          n = howmany,
          sigma = Sigma,
          df = fitrMD$df,
          delta = fitrMD$coefficients,
          type = "shifted",
          method = method
        )
      #      if(sum(names(fitrMD$coefficients)%in%colnames(bma.coeffs))!= length(names(fitrMD$coefficients))){
      #        stop("The names of your covariates may contain a non-standar charater ($,',:,ect.). Please check and replace.")
      #      }

      bma.coeffs[(max(cs.tmodels[iter - 1], 0) + 1):cs.tmodels[iter], names(fitrMD$coefficients)] <-
        rcoeff

    }
  }
  if (x$method == "gibbsWithFactors") {
    cat("\n")
    cat("Simulations obtained using the ",
        dim(x$modelslogBF)[1],
        " sampled models.\n")
    cat("Their frequencies are taken as the true posterior probabilities\n")
		
  	Xf<- get_rdX(x$lmfull)
		x$lmfull$coefficients<- rep(0, dim(Xf)[2])
		x$lmfull$x<- Xf
	  if (colnames(x$lmfull$x)[1] == "(Intercept)") {
	    colnames(x$lmfull$x)[1] <- "Intercept"
	  }
		
	  bma.coeffs <- matrix(0, nrow = n.sim, ncol = length(x$lmfull$coefficients))
	  colnames(bma.coeffs) <- colnames(x$lmfull$x)
    #draw n.sim models with replacement (weights are implicitly proportional
    #to their posterior prob since repetitions are included in the sample):
    models <-
      sample(x = dim(x$modelslogBF)[1],
             size = n.sim,
             replace = TRUE)
    #table of these models
    t.models <- table(models)
    cs.tmodels <- cumsum(t.models)

    X <- x$lmfull$x

    for (iter in 1:length(t.models)) {
			cat(iter,"\n")
      #rMD is model drawn (a number between 1 and n.keep)
      rMD <- as.numeric(names(t.models)[iter])
      howmany <- t.models[iter]

      #covs in that model rMD (apart from the fixed ones)
      covsrMD <- names(x$modelswllogBF[rMD, ])[x$modelswllogBF[rMD, ] == "1"]
      #the data in model drawn (dependent variable in first column)
      datarMD <- as.data.frame(cbind(x$lmfull$model[, name.y], X[, covsrMD]))

      colnames(datarMD) <- c(name.y, covsrMD)
      #now add the fixed.cov if any
      if (!is.null(x$lmnull)) {
        datarMD <- cbind(datarMD, x$lmnull$x)
      }

      #remove rare characters because a double application of lm command
      #(happens for instance in covariates with ":")
      colnames(datarMD) <- gsub("`", "", colnames(datarMD))
      #formula for that model
      formMD <- as.formula(paste(name.y, "~.-1", sep = ""))

      #fit
      fitrMD <-
        lm(formula = formMD,
           data = as.data.frame(datarMD),
           qr = TRUE)

      #simulated value for the beta:
      #simple version
      #rcoeff<- rmvnorm(n=howmany, mean=fitrMD$coefficients, sigma=vcov(fitrMD))
      #exact version
      Rinv <- qr.solve(qr.R(fitrMD$qr))
      iXtX <- Rinv %*% t(Rinv)
      Sigma <- sum(fitrMD$residuals * datarMD[, name.y]) * iXtX / fitrMD$df
      rcoeff <-
        rmvt(
          n = howmany,
          sigma = Sigma,
          df = fitrMD$df,
          delta = fitrMD$coefficients,
          type = "shifted",
          method = method
        )
      #      if(sum(names(fitrMD$coefficients)%in%colnames(bma.coeffs))!= length(names(fitrMD$coefficients))){
      #        stop("The names of your covariates may contain a non-standar charater ($,',:,ect.). Please check and replace.")
      #      }

      bma.coeffs[(max(cs.tmodels[iter - 1], 0) + 1):cs.tmodels[iter], gsub("`", "", names(fitrMD$coefficients))] <-
        rcoeff

    }
	
  }
	
	

  class(bma.coeffs) <- "bma.coeffs"
    return(bma.coeffs)

}
