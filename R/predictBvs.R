#' Bayesian Model Averaged predictions
#'
#' Samples of the model averaged objective predictive distribution
#'
#'
#' The distribution that is sampled from is the discrete mixture of the
#' (objective) predictive distribution with weights proportional to the
#' posterior probabilities of each model. That is, from
#'
#' \eqn{latex}{ sum_M f(y^* | data, newdata, M) Pr(M | data)}
#'
#' The models used in the mixture above are the retained best models (see the
#' argument \code{n.keep} in \link[BayesVarSel]{Bvs}) if \code{x} was generated
#' with \code{Bvs} and the sampled models with the associated frequencies if
#' \code{x} was generated with \code{GibbsBvs}. The formula for the objective
#' predictive distribution within each model \eqn{latex}{f(\beta | data, M)} is
#' taken from Bernardo and Smith (1994) page 442.
#'
#' @export
#' @param x An object of class \code{Bvs}
#' @param newdata A data frame in which to look for variables with which to
#' predict
#' @param n.sim Number of simulations to be produced
#' @return \code{predictBvs} returns a matrix with \code{n.sim} rows with the
#' simulations. Each column of the matrix corresponds to each of the
#' configurations for the covariates defined in \code{newdata}.
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#' Maintainer: <anabel.forte@@uv.es>
#' @seealso See \code{\link[BayesVarSel]{Bvs}}, \code{\link[BayesVarSel]{PBvs}}
#' and \code{\link[BayesVarSel]{GibbsBvs}} for creating objects of the class
#' \code{Bvs}.
#' @references Bernardo, J. M. and Smith, A. F. M.
#' (1994)<DOI:10.1002/9780470316870> Bayesian Theory. Chichester: Wiley.
#' @examples
#'
#' \dontrun{
#'
#' #Analysis of Crime Data
#' #load data
#' data(UScrime)
#'
#' crime.Bvs<- Bvs(formula="y~.", data=UScrime, n.keep=1000)
#' #predict a future observation associated with the first two sets of covariates
#' crime.Bvs.predict<- predictBvs(crime.Bvs, newdata=UScrime[1:2,], n.sim=10000)
#' #(Notice the best 1000 models are used in the mixture)
#'
#' #Here you can use standard summaries to describe the underlying predictive distribution
#' #summary(crime.Bvs.predict)
#' #
#' #To study more in deep the first set:
#' #plot(density(crime.Bvs.predict[,1]))
#' #Point prediction
#' #median(crime.Bvs.predict[,1])
#' #A credible 95% interval for the prediction:
#' #lower bound:
#' #quantile(crime.Bvs.predict[,1], probs=0.025)
#' #upper bound:
#' #quantile(crime.Bvs.predict[,1], probs=0.975)
#'
#' }
#'
predictBvs <- function(x, newdata, n.sim = 10000) {
  #x is an object of class Bvs

  #simulations of the model averaging mixture of the objective posterior predictive distributions.
  #Weights are taken from the previous Bvs or GibbsBvs run.

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

  if (!is.data.frame(newdata)) {
    stop("newdata must be a data.frame\n")
  }

  #Now add the intercept if needed
  if (colnames(x$lmfull$x)[1] == "Intercept") {
    newdata$Intercept <- 1
  }

  #results are given as a matrix rpredictions
  rpredictions <- matrix(0, nrow = n.sim, ncol = dim(newdata)[1])

  #differentiate if method="full" (enumeration) or method="Gibbs"
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
      #formula for that model
      formMD <- as.formula(paste(name.y, "~.-1", sep = ""))

      #fit
      fitrMD <-
        lm(
          formula = formMD,
          data = as.data.frame(datarMD),
          qr = TRUE,
          x = TRUE
        )

      #corresponding new design for that model:
      newdataMD <- newdata[, names(fitrMD$coefficients)]

      fnx <- apply(
        X = as.matrix(newdataMD),
        MARGIN = 1,
        FUN = function(x, DM) {
          Xstar <- rbind(DM, x)
          qrXstar <- qr(Xstar)
          Rinv <- qr.solve(qr.R(qrXstar))
          iXstartXstar <- Rinv %*% t(Rinv)
          1 - t(x) %*% iXstartXstar %*% x
        },
        DM = as.matrix(fitrMD$x)
      )

      sigmas <- sum(fitrMD$residuals * datarMD[, name.y]) / (fnx * fitrMD$df)
      #compute means:
      means <- as.matrix(newdataMD) %*% fitrMD$coefficients

      #simulated value for the new y's:
      if (dim(newdata)[1] == 1) {
        sigmaM <- as.matrix(sigmas, nr = 1, nc = 1)
      }
      if (dim(newdata)[1] > 1) {
        sigmaM <- diag(sigmas)
      }
      rpredictions[(max(cs.tmodels[iter - 1], 0) + 1):cs.tmodels[iter], ] <-
        rmvt(
          n = howmany,
          delta = means,
          sigma = sigmaM,
          df = fitrMD$df,
          type = "shifted"
        )

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
      #formula for that model
      formMD <- as.formula(paste(name.y, "~.-1", sep = ""))

      #fit
      fitrMD <-
        lm(
          formula = formMD,
          data = as.data.frame(datarMD),
          qr = TRUE,
          x = TRUE
        )

      #corresponding new design for that model:
      newdataMD <- newdata[, names(fitrMD$coefficients)]

      #compute scales (one for each configuration)
      #(notation fnx from Bernardo's book)
      fnx <- apply(
        X = as.matrix(newdataMD),
        MARGIN = 1,
        FUN = function(x, DM) {
          Xstar <- rbind(DM, x)
          qrXstar <- qr(Xstar)
          Rinv <- qr.solve(qr.R(qrXstar))
          iXstartXstar <- Rinv %*% t(Rinv)
          1 - t(x) %*% iXstartXstar %*% x
        },
        DM = as.matrix(fitrMD$x)
      )

      sigmas <- sum(fitrMD$residuals * datarMD[, name.y]) / (fnx * fitrMD$df)
      #compute means:
      means <- as.matrix(newdataMD) %*% fitrMD$coefficients

      #simulated value for the new y's:
      if (dim(newdata)[1] == 1) {
        sigmaM <- as.matrix(sigmas, nr = 1, nc = 1)
      }
      if (dim(newdata)[1] > 1) {
        sigmaM <- diag(sigmas)
      }
      rpredictions[(max(cs.tmodels[iter - 1], 0) + 1):cs.tmodels[iter], ] <-
        rmvt(
          n = howmany,
          delta = means,
          sigma = sigmaM,
          df = fitrMD$df,
          type = "shifted"
        )

    }
  }

  return(rpredictions)

}
