#' Convergence plot for the inclusion probabilities of a Bvs (Gibbs) object.
#'
#'Returns a picture of the trace of the inclusion probabilities of the potential variables
#'
#' @export
#'
#' @param object An object of class \code{Bvs} obtained with \code{GibbsBvs}
#' @param variables A vector containing the name of the variables that would be plotted (all of them by default)
#'
#' @return \code{convPlot} returns a picture of the trace of the inclusion probabilities of the variables included in \code{variable}
#' @author Anabel Forte and Gonzalo Garcia-Donato
#'
#' @examples
#' \dontrun{
#' #Analysis of Ozone35 data
#'
#' data(Ozone35)
#'
#' #We use here the (Zellner) g-prior for
#' #regression parameters and constant prior
#' #over the model space
#' #In this Gibbs sampling scheme, we perform 10100 iterations,
#' #of which the first 100 are discharged (burnin) and of the remaining
#' #only one each 10 is kept.
#' #as initial model we use the Full model
#' Oz35.GibbsBvs<- GibbsBvs(formula="y~.", data=Ozone35, prior.betas="gZellner",
#' prior.models="Constant", n.iter=10000, init.model="Full", n.burnin=100,
#' time.test = FALSE)
#'
#' #Convergence plot for all the variables
#' convPlot(Oz35.GibbsBvs)
#'
#' #Convergence plot for x4 x5 and x6
#' convPlot(Oz35.GibbsBvs, variables = c("x4", "x5", "x6"))
#' }
convPlot <- function(object, variables = "All") {
  if (!inherits(object, "Bvs"))
    warning("calling Jointness(<fake-Bvs-object>) ...")
  if (!object$method == "gibbs")
    stop("For convergence graphics use an object from GibbsBvs")

  nmodels <- dim(object$modelslogBF)[1]

  aux <- matrix(rep(1 / 1:nmodels, nmodels), ncol = nmodels, nrow = nmodels)
  aux2 <- matrix(0, ncol = nmodels, nrow = nmodels)
  aux2[lower.tri(aux2)] <- aux[lower.tri(aux)]
  #aux2 used for the calculations of the partial means

  if ("All" %in% variables) {
    nvar <- dim(object$modelslogBF)[2] - 1

    incprobiter <- data.frame(t(rep(NA, nvar)))
    incprobiter <- aux2 %*% object$modelslogBF[, 1:nvar]

    graphics::plot(
      1:nmodels,
      incprobiter[, 1],
      type =  "l",
      ylim = c(0, 1),
      ylab = "Inclusion Probability",
      xlab = "n.iter"
    )
    for (i in 2:nvar)
      graphics::lines(1:nmodels, incprobiter[, i], col = i)
  } else{
    if (sum(variables %in% names(as.data.frame(object$modelslogBF))) < length(variables)) {
      stop("At least one of the variables is not among the potential ones")
    }
    nvar <- length(variables)

    incprobiter <- data.frame(t(rep(NA, nvar)))
    incprobiter <- aux2 %*% object$modelslogBF[, variables]

    plot(
      1:nmodels,
      incprobiter[, 1],
      type =  "l",
      ylim = c(0, 1),
      ylab = "Inclusion Probability",
      xlab = "n.iter"
    )
    if (nvar >= 2) {
      for (i in 2:nvar)
        graphics::lines(1:nmodels, incprobiter[, i], col = i)
    }
    graphics::legend(0,
           1,
           legend = variables,
           lty = 1,
           col = 1:nvar)
  }
}
