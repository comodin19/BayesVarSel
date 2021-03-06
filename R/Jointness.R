#' Computation of Jointness measurements.
#'
#' \code{Jointness} computes the joint inclusion probabilitiy of two given
#' covariates as well as the jointness measurements of Ley and Steel (2007)
#'
#'
#' @export
#' @param x An object of class \code{Bvs}
#' @param covariates It can be either "All"(default) or a vector contaning the
#' name of two covariates.
#' @return An object of class \code{jointness} is returned.
#'
#' If \code{covariates} is "All" this object is a list with three matrices
#' containg different jointness measurements for all pairs of covariates is
#' returned.  In particular, for covariates i and j the jointness measurements
#' are:
#'
#' The Joint inclusion probabilities:
#'
#' \eqn{P(i and j)}
#'
#' And the two measurements of Ley and Steel (2007)
#'
#' \eqn{J*= P(i and j)/P(i or j)}
#'
#' \eqn{J*=P(i and j)/(P(i or j)-P(i and j))}
#'
#' If \code{covariates} is a vector of length 2, \code{Jointness} return a
#' list of four elements. The first three of them is a list of three values containing the
#' measurements above but just for the given pair of covariates. The fourth
#' element is the \code{covariates} vector.
#'
#' If method \code{print.jointness} is used a message with the meaning of the
#' measurement si printed.
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#' Maintainer: <anabel.forte@@uv.es>
#' @seealso \code{\link[BayesVarSel]{Bvs}} and
#' \code{\link[BayesVarSel]{GibbsBvs}} for performing variable selection and
#' obtaining an object of class \code{Bvs}.
#'
#' \code{\link[BayesVarSel]{plot.Bvs}} for different descriptive plots of the
#' results, \code{\link[BayesVarSel]{BMAcoeff}} for obtaining model averaged
#' simulations of regression coefficients and
#' \code{\link[BayesVarSel]{predict.Bvs}} for predictions.
#' @references
#'
#' Ley, E. and Steel, M.F.J. (2007)<DOI:10.1016/j.jmacro.2006.12.002>Jointness
#' in Bayesian variable selection with applications to growth regression.
#' Journal of Macroeconomics, 29(3):476-493.
#' @keywords package
#' @examples
#'
#' \dontrun{
#' #Analysis of Crime Data
#' #load data
#'
#' data(UScrime)
#'
#' crime.Bvs<- Bvs(formula= y ~ ., data=UScrime, n.keep=1000)
#'
#' #A look at the jointness measurements:
#' Jointness(crime.Bvs, covariates="All")
#'
#' Jointness(crime.Bvs, covariates=c("Ineq","Prob"))
#' #---------
#' #The joint inclusion probability for Ineq and Prob is:  0.65
#' #---------
#' #The ratio between the probability of including both
#' #covariates and the probability of including at least one of then is: 0.66
#' #---------
#' #The probability of including both covariates at the same times is 1.95 times
#' #the probability of including one of them alone
#'
#' }
Jointness <- function(x, covariates = "All") {
  if (!inherits(x, "Bvs"))
    warning("calling Jointness(<fake-Bvs-object>) ...")


  if (length(covariates) == 1 && covariates == "All") {
    joint <- x$jointinclprob
    joint_LS1 <- x$jointinclprob
    joint_LS2 <- x$jointinclprob
    for (i in 1:dim(x$jointinclprob)[1]) {
      for (j in i:dim(x$jointinclprob)[1]) {
        joint_LS1[i, j] <-
          x$jointinclprob[i, j] / (x$jointinclprob[i, i] + x$jointinclprob[j, j] -
                                     x$jointinclprob[i, j])

        joint_LS1[j, i] <-
          x$jointinclprob[i, j] / (x$jointinclprob[i, i] + x$jointinclprob[j, j] -
                                     x$jointinclprob[i, j])

        joint_LS2[i, j] <-
          x$jointinclprob[i, j] / (x$jointinclprob[i, i] + x$jointinclprob[j, j] -
                                     2 * x$jointinclprob[i, j])

        joint_LS2[j, i] <-
          x$jointinclprob[i, j] / (x$jointinclprob[i, i] + x$jointinclprob[j, j] -
                                     2 * x$jointinclprob[i, j])
      }
      joint_LS2[i, i] <- NA
    }


    jointness <-
      list("prob_joint" = joint,
           "joint_LS1" = joint_LS1,
           "joint_LS2" = joint_LS2)
    class(jointness) <- "jointness"

    return(jointness)
  }
  if (length(covariates) > 1) {
    if (length(covariates) > 2) {
      stop("The number of covariates to obtain jointness measurements should be 2\n")
    }
    if (length(covariates) < 2) {
      stop("The number of covariates to obtain jointness measurements should be 2\n")
    }
    i <- 0
    j <- 0

    i <- which(names(x$jointinclprob) == covariates[1])
    j <- which(names(x$jointinclprob) == covariates[2])

    if (i == 0 || j == 0) {
      stop("At least one of the covariates is not part of the analysis")
    }

    prob_joint <- x$jointinclprob[i, j]

    joint_LS1 <-
      x$jointinclprob[i, j] / (x$jointinclprob[i, i] + x$jointinclprob[j, j] -
                                 x$jointinclprob[i, j])

    joint_LS2 <-
      x$jointinclprob[i, j] / (x$jointinclprob[i, i] + x$jointinclprob[j, j] -
                                 2 * x$jointinclprob[i, j])



    jointness <-
      list("prob_joint" = prob_joint,
           "joint_LS1" = joint_LS1,
           "joint_LS2" = joint_LS2,
           "covariates" =  covariates)
    class(jointness) <- "jointness"

    return(jointness)
  }
}


#' Print an object of class \code{jointness}
#'
#' Print an object of class \code{jointness}. Show the different jointness measurements with a small explanation.
#'
#' @export
#' @param x An object of class \code{jointness}
#' @param ... Additional parameters to be passed
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#'   Maintainer: <anabel.forte@@uv.es>
#' @seealso See \code{\link[BayesVarSel]{Jointness}} for creating objects of the class
#'   \code{jointness}.
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
#' jointness(crime.Bvs)
#' }
#'
print.jointness <- function(x, ...) {
  if (length(x) == 4) {
    cat("---------\n")
    cat(
      paste(
        "The joint inclusion probability for",
        x[[4]][1],
        "and",
        x[[4]][2],
        "is: ",
        round(x[[1]], 2),
        "\n",
        sep = " "
      )
    )
    cat("---------\n")
    cat(
      paste(
        "The ratio between the probability of including both covariates and the probability of including at least one of then is: ",
        round(x[[2]], 2),
        "\n",
        sep = ""
      )
    )
    cat("---------\n")
    cat(
      paste(
        "The probability of including both covariates together is",
        round(x[[3]], 2),
        "times the probability of including one of them alone \n",
        sep = " "
      )
    )
  }
  if (length(x) == 3) {
    cat("---------\n")
    cat(paste(
      "The joint inclusion probability for All covariates are \n",
      sep = " "
    ))
    print(x[[1]])
    cat("---------\n")
    cat(
      paste(
        "The ratio between the probability of including two covariates together and the probability of including at least one of them is: \n",
        sep = ""
      )
    )
    print(x[[2]])
    cat("---------\n")
    cat(
      paste(
        "The ratio between the probability of including two covariates together and the probability of including one of them alone is: \n",
        sep = " "
      )
    )
    print(x[[3]])

  }
}
