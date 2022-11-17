#' A function for plotting summaries of an object of class \code{Bvs}
#'
#' Four different plots to summarize graphically the results in an object of
#' class \code{Bvs}.
#'
#' If \code{option}="dimension" this function returns a barplot of the
#' posterior distribution of the dimension of the true model. If
#' \code{option}="joint" an image plot of the joint inclusion probabilities is
#' returned. If \code{option}="conditional" an image plot of the conditional
#' inclusion probabilities. These should be read as the probabilty that the
#' variable in the column is part of the true model if the corresponding
#' variables on the row is. If \code{option}="not" the image plot that
#' is returned is that of the the probabilty that the variable in the column is
#' part of the true model if the corresponding variables on the row is not. Finally,
#' if \code{option}="trace", only available if x$method == "Gibbs", returns a plot of the trace of the inclusion probabilities to check for convergence.
#'
#' @export
#' @param x An object of class \code{Bvs}
#' @param option One of "dimension", "joint", "conditional", "not" or "trace"
#' @param ... Additional graphical parameters to be passed
#' @return If \code{option}="joint", "conditional" or "not" \code{plot} also
#' returns an object of class \code{matrix} with the numeric values of the
#' printed probabilities.
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#' Maintainer: <anabel.forte@@uv.es>
#' @seealso See \code{\link[BayesVarSel]{Bvs}}, \code{\link[BayesVarSel]{GibbsBvs}} for creating objects of the class
#' \code{Bvs}.
#' @examples
#'
#'
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
#' #An image plot of the joint inclusion probabilities:
#' plot(crime.Bvs, option="joint")
#'
#' #Two image plots of the conditional inclusion probabilities:
#' plot(crime.Bvs, option="conditional")
#' plot(crime.Bvs, option="not")
#'
#'
plot.Bvs <-
  function(x, option = "dimension", ...) {
    #if (!inherits(x, "Bvs"))
    #  warning("calling plotBvs(<fake-Bvs-object>) ...")
    p <- x$p
    k <- x$k
    auxtp <- substr(tolower(option), 1, 1)

    if (auxtp != "d" && auxtp != "j" && auxtp != "c" && auxtp != "n"&& auxtp != "t") {
      stop("I am very sorry: type of plot not specified\n")
    }

    #dimension probabilities
    if (auxtp == "d") {
      par(mar = c(5, 4, 4, 2) + 0.1, mfrow = c(1, 1))
      if (x$method == "gibbs") {
        barplot(
          x$postprobdim,
          main = "Estimated Posterior Dimension Probabilities",
          xlab = "Number of covariates in the true model",
          ylab = "Probability",
          names.arg = (0:p) + k,
          ...
        )
      } else{
        barplot(
          x$postprobdim,
          main = "Posterior Dimension Probabilities",
          xlab = "Number of covariates in the true model",
          ylab = "Probability",
          names.arg = (0:p) + k,
          ...
        )
      }
    }
    #Special function for printing the joint and conditional probabilities.


    myImagePlot <- function(x, scale, ...) {
      x <- as.matrix(x)


      #What do we do with the diagonal?
      #diag(x)<- diagonal
      #What do we do with the NA's?
      x[is.na(x)] <- 0

      if (sum(is.na(x)) > 0)
        warning("The matrix contains NA's. The corresponding values have been set to 0")

      min <- min(x)
      max <- max(x)
      yLabels <- rownames(x)
      xLabels <- colnames(x)
      title <- c()
      # check for additional function arguments
      if (length(list(...))) {
        Lst <- list(...)
        if (!is.null(Lst$zlim)) {
          min <- Lst$zlim[1]
          max <- Lst$zlim[2]
        }
        if (!is.null(Lst$yLabels)) {
          yLabels <- c(Lst$yLabels)
        }
        if (!is.null(Lst$xLabels)) {
          xLabels <- c(Lst$xLabels)
        }
        if (!is.null(Lst$title)) {
          title <- Lst$title
        }
      }
      # check for null values
      if (is.null(xLabels)) {
        xLabels <- c(1:ncol(x))
      }
      if (is.null(yLabels)) {
        yLabels <- c(1:nrow(x))
      }

      #layout(matrix(data=c(1,2), nrow=2, ncol=1), widths=c(1,1), heights=c(1,4.5))
      layout(
        matrix(
          data = c(1, 2, 3, 4),
          nrow = 2,
          ncol = 2
        ),
        widths = c(5, 1),
        heights = c(1, 4.5)
      )

      ColorRamp <- gray.colors(40)[40:0]
      ColorLevels <- seq(0, 1, length = length(ColorRamp))

      # Reverse Y axis
      reverse <- nrow(x):1
      yLabels <- yLabels[reverse]
      x <- x[reverse, ]

      # Inclusion probs:
      par(mar = c(3, 5, 2.5, 2))
      image(
        1:length(scale),
        1,
        matrix(
          data = scale,
          ncol = 1,
          nrow = length(scale)
        ),
        xlab = "",
        ylab = "",
        col = ColorRamp, breaks=(0:length(ColorRamp))/length(ColorRamp),
        yaxt = "n",
        xaxt = "n"
        )

      if (!is.null(title)) {
        title(main = title)

      }# Title on the top

      # Data Map
      par(mar = c(3, 5, 2.5, 2))
      image(
        1:length(xLabels),
        1:length(yLabels),
        t(x),
        col = ColorRamp, breaks=(0:length(ColorRamp))/length(ColorRamp),
        xlab = "",
        ylab = "",
        axes = FALSE,
        zlim = c(min, max)
        )

      axis(
        side = 3,
        at = 1:length(xLabels),
        labels = xLabels,
        cex.axis = 0.7,
        las = 2
      )
      axis(
        side = 2,
        at = 1:length(yLabels),
        labels = yLabels,
        las = 1,
        cex.axis = 0.7
      )

      #Nothing
      par(mar = c(3, 5, 2.5, 2))
      plot(
        1,
        1,
        type = "n",
        axes = 0,
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = ""
      )

      # color scale:
      par(mar = c(3, 5, 2.5, 2))
      par(mar = c(3, 2.5, 2.5, 2))
      image(
        1,
        ColorLevels,
        matrix(
          data = ColorLevels,
          ncol = length(ColorLevels),
          nrow = 1
        ),
        col = ColorRamp,breaks=(0:length(ColorRamp))/length(ColorRamp),
        xlab = "",
        ylab = "",
        xaxt = "n"
      )
    }

    #Joint posterior probabilities
    if (auxtp == "j") {
      if (x$method == "gibbs") {
        myImagePlot(x$jointinclprob,
                    scale = as.vector(x$inclprob),
                    title = "Estimated Joint Inclusion Probabilities")
      } else{
        myImagePlot(x$jointinclprob,
                    scale = as.vector(x$inclprob),
                    title = "Joint Inclusion Probabilities")
      }
      prob_joint <- as.matrix(x$jointinclprob)
      return(invisible(prob_joint))
    }



    #conditional posterior probabilities
    if (auxtp == "c") {
      prob_joint <- as.matrix(x$jointinclprob)
      prob_cond <- as.matrix(x$jointinclprob)
      for (i in 1:p) {
        prob_cond[i, ] <- prob_joint[i, ] / prob_joint[i, i]
        prob_cond[i, i] <- 1
      }

      if (x$method == "gibbs") {
        myImagePlot(
          prob_cond,
          scale = as.vector(x$inclprob),
          diagonal = 1,
          title = "Estimated Inclusion prob. of column var. given the row var. is included")
      } else{
        myImagePlot(
          prob_cond,
          scale = as.vector(x$inclprob),
          diagonal = 1,
          title = "Inclusion prob. of column var. given the row var. is included")
      }
      return(invisible(prob_cond))
    }
    #conditional posterior probabilities given Not a variable
    if (auxtp == "n") {
      #auxiliar function to compute de A Given not B matrix
      pAgivenNotB <- function(miobject) {
        result <- 0 * miobject$jointinclprob
        for (i in 1:dim(result)[1]) {
          if (miobject$inclprob[i] == 1) {
            warning(
              paste(
                "The inclusion probabilities of",
                miobject$variables[i],
                "is equal to 1. The conditional posterior probability of any other variable given it can not be computed. Cero is returned instead\n",
                sep = " "
              )
            )
            for (j in 1:dim(result)[1]) {
              result[i, j] <- 0
            }
          } else{
            for (j in 1:dim(result)[1]) {
              result[i, j] <-
                (miobject$inclprob[j] - miobject$jointinclprob[i, j]) / (1 -
                                                                           miobject$inclprob[i])#REVISAR  P(j|Not.i)=(1-P(i|j))P(j)/(1-P(i))=(P(j)-P(j,i))/(1-P(i))
            }
          }
          result[i, i] <- 0
        }
        colnames(result) <- colnames(miobject$jointinclprob)
        rownames(result) <-
          paste("Not.", colnames(miobject$jointinclprob), sep = "")
        result
      }
      AgivenNotB <- pAgivenNotB(x)

      if (x$method == "gibbs") {
        myImagePlot(AgivenNotB,
                    scale = as.vector(x$inclprob),
                    title = "Est. Incl. prob of column var. given the row var. is NOT included")
        } else{
        myImagePlot(AgivenNotB,
                    scale = as.vector(x$inclprob),
                    title = "Incl. prob of column var. given the row var. is NOT included")
      }
      prob_not <- AgivenNotB
      return (invisible(prob_not))
    }

    #conditional posterior probabilities given Not a variable
    if (auxtp == "t") {
      if (!x$method == "gibbs")
        stop("For convergence graphics use an object from GibbsBvs")

      nmodels <- dim(x$modelslogBF)[1]

      aux <- matrix(rep(1 / 1:nmodels, nmodels), ncol = nmodels, nrow = nmodels)
      aux2 <- matrix(0, ncol = nmodels, nrow = nmodels)
      aux2[lower.tri(aux2)] <- aux[lower.tri(aux)]
      #aux2 used for the calculations of the partial means
      nvar <- dim(x$modelslogBF)[2] - 1

      incprobiter <- data.frame(t(rep(NA, nvar)))
      incprobiter <- aux2 %*% x$modelslogBF[, 1:nvar]

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
    }
  }
