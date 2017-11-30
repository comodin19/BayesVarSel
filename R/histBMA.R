#' A function for histograms-like representations of objects of class
#' \code{bma.coeffs}
#'
#' The columns in \code{bma.coeffs} are simulations of the model averaged
#' posterior distribution. This normally is a mixture of a discrete (at zero)
#' and several continuous distributions. This plot provides a convenient
#' graphical summary of such distributions.
#'
#' This function produces a histogram but with the peculiarity that the zero
#' values in the simulation are represented as bar centered at zero. The area
#' of all the bars is one and of these, the area of the bar at zero (colored
#' with \code{gray.0}) is, conditionally on the retained models (see details in
#' \code{\link[BayesVarSel]{BMAcoeff}}), the probability of that coefficient be
#' exactly zero. This number is included in the top of the zero bar if
#' \code{text} is set to TRUE.
#'
#' @export
#' @param x An object of class \code{bma.coeffs}
#' @param covariate The name of an explanatory variable whose accompanying
#' coefficient is to be represented. This must be the name of one of the
#' columns in \code{x}
#' @param n.breaks The number of equally lentgh bars for the histogram
#' @param text If set to TRUE the probability of the coefficient being zero is
#' added in top of the bar at zero. Note: this probability is based on the
#' models used in \code{bma.coeffs} (see details in that function)
#' @param gray.0 A numeric value between 0 and 1 that specifies the darkness,
#' in a gray scale (0 is white and 1 is black) of the bar at zero
#' @param gray.no0 A numeric value between 0 and 1 that specifies the darkness,
#' in a gray scale (0 is white and 1 is black) of the bars different from zero
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#' Maintainer: <anabel.forte@@uv.es>
#' @seealso See \code{\link[BayesVarSel]{BMAcoeff}}. Also see
#' \code{\link[BayesVarSel]{Bvs}} and
#' \code{\link[BayesVarSel]{GibbsBvs}} for creating objects of the class
#' \code{BMAcoeff}.
#' @examples
#'
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
#' #Observe the bimodality of the coefficient associated with regressor M
#' histBMA(crime.Bvs.BMA, "M")
#'
#' #Note 1:
#' #The value in top of the bar at zero (0.251 in this case) is the probability of beta_M is
#' #zero conditional on a model space containing the 1000 models used in the mixture. This value
#' #should be closed to the exact value
#' #1-crime.Bvs$inclprob["M"]
#' #which in this case is 0.2954968
#' #if n.keep above is close to 2^15
#'
#' #Note 2:
#' #The BMA posterior distribution of beta_M has two modes approximately located at 0 and 10
#' #If we summarize this distribution using the mean
#' mean(crime.Bvs.BMA[ ,"M"])
#' #or median
#' median(crime.Bvs.BMA[ ,"M"])
#' #we obtain values around 7 (or 7.6) which do not represent this distribution.
#'
#' #With the Gibbs algorithms:
#' data(Ozone35)
#'
#' Oz35.GibbsBvs<- GibbsBvs(formula="y~.", data=Ozone35, prior.betas="gZellner",
#' prior.models="Constant", n.iter=10000, init.model="Full", n.burnin=100,
#' time.test = FALSE)
#' Oz35.GibbsBvs.BMA<- BMAcoeff(Oz35.GibbsBvs, n.sim=10000)
#'
#' histBMA(Oz35.GibbsBvs.BMA, "x6.x7")
#' #In this case (Gibbs sampling), the value in top of the bar at zero (0.366)
#' #basically should coincide (if n.sim is large enough)
#' #with the estimated complement of the inclusion probability
#' #1-Oz35.GibbsBvs$inclprob["x6.x7"]
#' #which in this case is 0.3638
#' }
#'
histBMA <-
  function(x,
           covariate,
           n.breaks = 100,
           text = TRUE,
           gray.0 = 0.6,
           gray.no0 = 0.8) {
    #x must be an object of class bma.coeffs
    if (!inherits(x, "bma.coeffs")) {
      stop("x must be an object of class bma.coeffs")
    }
    #covariate must be a variable in the analysis
    i <- which(colnames(x) == covariate)
    if (sum(colnames(x) == covariate) == 0) {
      stop(paste("covariate ", covariate, " is not in the analysis", sep = ""))
    }

    v <- x[, covariate]

    igual <- sum(v == 0)
    menor <- sum(v < 0)
    mayor <- sum(v > 0)
    total <- length(v)
    dist <- (max(v) - min(v)) / (n.breaks)
    if (menor != 0 & mayor != 0 & igual != 0) {
      cortes <-
        c(
          seq(
            from = min(v),
            to = 0 - dist / 2,
            by = dist
          ),
          0 - dist / 2,
          seq(0 + dist / 2, max(v) + dist, by = dist)
        )
      cortes <- unique(cortes)
      h <- hist(v, breaks = cortes, plot = F)
      cutoff <- cut(h$breaks, breaks = c(-Inf, 0 - dist / 2, 0 + dist / 2, +Inf))
      color <-
        c(rep(gray(gray.no0), sum(as.numeric(cutoff) == 1) - 1), rep(gray(gray.0), 1), rep(gray(gray.no0), sum(as.numeric(cutoff) ==
                                                                                                                 3)))
      plot(
        h,
        col = color,
        freq = FALSE,
        yaxt = "n",
        main = covariate,
        xlab = "",
        ylab = ""
      )
      if (text)
        text(
          x = 0,
          y = h$density[h$mids == 0] + 0.025 * h$density[h$mids == 0],
          labels = round(igual / total, 3)
        )
    }
    if (menor == 0 & igual != 0 & mayor != 0) {
      cortes <-
        c(0 - dist / 2,
          0 + dist / 2,
          seq(
            from = 0 + dist / 2,
            to = max(v) + dist,
            by = dist
          ))
      cortes <- unique(cortes)
      h <- hist(v, breaks = cortes, plot = F)
      cutoff <- cut(h$breaks, breaks = c(-Inf, 0 - dist / 2, 0 + dist / 2, +Inf))
      color <-
        c(rep(gray(gray.no0), sum(as.numeric(cutoff) == 1) - 1), rep(gray(gray.0), 1), rep(gray(gray.no0), sum(as.numeric(cutoff) ==
                                                                                                                 3)))
      plot(
        h,
        col = color,
        freq = FALSE,
        yaxt = "n",
        xlab = "",
        ylab = "",
        main = covariate
      )
      if (text)
        text(
          x = 0,
          y = h$density[h$mids == 0] + 0.025 * h$density[h$mids == 0],
          labels = round(igual / total, 3)
        )
    }
    if (mayor == 0 & igual != 0 & menor != 0) {
      cortes <-
        c(seq(
          from = min(v),
          to = 0 - dist / 2,
          by = dist
        ), 0 - dist / 2, 0 + dist / 2)
      cortes <- unique(cortes)
      h <- hist(v, breaks = cortes, plot = F)
      cutoff <- cut(h$breaks, breaks = c(-Inf, 0 - dist / 2, 0 + dist / 2, +Inf))
      color <-
        c(rep(gray(gray.no0), sum(as.numeric(cutoff) == 1) - 1), rep(gray(gray.0), 1), rep(gray(gray.no0), sum(as.numeric(cutoff) ==
                                                                                                                 3)))
      plot(
        h,
        col = color,
        freq = FALSE,
        yaxt = "n",
        xlab = "",
        ylab = "",
        main = covariate
      )
      if (text)
        text(
          x = 0,
          y = h$density[h$mids == 0] + 0.025 * h$density[h$mids == 0],
          labels = round(igual / total, 3)
        )
    }

    if (igual == 0) {
      hist(
        v,
        breaks = 100,
        freq = FALSE,
        col = gray(gray.no0),
        yaxt = "n",
        xlab = "",
        ylab = "",
        main = covariate
      )
    }
    if (igual == total) {
      print(
        paste(
          "The regression coefficient associated to",
          covariate,
          "is 0 in all the models used for BMA"
        )
      )
    }
  }
