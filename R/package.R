
#' @import parallel MASS mvtnorm
#' @useDynLib BayesVarSel
#' @importFrom grDevices gray.colors gray
#' @importFrom graphics axis barplot image layout par plot hist text
#' @importFrom stats as.formula lm rbinom runif density quantile
#' @importFrom utils read.table write.table
NULL


#' Bayes Factors, Model Choice And Variable Selection In Linear Models
#'
#' Hypothesis testing, model selection and model averaging are important
#' statistical problems that have in common the explicit consideration of the
#' uncertainty about which is the true model. The formal Bayesian tool to solve
#' such problems is the Bayes factor (Kass and Raftery, 1995) that reports the
#' evidence in the data favoring each of the entertained hypotheses/models and
#' can be easily translated to posterior probabilities.
#'
#' This package has been specifically conceived to calculate Bayes factors in
#' linear models and then to provide a formal Bayesian answer to testing and
#' variable selection problems. From a theoretical side, the emphasis in the
#' package is placed on the prior distributions (a very delicate issue in this
#' context) and BayesVarSel allows using a wide range of them:
#' Jeffreys-Zellner-Siow (Jeffreys, 1961; Zellner and Siow, 1980,1984) Zellner
#' (1986); Fernandez et al. (2001), Liang et al. (2008) and Bayarri et al.
#' (2012).
#'
#' The interaction with the package is through a friendly interface that
#' syntactically mimics the well-known lm command of R. The resulting objects
#' can be easily explored providing the user very valuable information (like
#' marginal, joint and conditional inclusion probabilities of potential
#' variables; the highest posterior probability model, HPM; the median
#' probability model, MPM) about the structure of the true -data generating-
#' model. Additionally, BayesVarSel incorporates abilities to handle problems
#' with a large number of potential explanatory variables through parallel and
#' heuristic versions (Garcia-Donato and Martinez-Beneito 2013) of the main
#' commands.
#'
#' \tabular{ll}{ Package: \tab BayesVarSel\cr Type: \tab Package\cr Version:
#' \tab 1.7.0\cr Date: \tab 2016-08-31\cr License: \tab GPL-2\cr }
#'
#' @name BayesVarSel-package
#' @aliases BayesVarSel-package BayesVarSel
#' @docType package
#' @author Gonzalo Garcia-Donato and Anabel Forte
#'
#' Maintainer: Anabel Forte \email{anabel.forte@@uv.es}
#' @seealso \code{\link[BayesVarSel]{Btest}}, \code{\link[BayesVarSel]{Bvs}},
#' \code{\link[BayesVarSel]{PBvs}}, \code{\link[BayesVarSel]{GibbsBvs}},
#' \code{\link[BayesVarSel]{BMAcoeff}}, \code{\link[BayesVarSel]{predictBvs}}
#' @references
#'
#' Bayarri, M.J., Berger, J.O., Forte, A. and Garcia-Donato, G.
#' (2012)<DOI:10.1214/12-aos1013> Criteria for Bayesian Model choice with
#' Application to Variable Selection. The Annals of Statistics. 40: 1550-1577
#'
#' Fernandez, C., Ley, E. and Steel, M.F.J.
#' (2001)<DOI:10.1016/s0304-4076(00)00076-2> Benchmark priors for Bayesian
#' model averaging. Journal of Econometrics, 100, 381-427.
#'
#' Garcia-Donato, G. and Martinez-Beneito, M.A.
#' (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian
#' variable selection problems with large model spaces. Journal of the American
#' Statistical Association. 108: 340-352.
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O.
#' (2008)<DOI:10.1198/016214507000001337> Mixtures of g-priors for Bayesian
#' Variable Selection. Journal of the American Statistical Association.
#' 103:410-423.
#'
#' Zellner, A. and Siow, A. (1980)<DOI:10.1007/bf02888369>. Posterior Odds
#' Ratio for Selected Regression Hypotheses. In Bayesian Statistics 1 (J.M.
#' Bernardo, M. H. DeGroot, D. V. Lindley and A. F. M. Smith, eds.) 585-603.
#' Valencia: University Press.
#'
#' Zellner, A. and Siow, A. (1984) Basic Issues in Econometrics. Chicago:
#' University of Chicago Press.
#'
#' Zellner, A. (1986)<DOI:10.2307/2233941> On Assessing Prior Distributions and
#' Bayesian Regression Analysis with g-prior Distributions. In Bayesian
#' Inference and Decision techniques: Essays in Honor of Bruno de Finetti (A.
#' Zellner, ed.) 389-399. Edward Elgar Publishing Limited.
#' @keywords package
#' @examples
#' demo(BayesVarSel.Hald)
#'
NULL


#' Hald data
#'
#' The following data relates to an engineering application that was interested
#' in the effect of the cement composition on heat evolved during hardening
#' (for more details, see Woods et al., 1932).
#'
#'
#' @name Hald
#' @docType data
#' @format A data frame with 13 observations on the following 5 variables.
#' \describe{ \item{y}{Heat evolved per gram of cement (in calories)}
#' \item{x1}{Amount of tricalcium aluminate} \item{x2}{Amount
#' of tricalcium silicate} \item{x3}{Amount of tetracalcium alumino
#' ferrite} \item{x4}{Amount of dicalcium silicate} }
#' @references Woods, H., Steinour, H. and Starke, H.
#' (1932)<DOI:10.1021/ie50275a002> Effect of Composition of Porland Cement on
#' Heat Evolved During Hardening. Industrial and Engineering Chemistry
#' Research, 24, 1207-1214.
#' @keywords datasets
#' @examples
#' data(Hald)
#'
"Hald"


#' Ozone35 dataset
#'
#' Polution data
#'
#' This dataset has been used by Garcia-Donato and Martinez-Beneito (2013) to
#' illustrate the potential of the Gibbs sampling method (in \code{BayesVarSel}
#' implemented in \code{\link[BayesVarSel]{GibbsBvs}}).
#'
#' This data were previously used by Casella and Moreno (2006) and Berger and
#' Molina (2005) and concern N = 178 measures of ozone concentration in the
#' atmosphere. Of the 10 main effects originally considered, we only make use
#' of those with an atmospheric meaning x4 to x10, as was done by Liang et al.
#' (2008). We then have 7 main effects which, jointly with the quadratic terms
#' and second order interactions, produce the above-mentioned p = 35 possible
#' regressors.
#'
#' @name Ozone35
#' @docType data
#' @format A data frame with 178 observations on the following 36 variables.
#' \describe{ \item{y}{Response = Daily maximum 1-hour-average ozone
#' reading (ppm) at Upland, CA} \item{x4}{500-millibar pressure height
#' (m) measured at Vandenberg AFB} \item{x5}{Wind speed (mph) at Los
#' Angeles International Airport (LAX)} \item{x6}{Humidity (percentage)
#' at LAX} \item{x7}{Temperature (Fahrenheit degrees) measured at
#' Sandburg, CA} \item{x8}{Inversion base height (feet) at LAX}
#' \item{x9}{Pressure gradient (mm Hg) from LAX to Daggett, CA}
#' \item{x10}{Visibility (miles) measured at LAX}
#' \item{x4.x4}{=x4*x4} \item{x4.x5}{=x4*x5}
#' \item{x4.x6}{=x4*x6} \item{x4.x7}{=x4*x7}
#' \item{x4.x8}{=x4*x8} \item{x4.x9}{=x4*x9}
#' \item{x4.x10}{=x4*x10} \item{x5.x5}{=x5*x5}
#' \item{x5.x6}{=x5*x6} \item{x5.x7}{=x5*x7}
#' \item{x5.x8}{=x5*x8} \item{x5.x9}{=x5*x9}
#' \item{x5.x10}{=x5*x10} \item{x6.x6}{=x6*x6}
#' \item{x6.x7}{=x6*x7} \item{x6.x8}{=x6*x8}
#' \item{x6.x9}{=x6*x9} \item{x6.x10}{=x6*x10}
#' \item{x7.x7}{=x7*x7} \item{x7.x8}{=x7*x8}
#' \item{x7.x9}{=x7*x9} \item{x7.x10}{=x7*x10}
#' \item{x8.x8}{=x8*x8} \item{x8.x9}{=x8*x9}
#' \item{x8.x10}{=x8*x10} \item{x9.x9}{=x9*x9}
#' \item{x9.x10}{=x9*x10} \item{x10.x10}{=x10*x10} }
#' @references Berger, J. and Molina, G. (2005)<DOI:j.1467-9574.2005.00275.x>
#' Posterior model probabilities via path-based pairwise priors. Statistica
#' Neerlandica, 59:3-15.
#'
#' Casella, G. and Moreno, E. (2006)<DOI:10.1198/016214505000000646> Objective
#' Bayesian variable selection. Journal of the American Statistical
#' Association, 101(473).
#'
#' Garcia-Donato, G. and Martinez-Beneito, M.A.
#' (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian
#' variable selection problems with large model spaces. Journal of the American
#' Statistical Association, 108: 340-352.
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O.
#' (2008)<DOI:10.1198/016214507000001337> Mixtures of g-priors for Bayesian
#' Variable Selection. Journal of the American Statistical Association.
#' 103:410-423.
#' @keywords datasets
#' @examples
#' data(Ozone35)
#'
"Ozone35"
