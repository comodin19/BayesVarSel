#' Bayesian Variable Selection for linear regression models using Gibbs
#' sampling.
#'
#' Approximate computation of summaries of the posterior distribution using a
#' Gibbs sampling algorithm to explore the model space and frequency of
#' "visits" to construct the estimates.
#'
#' This is a heuristic approximation to the function
#' \code{\link[BayesVarSel]{Bvs}} so the details there apply also here.
#'
#' The algorithm implemented is a Gibbs sampling-based searching algorithm
#' originally proposed by George and McCulloch (1997). Garcia-Donato and
#' Martinez-Beneito (2013) have shown that this simple sampling strategy in
#' combination with estimates based on frequency of visits (the one here
#' implemented) provides very reliable results.
#' @export
#' @param formula Formula defining the most complex regression model in the
#' analysis. See details.
#' @param data data frame containing the data.
#' @param null.model A formula defining which is the simplest (null) model.
#' It should be nested in the full model. By default, the null model is defined
#' to be the one with just the intercept.
#' @param prior.betas Prior distribution for regression parameters within each
#' model. Possible choices include "Robust", "Liangetal", "gZellner",
#' "ZellnerSiow" and "FLS" (see details).
#' @param prior.models Prior distribution over the model space. Possible
#' choices are "Constant", "ScottBerger" and "User" (see details).
#' @param n.iter The total number of iterations performed after the burn in
#' process.
#' @param init.model The model at which the simulation process starts. Options
#' include "Null" (the model only with the covariates specified in
#' \code{fixed.cov}), "Full" (the model defined by \code{formula}), "Random" (a
#' randomly selected model) and a vector with p (the number of covariates to
#' select from) zeros and ones defining a model.
#' @param n.burnin Length of burn in, i.e. number of iterations to discard at
#' the beginning.
#' @param n.thin Thinning rate. Must be a positive integer.  Set 'n.thin' > 1
#' to save memory and computation time if 'n.iter' is large. Default is 1. This
#' parameter jointly with \code{n.iter} sets the number of simulations kept and
#' used to construct the estimates so is important to keep in mind that a large
#' value for 'n.thin' can reduce the precision of the results
#' @param time.test If TRUE and the number of variables is large (>=21) a
#' preliminary test to estimate computational time is performed.
#' @param priorprobs A p+1 dimensional vector defining the prior probabilities
#' Pr(M_i) (should be used in the case where \code{prior.models}="User"; see
#' the details in \code{\link[BayesVarSel]{Bvs}}.)
#' @param seed A seed to initialize the random number generator
#' @return \code{GibbsBvs} returns an object of class \code{Bvs} with the
#' following elements: \item{time }{The internal time consumed in solving the
#' problem} \item{lmfull }{The \code{lm} class object that results when the
#' model defined by \code{formula} is fitted by \code{lm}} \item{lmnull }{The
#' \code{lm} class object that results when the model defined by
#' \code{fixed.cov} is fitted by \code{lm}} \item{variables }{The name of all
#' the potential explanatory variables} \item{n }{Number of observations}
#' \item{p }{Number of explanatory variables to select from} \item{k }{Number
#' of fixed variables} \item{HPMbin }{The binary expression of the most
#' probable model found.} \item{inclprob }{A \code{data.frame} with the
#' estimates of the inclusion probabilities of all the variables.}
#' \item{jointinclprob }{A \code{data.frame} with the estimates of the joint
#' inclusion probabilities of all the variables.} \item{postprobdim }{Estimates
#' of posterior probabilities of the dimension of the true model.}
#' \item{modelslogBF}{A matrix with both the binary representation of the
#' visited models after the burning period and the Bayes factor (log scale) of
#' that model to the null model.}\item{priorprobs}{A p+1 dimensional vector containing values proportionals
#' to the prior probability of a model of each dimension (from 0 to p)} \item{call }{The \code{call} to the
#' function.} 
#' \item{C}{An estimation of the normalizing constant (C=sum BiPr(Mi), for Mi in the model space)}
#' \item{method }{\code{gibbs}}
#' @author Gonzalo Garcia-Donato and Anabel Forte
#' @seealso \code{\link[BayesVarSel]{plot.Bvs}} for several plots of the result,
#' \code{\link[BayesVarSel]{BMAcoeff}} for obtaining model averaged simulations
#' of regression coefficients and \code{\link[BayesVarSel]{predict.Bvs}} for
#' predictions.
#'
#' \code{\link[BayesVarSel]{Bvs}} for exact
#' version obtained enumerating all entertained models (recommended when
#' p<20).
#' @references Garcia-Donato, G. and Martinez-Beneito, M.A.
#' (2013)<DOI:10.1080/01621459.2012.742443> On sampling strategies in Bayesian
#' variable selection problems with large model spaces. Journal of the American
#' Statistical Association, 108: 340-352.
#'
#' George E. and McCulloch R. (1997) Approaches for Bayesian variable
#' selection. Statistica Sinica, 7, 339:372.
#' @keywords package
#' @examples
#'
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
#' Oz35.GibbsBvs<- GibbsBvs(formula= y ~ ., data=Ozone35, prior.betas="gZellner",
#' prior.models="Constant", n.iter=10000, init.model="Full", n.burnin=100,
#' time.test = FALSE)
#'
#' #Note: this is a heuristic approach and results are estimates
#' #of the exact answer.
#'
#' #with the print we can see which is the most probable model
#' #among the visited
#' Oz35.GibbsBvs
#'
#' #The estimation of inclusion probabilities and
#' #the model-averaged estimation of parameters:
#' summary(Oz35.GibbsBvs)
#'
#' #Plots:
#' plot(Oz35.GibbsBvs, option="conditional")
#' }
#'
GibbsBvs <-
  function(formula,
           data,
           null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
           prior.betas = "Robust",
           prior.models = "ScottBerger",
           n.iter = 10000,
           init.model = "Full",
           n.burnin = 500,
           n.thin = 1,
           time.test = TRUE,
           priorprobs = NULL,
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
    if (!is.null(fixed.cov)) {
      #Eval the full model
      lmfull = lm(formula,
                  data = data,
                  y = TRUE,
                  x = TRUE)
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

      p <- dim(X)[2]#Number of covariates to select from

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
    }

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
    cat("The problem has a total of", 2 ^ (p), "competing models\n")
    iter <- n.iter
    cat("Of these,", n.burnin + n.iter, "are sampled with replacement\n")

    cat("Then,",
        floor(iter / n.thin),
        "are kept and used to construct the summaries\n")


    #prior for betas:
    pfb <- substr(tolower(prior.betas), 1, 1)
    if (pfb != "g" &&
        pfb != "r" &&
        pfb != "z" &&
        pfb != "l" &&
        pfb != "f")
      stop("I am very sorry: prior for betas not supported\n")
    #prior for model space:
    pfms <- substr(tolower(prior.models), 1, 1)
    if (pfms != "c" &&
        pfms != "s" &&
        pfms != "u")
      stop("I am very sorry: prior for model space not valid\n")
    if (pfms == "u" &&
        is.null(priorprobs)) {
      stop("A valid vector of prior probabilities must be provided\n")
    }
    if (pfms == "u" &&
        length(priorprobs) != (p + 1)) {
      stop("Vector of prior probabilities with incorrect length\n")
    }
    if (pfms == "u" &&
        sum(priorprobs < 0) > 0) {
      stop("Prior probabilities must be positive\n")
    }
    if (pfms == "u" &&
        priorprobs[1] == 0) {
      stop(
        "Vector of prior probabilities not valid: All the theory here implemented works with the implicit assumption that the null model could be the true model\n"
      )
    }
    if (pfms == "u" &&
        priorprobs[sum(init.model) + 1] == 0) {
      stop("The initial model has zero prior probability\n")
    }
    if (pfms == "u") {
      #The zero here added is for C compatibility
      write(
        priorprobs,
        ncolumns = 1,
        file = paste(wd, "/priorprobs.txt", sep = "")
      )
    }

    #Note: priorprobs.txt is a file that is needed only by the "User" routine. Nevertheless, in order
    #to mantain a common unified version the source files of other routines also reads this file
    #although they do not use. Because of this we create this file anyway.
    if (pfms == "c" | pfms == "s") {
      priorprobs <- rep(0, p + 1)
      write(
        priorprobs,
        ncolumns = 1,
        file = paste(wd, "/priorprobs.txt", sep = "")
      )
    }

    method <- paste(pfb, pfms, sep = "")

    #The previous test (for time)
    estim.time <- 0
    if (p <= 20) {
      warning(
        "The number of variables is small enough to visit every model. Consider Bvs (or pBvs for its parallel version).\n"
      )
    }
    if (time.test && p > 20) {
      cat("Time test. . . .\n")


      result <- switch(
        method,
        "gc" = .C(
          "GibbsgConst",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "gs" = .C(
          "GibbsgSB",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "gu" = .C(
          "GibbsgUser",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "rc" = .C(
          "GibbsRobustConst",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "rs" = .C(
          "GibbsRobustSB",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "ru" = .C(
          "GibbsRobustUser",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "lc" = .C(
          "GibbsLiangConst",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "ls" = .C(
          "GibbsLiangSB",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "lu" = .C(
          "GibbsLiangUser",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "zc" = .C(
          "GibbsZSConst",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "zs" = .C(
          "GibbsZSSB",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "zu" = .C(
          "GibbsZSUser",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "fc" = .C(
          "GibbsflsConst",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "fs" = .C(
          "GibbsflsSB",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        ),
        "fu" = .C(
          "GibbsflsUser",
          as.character(""),
          as.integer(n),
          as.integer(p),
          as.integer(49),
          as.character(wd),
          as.integer(1),
          as.double(estim.time),
          as.integer(knull),
          as.integer(1),
          as.integer(seed)
        )
      )

      estim.time <- result[[7]] * (n.burnin + n.iter) / (60 * 50)
      cat("The problem would take ",
          ceiling(estim.time),
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

    #Call the corresponding function:
    result <- switch(
      method,
      "gc" = .C(
        "GibbsgConst",
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
      "gs" = .C(
        "GibbsgSB",
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
      "gu" = .C(
        "GibbsgUser",
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
      "rc" = .C(
        "GibbsRobustConst",
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
      "rs" = .C(
        "GibbsRobustSB",
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
      "ru" = .C(
        "GibbsRobustUser",
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
      "lc" = .C(
        "GibbsLiangConst",
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
      "ls" = .C(
        "GibbsLiangSB",
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
      "lu" = .C(
        "GibbsLiangUser",
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
      "zc" = .C(
        "GibbsZSConst",
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
      "zs" = .C(
        "GibbsZSSB",
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
      "zu" = .C(
        "GibbsZSUser",
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
      "fc" = .C(
        "GibbsflsConst",
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
      "fs" = .C(
        "GibbsflsSB",
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
      "fu" = .C(
        "GibbsflsUser",
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
      )
    )

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


    #Highest probability model
    mod.mat <- as.data.frame(t(models))


    inclusion <- incl
    names(inclusion) <- namesx
    result <- list()
    #
    result$time <- time #The time it took the programm to finish
    result$lmfull <- lmfull # The lm object for the full model
    if(!is.null(fixed.cov)){
      result$lmnull <- lmnull # The lm object for the null model
    }

    result$variables <- namesx #The name of the competing variables
    result$n <- n #number of observations
    result$p <- p #number of competing variables
    result$k <- knull#number of fixed covariates
    result$HPMbin <- models#The binary code for the HPM model
    names(result$HPMbin) <- namesx
    #result$modelsprob <- mod.mat
    result$modelslogBF <-modelslBF#The binary code for all the visited models (after n.thin is applied) and the correspondent log(BF)
    result$inclprob <- inclusion #inclusion probability for each variable
    names(result$inclprob) <- namesx

    result$jointinclprob <- data.frame(joint[1:p,1:p],row.names=namesx)#data.frame for the joint inclusion probabilities
    names(result$jointinclprob) <- namesx
    #
    result$postprobdim <- dimen #vector with the dimension probabilities.
    names(result$postprobdim) <- (0:p)+knull #dimension of the true model
    #
    #result$betahat <- betahat
    #rownames(result$betahat)<-namesx
    #names(result$betahat) <- "BetaHat"
    result$call <- match.call()
    if (pfms == "c" ) priorprobs <- rep(1, p + 1)
    if (pfms == "s" ) priorprobs <- 1/choose(p,0:p)
    result$priorprobs <- priorprobs
		
		result$C<- calculaC(modelslBF, priorprobs, p)
		
    result$method <- "gibbs"
    class(result)<- "Bvs"
    result


  }


#A function to obtain an estimate of the normalizing constant
#based on the method by George&McCulloch(1997)
#from the output of Gibbs
#not to be exported as is expected to be used only within the GibbsBvs
calculaC<- function(modelslBF, priorprobs, p){
	n<- dim(modelslBF)[1]
	#The method of George&McCulloch uses to sets. Here these are obtained as 
	#subsets of the MCMC sample of length K
	K<- round(n/2)
	Aset<- sample(x=1:n, size=K, rep=F)
	Bset<- (1:n)[-Aset]
	#log-Bayes factors of the models in A
	lBFAset<- modelslBF[Aset, "logBFi0"]
	#dimension of these models
	dimAset<- rowSums(modelslBF[Aset, -dim(modelslBF)[2]])

	#Remove repetitions to have finally the definition of A
	notdup<- !duplicated(lBFAset)
	lBFAset<- unique(lBFAset)
	dimAset<- dimAset[notdup]

	#Prior probabilities of the models in A (suming one is because the first position is occupied by dimension=0)
	priorAset<- priorprobs[dimAset+1]/sum(priorprobs*choose(p, 0:p))
	#The sum of Bi0*Pr(Mi) (g(A) in G&McC notation)
	gAset<- sum(exp(lBFAset+log(priorAset)))
	#How many of the models in Bset are in A?
	sumIA<- sum(modelslBF[Bset,"logBFi0"]%in%modelslBF[Aset,"logBFi0"])

	#The estimation of the normalizing constant is then: (1/C in G&McC's notation)
	C<- gAset*K/sumIA
	return(C)
}
