#' @export
GibbsBvsF <-
  function(formula,
           data,
           null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),					 
           prior.betas = "Robust",
           prior.models = "SBSB2",
           n.iter = 10000,
           init.model = "Full",
           n.burnin = 500,
           n.thin = 1,
           time.test = TRUE,
           priorprobs = NULL,
           seed = runif(1, 0, 16091956)) {

		cat("Recall: at this point of development null model should contain intercept\n")
		
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
			#Factors:
			#before:X.full <- lmfull$x
    	X.full<- lmerTest:::get_rdX(lmfull) #rank defficient paramet

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

		#Factors:
		#positions is a matrix with number of rows equal to the number of regressors
		#(either factor or numeric) and number of columns the number of columns of X
		#Each row describes the position (0-1) in X of a regressor (several positions in case
		#this regressor is a factor)
		depvars<- attr(lmfull$terms, "term.labels")
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
		if (prior.models!="SBSB2" & prior.models!="ConstConst2" & prior.models!="SB2" & prior.models!="Const2" & prior.models!="SBConst2")
			{stop("Prior over the model space not supported\n")}
		
		
    if (prior.betas == "Unitary"){write(0, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}
    if (prior.betas == "Robust"){write(1, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}
    if (prior.betas == "Liangetal"){write(4, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}
    if (prior.betas == "gZellner"){write(2, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}
    if (prior.betas == "ZellnerSiow"){write(5, ncolumns=1, file=paste(wd, "/typeofBF.txt", sep = ""))}
	  if (prior.betas == "FLS"){stop("Prior FLS not yet supported\n")}

		if (prior.betas != "Unitary" & prior.betas != "Robust" & prior.betas != "Liangetal" & 
		    prior.betas != "gZellner" & prior.betas != "ZellnerSiow" & prior.betas != "FLS") {stop("Dont recognize the prior for betas\n")}	
		
		if (prior.models=="SBSB2"){method<- "rSBSB"}
		if (prior.models=="ConstConst2"){method<- "rConstConst"}
		if (prior.models=="SBConst2"){method<- "rSBConst"}	
		if (prior.models=="SB2"){method<- "rSB"}
		if (prior.models=="Const2"){method<- "rConst"}
		
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
		if (prior.models == "SBSB2"){modelslBFwR<- BayesVarSel:::resamplingSBSB(modelslBF, positions)}
		if (prior.models == "ConstConst2"){modelslBFwR<- BayesVarSel:::resamplingConstConst(modelslBF, positions)}
		if (prior.models == "SB2"){modelslBFwR<- BayesVarSel:::resamplingSB(modelslBF, positions)}
		if (prior.models == "Const2"){modelslBFwR<- BayesVarSel:::resamplingConst(modelslBF, positions)}
		if (prior.models == "SBConst2"){modelslBFwR<- BayesVarSel:::resamplingSBConst(modelslBF, positions)}
		
		#Now we convert the sampled models to vars and factors:
		modelslBFF<- t(apply(modelslBFwR[,-(p+1)], MARGIN=1, FUN=function(x,M){as.numeric(x%*%M>0)}, M=t(positions)))
		
    #Inclusion probabilities:
		inclusion <- colMeans(modelslBFF)
		names(inclusion)<- depvars
		
		#HPM with factors:
		HPMFbin<- modelslBFF[which(allBF==max(allBF))[1], ]
		names(HPMFbin)<- depvars
		
		#joint inclusion probs:
		jointinclprob<- matrix(0, ncol=dim(positions)[1], nrow=dim(positions)[1])

		for (i in 1:dim(modelslBFF)[1]){
			jointinclprob<- jointinclprob + matrix(modelslBFF[i,], nc=1)%*%matrix(modelslBFF[i,], nr=1)
		}
			jointinclprob<- jointinclprob/dim(modelslBFF)[1]
		colnames(jointinclprob)<- depvars; rownames(jointinclprob)<- depvars
		
		#Dimension of the true model:
		dimenF<- c(rowSums(modelslBFF), 0:dim(positions)[1])
		dimenF<- (table(dimenF)-1)/dim(modelslBFF)[1]
		names(dimenF)<- (0:dim(positions)[1])+knull
		
		#Attach the column with the log(BF)
		modelslBFF<- cbind(modelslBFF, log(allBF))
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
		
    result$inclprob <- inclusion #inclusion probability for each variable

    result$jointinclprob <- data.frame(jointinclprob) #data.frame for the joint inclusion probabilities
    #
    result$postprobdim <- dimenF #vector with the dimension probabilities.
		
		result$positions<- positions
		result$positionsx<- positionsx
		
		#
		#
		#####################
		

    result$priorprobs <- "En obras"
    result$method <- "gibbsWithFactors"
    class(result)<- "Bvs"
    result


  }
	

