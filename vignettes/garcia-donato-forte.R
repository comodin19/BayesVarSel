#Accompanying code to the paper
#Bayesian Testing, Variable Selection and model averaging in linear Models using R with BayesVarSel
#by Gonzalo Garcia-Donato and Anabel Forte


library(BayesVarSel)


#################################
#Subsection Testing
#################################

###
#Example weight gains

#Two-samples t-test
#H0:mu1=mu2 vs. H1:mu1\ne mu2 (sigmas unknown but equal)

#(data taken from Lee (1997), page 143)
#samples of weight gains (in grams) between the 28th and 84th days of age rats
#recieving a high proteinic diet (diet=1) or not receiving such diet (diet=0)
weight.gains <- c(134, 146, 104, 119, 124, 161, 107, 83, 113, 129, 97, 123,
	70, 118, 101, 85, 107, 132, 94)
diet <- as.factor(c(rep(1,12),rep(0,7)))
rats <- data.frame(weight.gains=weight.gains,diet=diet)
M0 <- weight.gains ~ 1
M1 <- weight.gains ~ diet

Btest(models=c(H0=M0, H1=M1), data=rats)
# ---------
# Models:
# $H0
# y ~ 1
#
# $H1
# y ~ diet
#
# ---------
# Bayes factors (expressed in relation to H0)
#  H0.to.H0  H1.to.H0
# 1.0000000 0.8040127
# ---------
# Posterior probabilities:
#    H0    H1
# 0.554 0.446

###
#Example Savings

#Testing all predictors
#H0:beta1=beta2=...betap=0 vs. H1:no H0.
#example taken from Faraway(2002), Practical Regression and Anova using R (page 29)
#old economic dataset on 50 different countries. These data are averages over 1960-1970 (to remove
#business cycle or other short-term fluctuations). dpi is per-capita disposable income in U.S. dollars; ddpi is the percent rate of change in
#per capita disposable income; sr is aggregate personal saving divided by disposable income. The percentage population under 15 (pop15) and
#over 75 (pop75) are also recorded. The data come from Belsley, Kuh, and Welsch (1980).


#install.packages("faraway")
data("savings", package = "faraway")

fullmodel <- sr ~ pop15 + pop75 + dpi + ddpi
nullmodel <- sr ~ 1

Btest(models=c(H0=nullmodel, H1=fullmodel), data=savings)
# ---------
# Models:
# $H0
# sr ~ 1
#
# $H1
# sr ~ pop15 + pop75 + dpi + ddpi
#
# ---------
# Bayes factors (expressed in relation to H0)
# H0.to.H0 H1.to.H0
#   1.0000  21.46007
# ---------
# Posterior probabilities:
#    H0    H1
# 0.045 0.955

###
#Example savings (cont.)

#Multiple hypotheses
#H0, H1, H2, etc
#Continuation of the example above taken from Faraway(2002)
#H_0:beta1=beta2=...betap=0; H1:all betai\ne 0; H2:betai=0
fullmodel <- sr ~ pop15 + pop75 + dpi + ddpi
reducedmodel <- sr ~ pop75 + dpi + ddpi
nullmodel <- sr ~ 1
Btest(models=c(H0=nullmodel, H1=fullmodel, H2=reducedmodel), data=savings)
# ---------
# Models:
# $H0
# sr ~ 1
#
# $H1
# sr ~ pop15 + pop75 + dpi + ddpi
#
# $H2
# sr ~ pop75 + dpi + ddpi
#
# ---------
# Bayes factors (expressed in relation to H0)
#   H0.to.H0   H1.to.H0   H2.to.H0
#  1.0000000 21.4600656  0.7017864
# ---------
# Posterior probabilities:
#    H0    H1    H2
# 0.043 0.927 0.030

#################################
#Subsection Variable Selection
#################################

###
#Example savings (cont.)
Bvs(formula = sr ~ pop15 + pop75 + dpi + ddpi, data = savings)
# Info. . . .
# Most complex model has 5 covariates
# From those 1 is fixed and we should select from the remaining 4
# pop15, pop75, dpi, ddpi
# The problem has a total of 16 competing models
# Of these, the  10 most probable (a posteriori) are kept
# Working on the problem...please wait.
#
# Call:
# Bvs(formula = "sr~pop15+pop75+dpi+ddpi", data = savings)
#
# The 10 most probable models and their probabilities are:
#    pop15 pop75 dpi ddpi        prob
# 1      *     *   *    * 0.297315642
# 2      *     *        * 0.243433493
# 3      *              * 0.133832367
# 4      *                0.090960327
# 5      *         *    * 0.077913429
# 6      *     *          0.057674755
# 7      *         *      0.032516780
# 8      *     *   *      0.031337639
# 9                       0.013854369
# 10           *        * 0.006219812
#

###################################################
#Section Hypothesis testing with BayesVarSel
###################################################

###
#Example savings (cont.)
fullmodel <-  sr ~ pop15 + pop75 + dpi + ddpi
reducedmodel<- sr ~ pop75 + dpi + ddpi
nullmodel <- sr ~ 1
Btest(models=c(H0=nullmodel, H1=fullmodel, H2=reducedmodel), data=savings, prior.models="User",
            priorprobs=c(H0=1/2, H1=1/4, H2=1/4))
# ---------
# Models:
# $H0
# sr ~ 1
#
# $H1
# sr ~ pop15 + pop75 + dpi + ddpi
#
# $H2
# sr ~ pop75 + dpi + ddpi
#
# ---------
# Bayes factors (expressed in relation to H0)
#   H0.to.H0   H1.to.H0   H2.to.H0
#  1.0000000 21.4600656  0.7017864
# ---------
# Posterior probabilities:
#    H0    H1    H2
# 0.083 0.888 0.029

###
#Example 2 (testing a subspace)

#Faraway, page 32 (testing linear combinations)
#H0:betapop15=betapop75
#H1:full
fullmodel <- sr ~ pop15 + pop75 + dpi + ddpi
equalpopmodel <- sr ~ I(pop15 + pop75) + dpi + ddpi
#Btest(models=c(H0 = equalpopmodel, H1 = fullmodel), data = savings)
#This sentence produces the following error:
#Error in Btest(models = c(H0 = equalpopmodel, H1 = fullmodel), data = savings) :
#  I suspect that perhaps the simplest (null) model is not nested in all the others.
# Define explicitly the simplest model if you are sure it is the case.

Btest(models=c(Heqp = equalpopmodel, H1 = fullmodel), data = savings, null.model = "Heqp")
#---------
#Bayes factors (expressed in relation to H0)
# H0.to.H0  H1.to.H0
#1.0000000 0.3336251
#---------
#Posterior probabilities:
#  H0   H1
#0.75 0.25

###################################################
#Section Variable Selection with BayesVarSel
###################################################

###
#Example Crime

#USCrime changing the null model
crime.Edfix <- Bvs(formula= y ~ . , data = UScrime, null.model = y ~ Ed)
# Info. . . .
# Most complex model has 16 covariates
# From those 2 are fixed and we should select from the remaining 14
# M, So, Po1, Po2, LF, M.F, Pop, NW, U1, U2, GDP, Ineq, Prob, Time
# The problem has a total of 16384 competing models
# Of these, the  10 most probable (a posteriori) are kept
# Working on the problem...please wait.
crime.Edfix
# Call:
#   Bvs(formula = y ~ ., data = UScrime, null.model = y ~ Ed)
#
# The 10 most probable models and their probabilities are:
#   M So Po1 Po2 LF M.F Pop NW U1 U2 GDP Ineq Prob Time        prob
# 1         *                                *           0.059218556
# 2  *      *                                *           0.023262296
# 3         *                                *    *      0.022282369
# 4  *      *                                *    *      0.022182995
# 5  *      *                       *        *    *      0.019136036
# 6  *      *                       *        *           0.012113449
# 7             *                            *           0.011247404
# 8  *      *                           *    *           0.008213641
# 9  *  *   *   *  *   *   *  *  *  *   *    *    *    * 0.007750578
# 10        *          *                     *           0.007305878

###
#Example crime (cont.)
# USCrime changing the prior model probabilities (theta=1/4)

theta <- 1/4; pgamma<- 0:15
crime.thQ <- Bvs(formula= y ~ . , data = UScrime, prior.models = "User",
            priorprobs = theta^pgamma*(1-theta)^(15-pgamma))
# Info. . . .
# Most complex model has 16 covariates
# From those 1 is fixed and we should select from the remaining 15
# M, So, Ed, Po1, Po2, LF, M.F, Pop, NW, U1, U2, GDP, Ineq, Prob, Time
# The problem has a total of 32768 competing models
# Of these, the  10 most probable (a posteriori) are kept
# Working on the problem...please wait.
crime.thQ
# Call:
#   Bvs(formula = y ~ ., data = UScrime, prior.models = "User", priorprobs = theta^pgamma *
#         (1 - theta)^(15 - pgamma))
#
# The 10 most probable models and their probabilities are:
#   M So Ed Po1 Po2 LF M.F Pop NW U1 U2 GDP Ineq Prob Time       prob
# 1        *   *                                *           0.06773093
# 2  *     *   *                                *           0.03790366
# 3        *   *                                *    *      0.03631971
# 4  *     *   *                                *    *      0.03345791
# 5            *          *                     *           0.02128005
# 6  *     *   *                       *        *    *      0.01891988
# 7  *     *   *                       *        *           0.01838720
# 8            *          *                     *    *      0.01671574
# 9        *       *                            *           0.01296033
# 10 *     *   *                           *    *           0.01251800

###
#Example SDM data with p=67
data(SDM)
set.seed(1234)
#Here the number of expected regressors, a priori, is 7 (as in Ley and Steel 09)
p <- ncol(SDM) - 1;wstar <- 7; b <- 	(67-wstar)/wstar; pgamma<- 0:67

growth.wstar7 <- GibbsBvs(formula= y ~ . , data = SDM, prior.models = "User",
            priorprobs = gamma(pgamma+1)*gamma(67-pgamma+b), n.iter = 10000, n.thin = 1, time.test = FALSE)
growth.wstar7
# Most complex model has 68 covariates
# From those 1 is fixed and we should select from the remaining 67
# ABSLATIT, AIRDIST, AVELF, BRIT, BUDDHA, CATH00, CIV72, COLONY, CONFUC, DENS60, DENS65C, DENS65I, DPOP6090, EAST, ECORG, ENGFRAC, EUROPE, FERTLDC1, GDE1, GDPCH60L, GEEREC1, GGCFD3, GOVNOM1, GOVSH61, GVR61, H60, HERF00, HINDU00, IPRICE1, LAAM, LANDAREA, LANDLOCK, LHCPC, LIFE060, LT100CR, MALFAL66, MINING, MUSLIM00, NEWSTATE, OIL, OPENDEC1, ORTH00, OTHFRAC, P60, PI6090, SQPI6090, PRIGHTS, POP1560, POP60, POP6560, PRIEXP70, PROT00, RERD, REVCOUP, SAFRICA, SCOUT, SIZE60, SOCIALIST, SPAIN, TOT1DEC1, TOTIND, TROPICAR, TROPPOP, WARTIME, WARTORN, YRSOPEN, ZTROPICS
# The problem has a total of 1.47574e+20 competing models
# Of these, 10500 are sampled with replacement
# Then, 10000 are kept and used to construct the summaries
# Working on the problem...please wait.



#############################################################
#Section Summaries of the posterior distribution
#############################################################

#Example crime (cont.)
summary(crime.Edfix)
# Call:
# Bvs(formula = y ~ ., data = UScrime, null.model = y ~ Ed)
#
# Inclusion Probabilities:
#      Incl.prob. HPM MPM
# M        0.6806       *
# So       0.2386
# Po1      0.8489   *   *
# Po2      0.3663
# LF       0.2209
# M.F      0.3184
# Pop      0.2652
# NW       0.2268
# U1       0.2935
# U2       0.4765
# GDP      0.3204
# Ineq     0.9924   *   *
# Prob     0.6174       *
# Time     0.2434
# ---
# Code: HPM stands for Highest posterior Probability Model and
#  MPM for Median Probability Model.

###
#plots of crime.Edfix

#pdf(file = "../plots/crimeEdfixjoint.pdf",width = 8, height=8)
mj <- plot(crime.Edfix, option = "joint")
#dev.off()

#pdf(file = "../plots/crimeEdfixconditional.pdf",width = 8, height=8)
mc <- plot(crime.Edfix, option = "conditional")
#dev.off()

#pdf(file = "../plots/crimeEdfixnot.pdf",width = 8, height=8)
mn <- plot(crime.Edfix, option = "not")
#dev.off()


###
#The probability of Po2 given not Po1 is
mn["Not.Po1" , "Po2"]
#[1] 0.9996444

###
#Jointness measures
joint_meassures <- Jointness(crime.Edfix, covariates = c("Po1" , "Po2"))
joint_meassures
#---------
#The joint inclusion probability for Po1 and Po2 is:  0.22
#---------
#The ratio between the probability of including both covariates and the probability of including at least one of then is: 0.22
#---------
#The probability of including both covariates together is 0.27 times the probability of including one of them alone


#The single number can be accessed using
joint_meassures$prob_joint
joint_meassures$joint_LS1
joint_meassures$joint_LS2


###
#dimension plot:
#pdf(file = "../plots/crimeEdfixdimension.pdf",width = 8, height=8)
plot(crime.Edfix, option = "dimension")
#dev.off()


#trace plot:
#pdf(file = "../plot/growthtrace.pdf",width = 10, height=8)
plot(growth.wstar7, option = "trace")
#dev.off()

###########################################################
#Section Model averaged estimations and predictions / Estimation
###########################################################

###
#Example Crime
#bma estimations: first run again the code for UScrime increasing the number
#of models saved to perform more accurate estimations
crime.Edfix <- Bvs(formula= y ~ . , data = UScrime, null.model = y ~ Ed, n.keep = 2000)
crime.Edfix
#now run the bma command:
set.seed(1234)
bma.crime.Edfix <- BMAcoeff(crime.Edfix)

###
#plots of the resulting posterior distributions:
#pdf(file = "../plots/Ineq.pdf",width = 8, height=6)
histBMA(bma.crime.Edfix, covariate = "Ineq", n.breaks = 50)
#dev.off()

#pdf(file = "../plots/Time.pdf",width = 8, height=6)
histBMA(bma.crime.Edfix, covariate = "Time", n.breaks = 50)
#dev.off()

#pdf(file = "../plots/Prob.pdf",width = 8, height=6)
histBMA(bma.crime.Edfix, covariate = "Prob", n.breaks = 50)
#dev.off()

###
#Example 3
#describing with percentiles the posterior distribution
quantile(bma.crime.Edfix[ , "Ineq"], probs = c(0.05, 0.5, 0.95))
#      5%       50%       95%
#4.075685  7.150184 10.326606


###
#Example SDM
#bma estimations in the SDM case:
set.seed(1234)
bma.growth.wstar7<- BMAcoeff(growth.wstar7)
histBMA(bma.growth.wstar7, covariate = "P60", n.breaks = 50)

###
#what is the probability that the effect of P60 over savings is greater than one?
mean(bma.growth.wstar7[ , "P60"] > 1)
#[1] 0.7511


###########################################################
#Section Model averaged estimations and predictions / Prediction
###########################################################

###
#Let us predict the y's associated with the "mean" case in the SDM example
set.seed(1234)
pred.growth.wstar7 <- predict(object = growth.wstar7, newdata = data.frame(t(colMeans(SDM))))
#pdf(file = "../plots/PredAverage.pdf",width = 8, height=6)
hist(pred.growth.wstar7[ , 1], main = "SDM", border = gray(0.6), col = gray(0.8), xlab = "y")
#dev.off()
