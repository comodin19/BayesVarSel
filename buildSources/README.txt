#12-6-14 Note: this version puts all sources files needed in
#BayesVarSel/src in ./sources
#Only changes main.c and mainGibbs.c
#If priorprob.{c,h}, allBF.{c,h} auxiliaryfuncs.{c,h} or Gibbsauxiliaryfuncs.{c,h}
#need to be changes this should be performed by hand.
#Note in this regard that auxiliaryfuncs has different functions for priorprob=SB
#or priorprob=Constant
#This version contemplates: priorprob=ScottBerger and Constant and
#Bayes Factors=gprior, Robust, Liang and ZS 
#Inputs are for main.c: maiHeader.c and mainTemplateVnew.c (main body) where Vnew is the version
#and for mainGibbs.c: mainGibbsHeader.c and mainGibbsTemplateVnew.c (main body) 
#The file templates sould be versions with gConst (ie gBayesFactor and constant priorprobs)


