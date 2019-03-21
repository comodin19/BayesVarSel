resamplingSBSB<- function(modelslBF, positions){
	#Given the matrix created within BvsF containing a binary expression of models and the log(BF) in the last columnbtained 
	#with GibbsBvs(...prior.models = "SBSB"...)
	#(meaning that prior over the model space is obtained with the hierarchical
	#Pr(Mg)=Pr(Mg|those factors and x's)*Pr(those factors and x's)
	#and Pr(Mg|those factors and x's)\proto 1/(number of models of that dimension with those factors and x's)
	#and Pr(those factors and x's)\propto 1/(number of models of that dimension))
	#the function makes a resampling with the new priolevelsfullr being
	#Pr(Mg|those factors and x's)\propto 1/(number of different models of that rank with those factors and x's)
	#returning a matrix of the same size with the resampled models
	levelsfull<- rowSums(positions)
	levelsf<- levelsfull[levelsfull>1]

	kplusp<- dim(positions)[1]
	#A matrix containing, for each sampled models, the number of active "levels" for each regressor:
	if (kplusp==1)
		m.actlevels<- matrix(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}), nc=1)
	else
		m.actlevels<- t(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}))
	#for each sampled model obtain the SBSB2 prior prob:
	allmodelsprior2<- t(apply(m.actlevels, MARGIN=1, FUN=priorSBSB2, levelsfull=levelsfull, levelsf=levelsf, kplusp=kplusp))
	#for each sampled model obtain the SBSB1 prior prob:
	allmodelsprior1<- t(apply(m.actlevels, MARGIN=1, FUN=priorSBSB1, levelsfull=levelsfull, kplusp=kplusp))
	#now the resampling:
	resamp<- sample(x=1:dim(modelslBF)[1], size=dim(modelslBF)[1], rep=T, prob=exp(allmodelsprior2-allmodelsprior1))
	return(modelslBF[resamp,])	
}

resamplingConstConst<- function(modelslBF, positions){
	#Given the matrix created within BvsF containing a binary expression of models and the log(BF) in the last column obtained 
	#with GibbsBvs(...prior.models = "ConstConst"...)
	#(meaning that prior over the model space is obtained with the hierarchical
	#Pr(Mg)=Pr(Mg|those factors and x's)*Pr(those factors and x's)
	#and Pr(Mg|those factors and x's)\proto 1/(number of models with those factors and x's)
	#and Pr(those factors and x's)\propto 1/(number of models))
	#the function makes a resampling with the new prior being the same but only
	#keeping unique models (in a same class the full is kept and the others are not)
	#returning a matrix of the same size with the resampled models
	levelsfull<- rowSums(positions)
	levelsf<- levelsfull[levelsfull>1]
	kplusp<- dim(positions)[1]
	#A matrix containing, for each sampled models, the number of active "levels" for each regressor:
	if (kplusp==1)
		m.actlevels<- matrix(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}), nc=1)
	else
		m.actlevels<- t(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}))
	#for each sampled model obtain the ConstConst2 prior prob:
	allmodelsprior2<- t(apply(m.actlevels, MARGIN=1, FUN=priorConstConst2, levelsfull=levelsfull, levelsf=levelsf, kplusp=kplusp))
	#for each sampled model obtain the SBSB1 prior prob:
	allmodelsprior1<- t(apply(m.actlevels, MARGIN=1, FUN=priorConstConst1, kplusp=kplusp, levelsfull=levelsfull))
	#now the resampling:
	resamp<- sample(x=1:dim(modelslBF)[1], size=dim(modelslBF)[1], rep=T, prob=exp(allmodelsprior2-allmodelsprior1))
	return(modelslBF[resamp,])	
}


resamplingSBConst<- function(modelslBF, positions){
	#Given the matrix created within BvsF containing a binary expression of models and the log(BF) in the last column obtained 
	#with GibbsBvs(...prior.models = "SBConst"...)
	#(meaning that prior over the model space is obtained with the hierarchical
	#Pr(Mg)=Pr(Mg|those factors and x's)*Pr(those factors and x's)
	#and Pr(Mg|those factors and x's)\proto 1/(number of models with those factors and x's)
	#and Pr(those factors and x's)\propto 1/((number of models of that dimension))
	#the function makes a resampling with the new prior being the same but only
	#keeping unique models (in a same class the full is kept and the others are not)
	#returning a matrix of the same size with the resampled models
	levelsfull<- rowSums(positions)
	levelsf<- levelsfull[levelsfull>1]
	kplusp<- dim(positions)[1]
	#A matrix containing, for each sampled models, the number of active "levels" for each regressor:
	if (kplusp==1)
		m.actlevels<- matrix(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}), nc=1)
	else
		m.actlevels<- t(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}))
	#for each sampled model obtain the ConstConst2 prior prob:
	allmodelsprior2<- t(apply(m.actlevels, MARGIN=1, FUN=priorSBConst2, levelsfull=levelsfull, levelsf=levelsf, kplusp=kplusp))
	#for each sampled model obtain the SBSB1 prior prob:
	allmodelsprior1<- t(apply(m.actlevels, MARGIN=1, FUN=priorSBConst1, kplusp=kplusp, levelsfull=levelsfull, levelsf=levelsf))
	#now the resampling:
	resamp<- sample(x=1:dim(modelslBF)[1], size=dim(modelslBF)[1], rep=T, prob=exp(allmodelsprior2-allmodelsprior1))
	return(modelslBF[resamp,])	
}


resamplingSB<- function(modelslBF, positions){
	#Given the matrix created within BvsF containing a binary expression of models and the log(BF) in the last columnbtained 
	#with GibbsBvs(...prior.models = "SB"...)
	#(meaning that prior over the model space is obtained with the standard SB prior
	#Pr(Mg)=\proto 1/(number of models of that dimension)
	#the function makes a resampling with the new prior being
	#Pr(Mg)\propto 1/(number of different models of that rank)
	#returning a matrix of the same size with the resampled models
	levelsfull<- rowSums(positions)
	levelsf<- levelsfull[levelsfull>1]
	#ff contains the matrix with the count of models of each rank for each combination of levels
	ff<- rank.levels2(pmax(levelsfull,2))
	kplusp<- dim(positions)[1]
	#A matrix containing, for each sampled models, the number of active "levels" for each regressor:
	if (kplusp==1)
		m.actlevels<- matrix(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}), nc=1)
	else
		m.actlevels<- t(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}))
	#for each sampled model obtain the SBSB2 prior prob:
	allmodelsprior2<- t(apply(m.actlevels, MARGIN=1, FUN=priorSB2, levelsfull=levelsfull, levelsf=levelsf, ff=ff))
	#for each sampled model obtain the SBSB1 prior prob:
	allmodelsprior1<- t(apply(m.actlevels, MARGIN=1, FUN=priorSB1, levelsfull=levelsfull))
	#now the resampling:
	resamp<- sample(x=1:dim(modelslBF)[1], size=dim(modelslBF)[1], rep=T, prob=exp(allmodelsprior2-allmodelsprior1))
	return(modelslBF[resamp,])	
}

resamplingConst<- function(modelslBF, positions){
	#Given the matrix created within BvsF containing a binary expression of models and the log(BF) in the last columnbtained 
	#with GibbsBvs(...prior.models = "Const"...)
	#(meaning that prior over the model space is obtained with standard constant prior
	#Pr(Mg)\proto 1/(number of models)
	#the function makes a resampling with the new prior being
	#Pr(Mg)\propto 1/(number of unique models)
	#returning a matrix of the same size with the resampled models
	levelsfull<- rowSums(positions)
	levelsf<- levelsfull[levelsfull>1]
	#ff contains the matrix with the count of models of each rank for each combination of levels
	ff<- matrix.rank.levels(all.levelsf=levelsf)
	kplusp<- dim(positions)[1]
	#A matrix containing, for each sampled models, the number of active "levels" for each regressor:
	if (kplusp==1)
		m.actlevels<- matrix(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}), nc=1)
	else
		m.actlevels<- t(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}))
	#for each sampled model obtain the SBSB2 prior prob:
	allmodelsprior2<- t(apply(m.actlevels, MARGIN=1, FUN=priorConst2, levelsfull=levelsfull, levelsf=levelsf))
	#for each sampled model obtain the SBSB1 prior prob:
	allmodelsprior1<- t(apply(m.actlevels, MARGIN=1, FUN=priorConst1, levelsfull=levelsfull))
	#now the resampling:
	resamp<- sample(x=1:dim(modelslBF)[1], size=dim(modelslBF)[1], rep=T, prob=exp(allmodelsprior2-allmodelsprior1))
	return(modelslBF[resamp,])	
}


priorSBSB1<- function(act.levels, levelsfull, kplusp){
	#The original SBSB prior
	lprMgamma<- -sum(log(levelsfull[act.levels>0])+lchoose(levelsfull[act.levels>0], act.levels[act.levels>0]))-log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0))
	return(lprMgamma)
}
	
priorSBSB2<- function(act.levels, levelsfull, levelsf, kplusp){
	#The corrected SB-SB prior prob (conditionally: inversely proportional to the number of models of that rank)
	#For copies of the same model, we only keep one representative (the full on that class)
	#(four different cases)
	
	#of the active levels take only those that correspond to factors:	
	act.levelsf<- act.levels[levelsfull>1]	
	
	#If the model does not contain any factor:
	if (sum(act.levelsf)==0) return(-log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0)))
	
	#Obtain the vector with the number of models of each rank:
	numberof<- rank.levels(levelsfull[levelsfull>1][act.levelsf>0])
	
	#if the model is not saturated nor oversaturated
	if (sum(act.levelsf >= (levelsf-1)) == 0){
		aux<- sum(act.levelsf)
		l2prMgamma<- -log(numberof[as.character(aux)])-
								log(sum(levelsf[act.levelsf>0])-2*length(levelsf[act.levelsf]>0)+1)-
								log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0))	
		return(l2prMgamma)	
			  }
  #if the model contains at least one level saturated then return -Inf:
	if (sum(act.levelsf[act.levelsf>0] == (levelsf[act.levelsf>0]-1)) >= 1){
		return(-Inf)
	}
	else { #keep the rest
	aux<- sum(act.levelsf[act.levelsf==levelsf]-1)+sum(act.levelsf[act.levelsf<levelsf])
	l2prMgamma<- -log(numberof[as.character(aux)])-
	              log(sum(levelsf[act.levelsf>0])-2*length(levelsf[act.levelsf>0])+1)-log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0))
	return(l2prMgamma)	
  }
}


priorConstConst1<- function(act.levels, kplusp, levelsfull){
	#The original Const-Const prior (conditionally, inversely proportional to the number of models)
	
	#of the active levels take only those that correspond to factors:	
	act.levelsf<- act.levels[levelsfull>1]	
	
	lprMgamma<- -sum(log(2^levelsfull[act.levelsf>0]-1))-kplusp*log(2)
	return(lprMgamma)
}


priorConstConst2<- function(act.levels, levelsfull, levelsf, kplusp){
	#The corrected Const-Const prior prob (conditionally: proportional to a constant for unique models)
	#For copies of the same model, we only keep one representative (the full on that class)
	#(four different cases)
	#of the active levels take only those that correspond to factors:	
	act.levelsf<- act.levels[levelsfull>1]	
	
	#If the model does not contain any factor:
	if (sum(act.levelsf)==0) return(-kplusp*log(2))
	
	#if the model is not saturated nor oversaturated
	if (sum(act.levelsf >= (levelsf-1)) == 0){
		l2prMgamma<- -sum(log(2^levelsf[act.levelsf>0]-levelsf[act.levelsf>0]-1))-kplusp*log(2)
		return(l2prMgamma)	
			  }
  #if the model contains at least one level saturated then return -Inf:
	if (sum(act.levelsf[act.levelsf>0] == (levelsf[act.levelsf>0]-1)) >= 1){
		return(-Inf)
	}
	else { #keep the rest
		l2prMgamma<- -sum(log(2^levelsf[act.levelsf>0]-levelsf[act.levelsf>0]-1))-kplusp*log(2)
		return(l2prMgamma)	
  }
}
	

priorSB1<- function(act.levels, levelsfull){
	#The standard SB prior where prob over models is inversely proportional to the number of models
	#of that dimension
	dimMg<- sum(act.levels)
	return(-lchoose(sum(levelsfull), dimMg)-log(sum(levelsfull)+1))
}
	
	
priorSB2<- function(act.levels, levelsfull, levelsf, ff){
	#The corrected SB prior prob (inversely proportional to the number of models of that rank)
	#For copies of the same model, we only keep one representative (the full on that class)
	#(four different cases)

  #this prior is only based on the rank of the model:
  rankMg<- sum(pmin(act.levels[levelsfull>1], levelsfull[levelsfull>1]-1))+sum(act.levels[levelsfull==1])
	
	#for the null
	if (sum(rankMg==0)==length(rankMg)) return(-log(length(ff)))
		
	act.levelsf<- act.levels[levelsfull>1]	

	#if the model is not saturated nor oversaturated
	if (sum(act.levelsf >= (levelsf-1)) == 0){
		l2prMgamma<- -log(ff[as.numeric(names(ff))==rankMg])-log(length(ff))
		return(l2prMgamma)	
			  }
  #if the model contains at least one level saturated then return -Inf:
	if (sum(act.levelsf[act.levelsf>0] == (levelsf[act.levelsf>0]-1)) >= 1){
		return(-Inf)
	}
	else { #keep the rest
	l2prMgamma<- -log(ff[as.numeric(names(ff))==rankMg])-log(length(ff))
	return(l2prMgamma)	
  }
}

priorConst1<- function(act.levels, levelsfull){
	#The standard constant prior, inversely proportional to the number of models
	return(-sum(levelsfull)*log(2))
}

priorConst2<- function(act.levels, levelsfull, levelsf){
	#The constant prior over the unique models
	#For copies of the same model, we only keep one representative (the full on that class)
  #this prior is only based on the rank of the model:
  rankMg<- sum(pmin(act.levels[levelsfull>1], levelsfull[levelsfull>1]-1))+sum(act.levels[levelsfull==1])
	
	#Rui's number (number of unique models):
	lrui.number<- -sum(log(2^levelsfull-levelsfull))
	
	#for the null
	if (sum(rankMg==0)==length(rankMg)) return(lrui.number)
		
	act.levelsf<- act.levels[levelsfull>1]	

	#if the model is not saturated nor oversaturated
	if (sum(act.levelsf >= (levelsf-1)) == 0){
		return(lrui.number)	
			  }
  #if the model contains at least one level saturated then return -Inf:
	if (sum(act.levelsf[act.levelsf>0] == (levelsf[act.levelsf>0]-1)) >= 1){
		return(-Inf)
	}
	else { #keep the rest
	return(lrui.number)	
  }
}

priorSBConst1<- function(act.levels, levelsfull, kplusp, levelsf){
	#The original SBConst prior
	act.levelsf<- act.levels[levelsfull>1]	
	
	lprMgamma<- -sum(log(2^levelsf[act.levelsf>0]-1))-kplusp*log(2)-log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0))
	return(lprMgamma)
}
	

priorSBConst2<- function(act.levels, levelsfull, levelsf, kplusp){
	#The corrected SB-Const prior prob (conditionally: proportional to a constant for unique models)
	#For copies of the same model, we only keep one representative (the full on that class)
	#(four different cases)
	#of the active levels take only those that correspond to factors:	
	act.levelsf<- act.levels[levelsfull>1]	

	#If the model does not contain any factor:
	if (sum(act.levelsf)==0) return(-log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0)))
	
	#if the model is not saturated nor oversaturated
	if (sum(act.levelsf >= (levelsf-1)) == 0){
		l2prMgamma<- -sum(log(2^levelsf[act.levelsf>0]-levelsf[act.levelsf>0]-1))-log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0))
		return(l2prMgamma)	
			  }
  #if the model contains at least one level saturated then return -Inf:
	if (sum(act.levelsf[act.levelsf>0] == (levelsf[act.levelsf>0]-1)) >= 1){
		return(-Inf)
	}
	else { #keep the rest
		l2prMgamma<- -sum(log(2^levelsf[act.levelsf>0]-levelsf[act.levelsf>0]-1))-log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0))
		return(l2prMgamma)	
  }
}
	
my.choose<- function(n, k){
	if (k>n | k<0 | n<=0) stop("Bad arguments on my.choose\n")
	if (k<(n-1)) return(choose(n, k)) 
		else return (1) 
}

matrix.rank.levels<- function(all.levelsf){
	#Given a vector with all possible levels (l1,l2,...,lR), this function returns a matrix that
	#for each combination of levels contains how many
	#models there are of all possible ranks R<=r<=sum(l_i)-R
	#Example:
	#ff<- matrix.rank.levels(all.levelsf=c(3,2,4))
	#ff["f2f3", 2] will return how many models are of rank 2 when the factors 2 and 3 are active (each
		#with 2 and 4 levels respectively)	
	#
	#
	mmm<- matrix(0, nr=2^length(all.levelsf)-1, ncol=sum(all.levelsf)-length(all.levelsf))
	rownames(mmm)<- 1:(2^length(all.levelsf)-1)
	for (i in 1:(2^length(all.levelsf)-1)){
		whichactive<- integer.base.b_C(i, length(all.levelsf))==1
		resrank<- rank.levels(all.levelsf[whichactive])
		mmm[i,as.numeric(names(resrank))]<- resrank
		rownames(mmm)[i]<- paste("f",(1:length(all.levelsf))[whichactive], collapse="", sep="")		
	}
	return(mmm)	
}

rank.levels<- function(levelsf){
	#Expected to be used in combination with SBSB
	#Given a vector of ACTIVE factors with levels (l1,l2,...,lp), l_i>=2, this function computes how many
	#models --not saturated, so none of the active levels=l-1-- there are with the same number of active levels (r) 
	#p<=r<=sum(l_i)-p
	if (length(levelsf)==1){
		if (levelsf==2){
			result<- c(1)
			names(result)<- "1"
			return(result)
			}
	}
		
	p.levels<- prod(levelsf-1); n.levels<- length(levelsf)
	mm<- matrix(0, nrow=p.levels, ncol=n.levels)
	colnames(mm)<- paste("F", 1:length(levelsf), sep="")
	if (n.levels>1){
		for (i in 1:(n.levels-1)){
			mm[,i]<- rep(1:(levelsf[i]-1), each=prod(levelsf[(i+1):n.levels]-1), length.out=p.levels)
		}
			i<- i+1
			mm[,i]<- rep(1:(levelsf[i]-1), each=1, length.out=p.levels)
		}
	else mm[,1]<- rep(1:(levelsf[1]-1), each=1, length.out=p.levels)
	
	mm<- cbind(mm, rowSums(mm))
	mm<- cbind(mm, 0)
	colnames(mm)[length(levelsf)+1:2]<- c("sum.act.levels", "combin.prod")
	s<- 1
	for (i in 1:p.levels){
		for (j in 1:n.levels){
			s<- s*my.choose(levelsf[j], mm[i,j])
		}
		mm[i, n.levels+2]<- s; s<- 1
	}

	possible.values<- n.levels:(sum(levelsf)-n.levels)
	result<- rep(0, length(possible.values))
	names(result)<- possible.values
	for (i in possible.values){
		these<- which(mm[,"sum.act.levels"]==i)
		result[as.character(i)]<- sum(mm[these,"combin.prod"])		
	}
	#rectify for the saturated or oversaturated:
	result[as.character(i)]<- 1
	return(result)	
}


rank.levels2<- function(levelsfull){
	#Expected to be used in combination with SB functions
	#Given a vector of p factors and k variables levelsfull=(l1,l2,...,lp,1,...,1), 
	#(for variables we have the ones), this function computes how many
	#models there are of all possible ranks 0<=r<=sum(l_i-1)+k
	#(Because of a trick This works even if some of the levelsfull=1 (so numerical vars)
	#but I created the rank.levels3 exactly for that purpose
	#but it is much more complex and it does the same)
	
	#the resulting table is the same independently of li=1 or li=2
	levelsfull[levelsfull==1]<- 2

	if (length(levelsfull)==1){
		if (levelsfull==2){
			result<- c(1,1)
			names(result)<- c("0","1")
			return(result)
			}
	}
	p.levels<- prod(levelsfull); n.levels<- length(levelsfull)
	mm<- matrix(0, nrow=p.levels, ncol=n.levels)
	colnames(mm)<- paste("F", 1:length(levelsfull), sep="")
	if (n.levels>1){
		for (i in 1:(n.levels-1)){
			mm[,i]<- rep(0:(levelsfull[i]-1), each=prod(levelsfull[(i+1):n.levels]), length.out=p.levels)
		}
			i<- i+1
			mm[,i]<- rep(0:(levelsfull[i]-1), each=1, length.out=p.levels)
		}
	else mm[,1]<- rep(0:(levelsfull[1]-1), each=1, length.out=p.levels)
	
	mm<- cbind(mm, rowSums(mm))
	mm<- cbind(mm, 0)
	colnames(mm)[length(levelsfull)+1:2]<- c("sum.act.levels", "combin.prod")
	s<- 1
	for (i in 1:p.levels){
		for (j in 1:n.levels){
			s<- s*my.choose(levelsfull[j], mm[i,j])
		}
		mm[i, n.levels+2]<- s; s<- 1
	}

	possible.values<- min(mm[,"sum.act.levels"]):max(mm[,"sum.act.levels"])
	result<- rep(0, length(possible.values))
	names(result)<- possible.values
	for (i in possible.values){
		these<- which(mm[,"sum.act.levels"]==i)
		result[as.character(i)]<- sum(mm[these,"combin.prod"])		
	}
	#rectify for the saturated or oversaturated:
	result[as.character(i)]<- 1
	return(result)	
}

rank.levels3<- function(levelsfull){
	#Given a vector of ACTIVE factors and variables (l1,l2,...,lR), this function computes how many
	#models there are of all possible ranks 0<=r<=sum(l_i)-R
	
	#Same as rank.levels2 but for possible l_=1 (so numerical vars)
	#...at the end rank.levels2 also works in this case because of a trick there
	
	if (length(levelsfull)==1){
		if (levelsfull==2){
			result<- c(1)
			names(result)<- "1"
			return(result)
			}
	}
	
	#first the matrix for the factors:
	levelsfullF<- levelsfull[levelsfull>1]
	
	p.levels<- prod(levelsfullF); n.levels<- length(levelsfullF)
	mm<- matrix(0, nrow=p.levels, ncol=n.levels)
	colnames(mm)<- paste("F", 1:length(levelsfullF), sep="")
	if (n.levels>1){
		for (i in 1:(n.levels-1)){
			mm[,i]<- rep(0:(levelsfullF[i]-1), each=prod(levelsfullF[(i+1):n.levels]), length.out=p.levels)
		}
			i<- i+1
			mm[,i]<- rep(0:(levelsfullF[i]-1), each=1, length.out=p.levels)
		}
	else mm[,1]<- rep(0:(levelsfullF[1]-1), each=1, length.out=p.levels)
	
	#Second the numeric variables:
	k<- length(levelsfull[levelsfull==1])
	mm0<- mm
	if (k>0){
		mm<- cbind(mm0, rep(0,each=dim(mm0)[1]))
		for (i in 1:k){
			mm<- rbind(mm, cbind(mm0, rep(i, each=dim(mm0)[1])))
		}
	}
		
	mm<- cbind(mm, rowSums(mm))
	mm<- cbind(mm, 0)
	colnames(mm)[length(levelsfullF)+as.numeric(k>0)+1:2]<- c("sum.act.levels", "combin.prod")
	
	s<- 1
	for (i in 1:dim(mm)[1]){
		for (j in 1:n.levels){
			s<- s*my.choose(levelsfullF[j], mm[i,j])
		}
		if (k>0) s<- s*choose(k, mm[i, (n.levels+1)])
		mm[i, n.levels+3]<- s; s<- 1
	}

	possible.values<- min(mm[,"sum.act.levels"]):max(mm[,"sum.act.levels"])
	result<- rep(0, length(possible.values))
	names(result)<- possible.values
	for (i in possible.values){
		these<- which(mm[,"sum.act.levels"]==i)
		result[as.character(i)]<- sum(mm[these,"combin.prod"])		
	}
	#rectify for the saturated or oversaturated:
	result[as.character(i)]<- 1
	return(result)	
}



#Taken from Bvs.R:
integer.base.b_C <- function(x, k) {
  #x is the number we want to express in binary
  #k is the number positions we need
  if (x == 0)
    return(rep(0, k))
  else{
    ndigits <- (floor(logb(x, base = 2)) + 1)
    res <- rep(0, ndigits)
    for (i in 1:ndigits) {
      #i <- 1
      res[i] <- (x %% 2)
      x <- (x %/% 2)
    }
    return(c(res, rep(0, k - ndigits)))
  }
}

