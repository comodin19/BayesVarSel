resamplingSBSB<- function(modelslBF, positions){
	#Given the matrix created within BvsF containing a binary expression of models and the log(BF) in the last columnbtained 
	#with GibbsBvs(...prior.models = "SBSB"...)
	#(meaning that prior over the model space is obtained with the hierarchical
	#Pr(Mg)=Pr(Mg|those factors and x's)*Pr(those factors and x's)
	#and Pr(Mg|those factors and x's)\proto 1/(number of models of that dimension with those factors and x's)
	#and Pr(those factors and x's)\propto 1/(number of models of that dimension))
	#the function makes a resampling with the new prior being
	#Pr(Mg|those factors and x's)\propto 1/(number of different models of that rank with those factors and x's)
	#returning a matrix of the same size with the resampled models
	levels<- rowSums(positions)
	levelsf<- levels[levels>1]
	#ff contains the matrix with the count of models of each rank for each combination of levels
	ff<- matrix.rank.levels(all.levelsf=levelsf)
	kplusp<- dim(positions)[1]
	#A matrix containing, for each sampled models, the number of active "levels" for each regressor:
	m.actlevels<- t(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}))
	#for each sampled model obtain the SBSB2 prior prob:
	allmodelspriorSBSB2<- t(apply(m.actlevels, MARGIN=1, FUN=priorSBSB2))
	#for each sampled model obtain the SBSB1 prior prob:
	allmodelspriorSBSB1<- t(apply(m.actlevels, MARGIN=1, FUN=priorSBSB1))
	#now the resampling:
	resamp<- sample(x=1:dim(modelslBF)[1], size=dim(modelslBF)[1], rep=T, prob=exp(allmodelspriorSBSB2-allmodelspriorSBSB1))
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
	levels<- rowSums(positions)
	levelsf<- levels[levels>1]
	kplusp<- dim(positions)[1]
	#A matrix containing, for each sampled models, the number of active "levels" for each regressor:
	m.actlevels<- t(apply(modelslBF[,-dim(modelslBF)[2]], MARGIN=1, FUN=function(x){positions%*%x}))
	#for each sampled model obtain the ConstConst2 prior prob:
	allmodelspriorSBSB2<- t(apply(m.actlevels, MARGIN=1, FUN=ConstConst2))
	#for each sampled model obtain the SBSB1 prior prob:
	allmodelspriorSBSB1<- t(apply(m.actlevels, MARGIN=1, FUN=ConstConst1))
	#now the resampling:
	resamp<- sample(x=1:dim(modelslBF)[1], size=dim(modelslBF)[1], rep=T, prob=exp(allmodelspriorSBSB2-allmodelspriorSBSB1))
	return(modelslBF[resamp,])	
}




priorSBSB1<- function(act.levels){
	#The original SBSB prior
	lprMgamma<- -sum(log(levels[act.levels>0])+lchoose(levels[act.levels>0], act.levels[act.levels>0]))-log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0))
	return(lprMgamma)
}
	
priorSBSB2<- function(act.levels){
	#The corrected SB-SB prior prob (inversely proportional to the number of models of that rank)
	#For copies of the same model, we only keep one representative (the full on that class)
	#(four different cases)
	
	#of the active levels take only those that correspond to factors:	
	act.levelsf<- act.levels[levels>1]	
	
	#If the model does not contain any factor:
	if (sum(act.levelsf)==0) return(-log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0)))
	
	#if the model is not saturated nor oversaturated
	if (sum(act.levelsf >= (levelsf-1)) == 0){
		whomactive<- paste("f",(1:length(levelsf))[act.levelsf!=0], collapse="", sep="")
		l2prMgamma<- -log(ff[whomactive, sum(act.levelsf)])-
	              #log(sum(levelsf)-2*length(levelsf)+1)-
								log(sum(levelsf[act.levelsf>0])-2*length(levelsf[act.levelsf]>0)+1)-
								log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0))	
		return(l2prMgamma)	
			  }
  #if the model contains at least one level saturated then return -Inf:
	if (sum(act.levelsf[act.levelsf>0] == (levelsf[act.levelsf>0]-1)) >= 1){
		return(-Inf)
	}
	else { #keep the rest
	whomactive<- paste("f",(1:length(levelsf))[act.levelsf!=0], collapse="", sep="")
	l2prMgamma<- -log(ff[whomactive, sum(act.levelsf[act.levelsf==levelsf]-1)+sum(act.levelsf[act.levelsf<levelsf])])-
	              log(sum(levelsf[act.levelsf>0])-2*length(levelsf[act.levelsf>0])+1)-log(kplusp+1)-lchoose(kplusp, sum(act.levels!=0))
	return(l2prMgamma)	
  }
}


priorConstConst1<- function(act.levels){
	#The original Const-Const prior (conditionally, inversely proportional to the number of models)
	
	#of the active levels take only those that correspond to factors:	
	act.levelsf<- act.levels[levels>1]	
	
	lprMgamma<- -sum(log(2^levelsf[act.levelsf>0]-1))-kplusp*log(2)
	return(lprMgamma)
}


priorConstConst2<- function(act.levels){
	#The corrected Const-Const prior prob (proportional to a constant for unique models)
	#For copies of the same model, we only keep one representative (the full on that class)
	#(four different cases)
	#of the active levels take only those that correspond to factors:	
	act.levelsf<- act.levels[levels>1]	
	
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
	#Given a vector of ACTIVE levels (l1,l2,...,lR), this function computes how many
	#models there are of all possible ranks R<=r<=sum(l_i)-R
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

