library('Matrix')
library(MASS)
library(data.table)
##generate related individuals
og_mat <- matrix(c(1, 0, 0, 0.5, 0.5, 0, 0, 1, 0, 0.25, 0.5, 0.25, 0, 0.5, 0.5, 0, 0, 1), ncol=3, byrow=T)

mt_to_sib <- t(apply(og_mat, 1, function(x) c(x[1]^2, x[1]*x[2], x[1]*x[2], x[1]*x[3], x[1]*x[3], x[2]^2, x[2]*x[3], x[2]*x[3], x[3]^2)))

mtType <- function(x){
	s <- sum(x)
	if(s==4){
		return(1)
	}
	if(s==3){
		return(2)
	}
	if(s==2){
		pro <- x[1]*x[2]
		if(pro==0){
			return(3)
		}
		if(pro==1){
			return(4)
		}
	}
	if(s==1){
		return(5)
	}
	if(s==0){
		return(6)
	}
}

##simulate sibpair genotypes

##check the probability matrix and corresponding indices
simSib <- function(nRep, f, p, p2){
	p1 <- 2*(p-p2)
	p0 <- 1-p1-p2
	parent <- rmultinom(nRep, size=f, prob=c(p2^2, 2*p1*p2, 2*p0*p2, p1^2, 2*p0*p1, p0^2))
#	print(parent)
	sibgeno <- apply(parent, 2, a_fn)
	#print(sibgeno)
	ramsib <- apply(sibgeno, 2, function(x) rep(1:9, x)[sample(1:f, f)])
	sib_mat <- apply(ramsib, 2, function(x) as.numeric(sapply(x, b_fn)))
	return(sib_mat)                                            
	}

a_fn <- function(x){
	mat1 <- rmultinom(1, size=x[1], prob=mt_to_sib[1,])
	mat2 <- rmultinom(1, size=x[2], prob=mt_to_sib[2,])
	mat3 <- rmultinom(1, size=x[3], prob=mt_to_sib[3,])
	mat4 <- rmultinom(1, size=x[4], prob=mt_to_sib[4,])
	mat5 <- rmultinom(1, size=x[5], prob=mt_to_sib[5,])
	mat6 <- rmultinom(1, size=x[6], prob=mt_to_sib[6,])
	return(mat1+mat2+mat3+mat4+mat5+mat6)
}

b_fn <- function(x){
	(x==1)*(c(2,2))+(x==2)*c(2,1)+(x==3)*c(1,2)+(x==4)*c(2,0)+(x==5)*c(0,2)+(x==6)*c(1,1)+(x==7)*c(0,1)+(x==8)*c(1,0)+(x==9)*c(0,0)
}


##simulate parent-offspring genotypes
simPO <- function(f, p, p2, HWE=T){
	if(HWE==T){
		p2 <- p^2
	}
	p1 <- 2*(p-p2)
	p0 <- 1-p1-p2                                               
	parent <- sample(0:2, size=(2*f), replace=T, prob=c(p0, p1, p2))
	parent_mat <- matrix(parent, ncol=2, byrow=T)
	pomat <- matrix(nrow=f, ncol=2)
	pomat[,1] <- parent_mat[,1]
	for(i in 1:f){
		pomat[i,2] <- sample(2:0, size=1, prob=og_mat[mtType(parent_mat[i,]),])
	}
	po_G <- as.numeric(t(pomat))
	return(po_G)
	#return(pomat)
}

##compute covariance matrix of Parent-offspring
covPO <- function(g){
	v1 <- var(g[,1])
	v2 <- var(g[,2])
	cov1 <- mean(g[,1]*g[,2])-mean(g[,1])*mean(g[,2])
	matrix(c(v1, cov1, cov1, v2)/v2, ncol=2, byrow=T)
}



##simulate singletons
simRamSingle <- function(nS, p, p2, HWE=T){
	if(HWE==T){
		p2 <- p^2
	}
	p1 <- 2*(p-p2)
	p0 <- 1-p1-p2
	sample(0:2, size=nS, replace=T, prob=c(p0, p1, p2))
}

simSingle <- function(nRep, nS, p, p2, HWE=T){
	if(HWE==T){
		p2 <- p^2
	}
	p1 <- 2*(p-p2)
	p0 <- 1-p1-p2
	count_mat <- rmultinom(n=nRep, size=nS, prob=c(p0, p1, p2))
	return(apply(count_mat, 2, function(x) rep(c(0,1,2), c(x[1], x[2], x[3]))))
}


##simulate phenotypes
simPheno <- function(gCausal, alpha, beta, h2, p, type='causal', corZ){
	sigma <- sqrt(2*p*(1-p)*beta^2*(1-h2)/h2)
	if(type=='causal'){
		y <- alpha + beta*gCausal + replicate(ncol(gCausal), rnorm(nrow(gCausal), mean=0, sd=sigma))
		return(y)
	}
	if(type=='theory'){
		y <- replicate(ncol(gCausal), rnorm(nrow(gCausal), mean=0, sd=sigma))
		return(y)
	}
	if(type=='GE'){
		f <- nrow(gCausal)/2
		Z <- replicate(ncol(gCausal), as.numeric(t(mvrnorm(f, c(0, 0), sigma*matrix(c(1, corZ, corZ, 1), ncol=2)))))
		y <- alpha + beta*gCausal + Z
		return(y)
	}
	if(type=='env'){
		f <- nrow(gCausal)/2
		Z <- replicate(ncol(gCausal), as.numeric(t(mvrnorm(f, c(0, 0), sigma*matrix(c(1, corZ, corZ, 1), ncol=2)))))
		y <- alpha + Z
		return(y)
	}	
}


##convert mating types to sibpair genotypes


##T1E of LMM for sibpairs
lmm_anova <- function(y, g, h2=1, cType='h2'){
	y_mat <- matrix(y, ncol=1); g2_mat <- matrix(g, ncol=1)
	f <- length(y)/2
	if(cType=='est'){
		ymat <- matrix(y, ncol=2, byrow=T)
		h2 <- mean((ymat[,1]-mean(y))*(ymat[,2]-mean(y)))/mean((y-mean(y))^2)
	}
	Sys.time()
	small_mat <- chol(solve(matrix(c(1, 0.5*h2, 0.5*h2, 1), ncol=2)))
	Sys.time()
#	print(small_mat)
	covmat <- as.matrix(bdiag(replicate(f, list(small_mat))))
	Sys.time()
	cholY <- covmat%*%y_mat
	cholG <- covmat%*%g2_mat
	cholone <- covmat%*%matrix(rep(1, 2*f), ncol=1)
	Sys.time()
	fit1 <- lm(cholY~cholG+cholone-1)
	fit0 <- lm(cholY~cholone-1)
	t1e <- anova(fit1, fit0, test='Rao')$'Pr(>Chi)'[2]
	return(t1e)
}


rev_anova <- function(y, g, h2=1){
	y_mat <- matrix(y, ncol=1); g2_mat <- matrix(g, ncol=1)
	f <- length(y)/2
	xmat <- matrix(g, ncol=2, byrow=T)
	num <- sum((xmat[,1]-mean(g))*(xmat[,2]-mean(g)))
	den <- sum((g-mean(g))^2)
	corr <- num/den*2
	small_mat <- chol(solve(matrix(c(1, corr, corr, 1), ncol=2)))
#	print(small_mat)
	covmat <- as.matrix(bdiag(replicate(f, list(small_mat))))
	cholY <- covmat%*%y_mat
	cholG <- covmat%*%g2_mat
	cholone <- covmat%*%matrix(rep(1, 2*f), ncol=1)
	fit1 <- lm(cholY~cholG+cholone-1)
	fit0 <- lm(cholY~cholone-1)
	t1e <- anova(fit1, fit0, test='Rao')$'Pr(>Chi)'[2]
	return(t1e)
}


lmm_la <- function(x, y, corr, cType='h2'){
	f <- length(x)/2
	if(cType=='est'){
		ymat <- matrix(y, ncol=2, byrow=T)
		corr <- mean((ymat[,1]-mean(y))*(ymat[,2]-mean(y)))/mean((y-mean(y))^2)
	}
	small <- solve(matrix(c(1, corr, corr, 1), ncol=2))
	sigma <- as.matrix(bdiag(replicate(f, list(small))))
	xM <- matrix(x, ncol=1); yM <- matrix(y, ncol=1)
	oneV <- matrix(rep(1, 2*f), ncol=1)
	xBar <- solve(t(oneV)%*%sigma%*%oneV)%*%t(oneV)%*%sigma%*%xM
	yBar <- solve(t(oneV)%*%sigma%*%oneV)%*%t(oneV)%*%sigma%*%yM
	varX <- t(xM-oneV%*%xBar)%*%sigma%*%(xM-oneV%*%xBar)
	varY <- t(yM-oneV%*%yBar)%*%sigma%*%(yM-oneV%*%yBar)
	covXY <- t(xM)%*%sigma%*%(yM - oneV%*%yBar)
	t_stat <- (covXY[1,1])^2/varX[1,1]/varY[1,1]*2*f
	return(pchisq(t_stat, df=1, lower.tail=F))
}


rev_la <- function(x, y){
	f <- length(x)/2
	xmat <- matrix(x, ncol=2, byrow=T)
	num <- sum((xmat[,1]-mean(x))*(xmat[,2]-mean(x)))
	den <- sum((x-mean(x))^2)
	corr <- num/den*2
	small <- solve(matrix(c(1, corr, corr, 1), ncol=2))
	sigma <- as.matrix(bdiag(replicate(f, list(small))))
	xM <- matrix(x, ncol=1); yM <- matrix(y, ncol=1)
	oneV <- matrix(rep(1, 2*f), ncol=1)
	xBar <- solve(t(oneV)%*%sigma%*%oneV)%*%t(oneV)%*%sigma%*%xM
	yBar <- solve(t(oneV)%*%sigma%*%oneV)%*%t(oneV)%*%sigma%*%yM
	varX <- t(xM-oneV%*%xBar)%*%sigma%*%(xM-oneV%*%xBar)
	varY <- t(yM-oneV%*%yBar)%*%sigma%*%(yM-oneV%*%yBar)
	covXY <- t(xM)%*%sigma%*%(yM - oneV%*%yBar)
	t_stat <- (covXY[1,1])^2/varX[1,1]/varY[1,1]*2*f
	return(pchisq(t_stat, df=1, lower.tail=F))
}

lmmT1E_la <- function(nSim, f, alpha, beta, pCausal, pTested, h2, yType, gType){
	gCausal <- simSib(1, f, p=pCausal[1], p2=pCausal[2])
	gTested <- simSib(nSim, f, p=pTested[1], p2=pTested[2])
	y <- simPheno(gCausal, alpha, beta, h2, pCausal[1], type=yType)
	pval_lmm <- apply(gTested, 2, lmm_la, y=y, corr=0.5*h2)
	pval_h2est <- apply(gTested, 2, lmm_la, y=y, cType='est')
	pval_rev <- apply(gTested, 2, rev_la, y=y)
	return(data.frame(pval_lmm, pval_h2est, pval_rev))
}

lmmT1E_anova <- function(nSim, f, alpha, beta, pCausal, pTested, h2, yType, gType){
	gCausal <- simSib(1, f, p=pCausal[1], p2=pCausal[2])
	gTested <- simSib(nSim, f, p=pTested[1], p2=pTested[2])
	y <- simPheno(gCausal, alpha, beta, h2, pCausal[1], type=yType)
	pval_lmm <- apply(gTested, 2, lmm_anova, y=y, h2)
	pval_h2est <- apply(gTested, 2, lmm_anova, y=y, h2, cType='est')
	pval_rev <- apply(gTested, 2, rev_anova, y=y, h2)
	return(data.frame(pval_lmm, pval_h2est, pval_rev))
}

var_fast <- function(x, corr){
	xbar <- mean(x)
	xmat <- matrix(x, ncol=2, byrow=T)
	sum((x-xbar)^2) - 2*corr*sum((xmat[,1]-xbar)*(xmat[,2]-xbar))
}

lmm_fast <- function(y, g, h2=1, cType='h2'){
	ybar <- mean(y); gmat <- matrix(g, ncol=2, byrow=T)
	if(cType=='h2'){
		corr <- h2*0.5
	}
	if(cType=='est'){
		ymat <- matrix(y, ncol=2, byrow=T)
		corr <- mean((ymat[,1]-ybar)*(ymat[,2]-ybar))/mean((y-ybar)^2)
	}
	if(cType=='rev'){
		gbar <- mean(g)
		corr <- mean((gmat[,1]-gbar)*(gmat[,2]-gbar))/mean((g-gbar)^2)
	}
	gcovmat <- apply(gmat, 1, function(x) c(x[1]-corr*x[2], x[2]-corr*x[1]))
	top <- sum(as.numeric(gcovmat)*(y-ybar))
	tStat <- top^2/var_fast(g, corr)/var_fast(y, corr)*length(g)
	return(pchisq(tStat, df=1, lower.tail=F))
}


lmmT1E_fast <- function(nSim, f, alpha, beta, pCausal, pTested, h2, yType){
	#print(Sys.time()); print('start')
	gCausal <- simSib(nSim, f, p=pCausal[1], p2=pCausal[2])
	gTested <- simSib(nSim, f, p=pTested[1], p2=pTested[2])
	y <- simPheno(gCausal, alpha, beta, h2, pCausal[1], type=yType)
	#print(Sys.time()); print('assoc')
	pval <- matrix(nrow=nSim, ncol=3)
	for(i in 1:nSim){
		pval[i,1] <- lmm_fast(y[,i], gTested[,i], h2=h2, cType='h2')
		pval[i,2] <- lmm_fast(y[,i], gTested[,i], h2=h2, cType='est')
		pval[i,3] <- lmm_fast(y[,i], gTested[,i], h2=h2, cType='rev')
	}	
	#print(Sys.time()); print('done')
	return(pval)
}

perm_fast <- function(yOld, g, permID, h2=1, cType='h2'){
	yOldMat <- matrix(yOld, ncol=2, byrow=T)
	ymat <- yOldMat[permID,]; y <- as.numeric(t(ymat))
	ybar <- mean(y); gmat <- matrix(g, ncol=2, byrow=T)
	if(cType=='h2'){
		corr <- h2*0.5
	}
	if(cType=='est'){
		corr <- mean((ymat[,1]-ybar)*(ymat[,2]-ybar))/mean((y-ybar)^2)
	}
	if(cType=='rev'){
		gbar <- mean(g)
		corr <- mean((gmat[,1]-gbar)*(gmat[,2]-gbar))/mean((g-gbar)^2)
	}
	gcovmat <- apply(gmat, 1, function(x) c(x[1]-corr*x[2], x[2]-corr*x[1]))
	top <- sum(as.numeric(gcovmat)*(y-ybar))
	tStat <- top^2/var_fast(g, corr)/var_fast(y, corr)*length(g)
	return(pchisq(tStat, df=1, lower.tail=F))
}

lmmT1E_perm <- function(nSim, nPerm, f, alpha, beta, pCausal, pTested, h2, yType){
	print(Sys.time()); print('start')
	gCausal <- simSib(nSim, f, p=pCausal[1], p2=pCausal[2])
	y <- simPheno(gCausal, alpha, beta, h2, pCausal[1], type=yType)
	pList <- list()
	for(i in 1:nSim){
		print(Sys.time()); print('start loop')
		Onepval <- matrix(nrow=nPerm, ncol=3)
		Cpval <- matrix(nrow=nPerm, ncol=3)
		gTested <- simSib(nPerm, f, p=pTested[1], p2=pTested[2])
		yPermID <- replicate(nPerm, sample(1:f, size=f))#; print(dim(yPermID))
		for(j in 1:nPerm){
			Onepval[j,1] <- lmm_fast(y[,i], gTested[,j], h2=h2, cType='h2')
			Onepval[j,2] <- lmm_fast(y[,i], gTested[,j], h2=h2, cType='est')
			Onepval[j,3] <- lmm_fast(y[,i], gTested[,j], h2=h2, cType='rev')
			Cpval[j,1] <- perm_fast(y[,i], gTested[,j], permID=yPermID[,j], h2=h2, cType='h2')
			Cpval[j,2] <- perm_fast(y[,i], gTested[,j], permID=yPermID[,j], h2=h2, cType='est')
			Cpval[j,3] <- perm_fast(y[,i], gTested[,j], permID=yPermID[,j], h2=h2, cType='rev')
		}
		pList[[i]] <- list(Onepval, Cpval)
	}	
	print(Sys.time()); print('done')
	return(pList)
}

##simulation and application on CF
sib_cor <- function(geno){
	if(sum(is.na(geno)) >0){
		return(c(mean(geno, na.rm=T)/2, mean(geno==2, na.rm=T), NA, NA))
	}
	if(var(geno) <= 0){
		return(c(mean(geno, na.rm=T)/2, mean(geno==2, na.rm=T), NA, NA))
	}
	g <- matrix(geno, ncol=2, byrow=T)
	alpha <- mean(g); p <- alpha/2; p2 <- mean(geno==2)
	num <- mean((g[,1]-alpha)*(g[,2]-alpha))#; print(num)
	den <- mean((geno-alpha)^2)#; print(den)
	A <- num/den; varA <- (9*p*(1-p)+3)/(32*p*(1-p))
	rho <- A/0.5 - 1
	n <- length(geno); f <- n/2
	HWE_test <- (A-0.5)/sqrt(varA/f)
	return(c(p, p2, rho, HWE_test, -p/(1-p)))
}

##MAF estimation
maf_mcpeek <- function(g, f, n, phi=0.25, mode='onetype'){
	smallcov <- solve(matrix(c(1, 2*phi, 2*phi, 1), ncol=2, byrow=T))
	one_mat <- matrix(rep(1, length(g)), ncol=1)
	if(mode == 'onetype'){
		bigcov <- as.matrix(bdiag(replicate(f, list(smallcov))))
	}
	if(mode == 'mix'){
		smallCov <- as.matrix(bdiag(replicate(f, list(smallcov))))
		bigcov <- as.matrix(bdiag(list(smallCov, diag(n))))
	}
	maf <- solve(t(one_mat)%*%bigcov%*%one_mat)%*%t(one_mat)%*%bigcov%*%g
	return(maf)
}

maf_reverse <- function(g, mode='sibpair'){
##g: numeric, f pairs of individuals. If parent-child, parent genotype comes first
	if(mode == 'sibpair'){
		maf <- mean(g)
	}
	if(mode == 'parent-child'){
		PCmat <- matrix(g, ncol=2, byrow=T)
		gP <- PCmat[,1]; gC <- PCmat[,2]
		maf <- (1-0.5*a)*mean(gP) + 0.5*a*mean(gC)
	}
	if(mode == 'mix-PC'){
		maf <- NA
	}
	if(mode == 'mix-Sib'){
		maf <- NA
	}
	return(maf)
}

##estimate parent-offspring maf
find_D <- function(gmat, a, g1bar, g2bar){
	alpha <- (1-0.5*a)*g1bar + 0.5*a*g2bar
	A <- sum((gmat[,1]-alpha)^2)
	B <- sum((gmat[,2]-alpha)^2)
	C <- sum((gmat[,1]-alpha)*(gmat[,2]-alpha))
	D <- A/(0.5*A + B - C)
	delta <- abs(D-a)
	return(delta)
}

solve_mle2 <- function(g){
	geno <- matrix(g, ncol=2, byrow=T)
	g1bar <- mean(geno[,1]); g2bar <- mean(geno[,2])
	s <- 0; e <- 2
	fs <- find_D(geno, s, g1bar, g2bar); fe <- find_D(geno, e, g1bar, g2bar)	
	while(min(fs, fe) > 0.00001){
		if(fs > fe){
			s <- (s+e)/2
			fs <- find_D(geno, s, g1bar, g2bar)
		}else{
			e <- (s+e)/2
			fe <- find_D(geno, e, g1bar, g2bar)
		}
	}
	return(ifelse(fs<fe, s, e))
}

