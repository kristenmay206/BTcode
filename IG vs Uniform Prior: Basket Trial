# July 11, 2017
# Manuscript: Prior specification of the variance parameters for a basket trial design using Bayesian hierarchical modeling
# To run: you must create three files from the code below (a '.R' script and two '.txt' scripts)
# The two '.txt' scripts contain the JAGS code, which is embedded at bottom
# Save in separate files and name the files 'model1.txt' (for inverse-gamma prior) and 'model2.txt' (for uniform prior)

library(parallel)
library(rjags)
no.cores = 9 # no. of cpu cores used in mclapply
# note: as specified IG prior simulation study runs 4.5 hours
# note: as specified Uniform prior simulation study runs 6.5 hours

# output results (csv files)
#setwd("G:/My Documents/Basket Trials/Paper Code")
########################################################################################################
# define simulated trial function 
trial = function(it, h, K, gammaD, max.n, h0, h1, gammaF, gammaE, int1, intk, no.interim, AR, model,n.int.min,n.final.min, hyp1, hyp2){

# true response rates
if(h == 0){ p.sim = h0
	}else if(h == K){ p.sim = h1
	}else{ p.sim = c(h1[1:h],h0[(h+1):K]) }

p.mid = (h0+h1)/2 # mid pt
mu.h = mean(log(h1/(1-h1)) - log(h0/(1-h0))) # change in log odds from the null - lower level mean hyperparameter for mu
n.iter = 10000 # number of mcmc samples for inference
burn.in = 1000 # additional burn in + default 1000 in jags.model 
keep = matrix(NA, nrow = no.interim*3, ncol = K) # keep track of when baskets drop out
n.interim = matrix(NA, nrow = no.interim*3 + 1, ncol = K) # number/basket enrolled at each stage
N.basket = {}
yes.keep = {}

# define initial values for mcmc sampler for different models
if(model == 1){
	a = hyp2/2; b = hyp1**2*hyp2/2
	inits.list = list("prec" = 1, "mu" = mu.h)
}else if(model == 2){
	a = hyp1; b = hyp2
	inits.list = list("sigma" = 1, "mu" = mu.h)
}

# generate enrollment times for each basket
Baskets = matrix(NA,nrow = max.n*K*2,ncol = 2)
for (i in 1:K){
	Baskets[(i+(i-1)*(max.n*2-1)):(i*max.n*2),] = cbind(rexp(max.n*2,AR[i]),rep(i,max.n*2))
}
	ord = order(Baskets[,1])
	keep.bask = Baskets[ord,][1:(int1*K),]
	tab = as.numeric(table(factor(keep.bask[,2],levels = 1:K)))
	n.interim[1,] = tab
	N.basket = colSums(n.interim,na.rm = T)
	exclude.int = which(N.basket < n.int.min) # exclude these baskets from stopping rules at next interim
	
## generate data ##
yes = rbinom(K,n.interim[1,],p.sim); yes.keep = rbind(yes.keep,yes)

## fit model ##
fit = jags.model(file = paste("model",model,".txt",sep = ""), data = list("resp" = yes, "n" = N.basket, "K" = K, "p.0" = h0,"mu.h" = mu.h, "a" = a, "b" = b), inits = inits.list ,quiet = T) # default: 1000 burn in 
update(fit,burn.in) # additional burn in
p.post = coda.samples(fit,variable.names = "p",n.iter = n.iter)[[1]]

cont.bask.t = {}; dec = rep(NA,K); cont.bask.e = {}
## interim analyses ##
	int.count = 1
	for (i in 1:K){
		if ( !(i %in% exclude.int) ) { 
			if (mean(p.post[,i] > p.mid[i]) < gammaF){ cont.bask.t[i] = 0 ; dec[i] = 0 # 0 means stop
			}else{  cont.bask.t[i] = 1 }
		if (mean(p.post[,i] > p.mid[i]) > gammaE){ cont.bask.e[i] = 0 ; dec[i] = 1
			}else{  cont.bask.e[i] = 1 }
		keep[int.count,i] = ifelse(cont.bask.t[i] == 1 & cont.bask.e[i] == 1, 1, NA)
		}else{ keep[int.count,i] = 1 } 
	}

	while(length(na.omit(dec)) != K){ 
		dropped = which(is.na(keep[int.count,]))
		current.size = sum(N.basket)
		K.star = which(!(1:K %in% dropped))
		temp = Baskets[ord,][-(1:current.size),]
		if(length(dropped) > 0){ new.Baskets = Baskets[ord,][-(current.size + which(temp[,2] %in% dropped)),]
		}else{ new.Baskets = Baskets[ord,]}
		keep.bask = new.Baskets[(current.size+1):(min(current.size+intk*length(K.star),dim(new.Baskets)[1])),]
		tab = as.numeric(table(factor(keep.bask[,2],levels = 1:K)))
		n.interim[(int.count+1),] = tab
		N.basket = colSums(n.interim,na.rm = T)
		exclude.int = which(N.basket < n.int.min) # exclude these baskets from stopping rules at next interim

		for (i in which(!is.na(keep[int.count,])) ){yes[i] = yes[i] + rbinom(1,n.interim[(int.count+1),i],p.sim[i])}
		yes.keep = rbind(yes.keep,yes)
		fit = jags.model(file = paste("model",model,".txt",sep = ""), data = list("resp" = yes, "n" = N.basket, "K" = K, "p.0" = h0,"mu.h" = mu.h, "a" = a, "b" = b), inits = inits.list ,quiet = T) # default: 1000 burn in 
		update(fit,burn.in) # additional burn in
		p.post = coda.samples(fit,variable.names = "p",n.iter = n.iter)[[1]]

		int.count = int.count+1
		for (i in 1:K){
		if (is.na(dec[i]) == T){ 
			# if max sample size not reached and can evaluate
			if ( N.basket[i] < max.n & !(i %in% exclude.int) ){
			 	if (mean(p.post[,i] > p.mid[i]) < gammaF){ cont.bask.t[i] = 0 ; dec[i] = 0 # 0 means stop
					}else{  cont.bask.t[i] = 1 }
				if (mean(p.post[,i] > p.mid[i]) > gammaE){ cont.bask.e[i] = 0 ; dec[i] = 1
					}else{  cont.bask.e[i] = 1 }
				keep[int.count,i] = ifelse((cont.bask.t[i] == 1 & cont.bask.e[i] == 1), 1, NA) # NA if stop
			}else if ( N.basket[i] < max.n & (i %in% exclude.int) ){
				keep[int.count,i] = 1 # continue with not enough SS to evaluate
			}else if ( N.basket[i] > n.final.min & N.basket[i] >= max.n ){ # final analysis 
				if (int.count < no.interim*2){keep[int.count,i] = NA}
				dec[i] = ifelse(mean(p.post[,i] > h0[i]) > gammaD,1,0)
			}
		} } 

	} # end while

	samples = coda.samples(fit, variable.names = c("mu","prec","sigma","theta"),n.iter = n.iter)[[1]] 
	
	return(c(dec,yes,yes/N.basket,colMeans(p.post),N.basket,sum(N.basket),colMeans(samples),keep,n.interim))

} # end trial

########################################################################################################
# define simulation study function 
sim = function(h, n.sim, K, gammaD, max.n, h0, delta, gammaF, gammaE, int1, intk, no.interim, AR, model, write.csv.ind,n.int.min,n.final.min, hyp1, hyp2){
	h1 = h0 + delta
	set.seed(4206)
	results.list = lapply(1:n.sim, trial, h, K, gammaD, max.n, h0, h1, gammaF, gammaE, int1, intk, no.interim, AR, model,n.int.min,n.final.min, hyp1, hyp2) 
	results = matrix(data = unlist(results.list),ncol = n.sim, nrow = 6*K + 4 + (no.interim*6+1)*K)
	if (write.csv.ind == 1) {write.csv(results,paste("Raw_results_model",model,"_sc",sc,"_KdeltagammaD",K,delta*100,gammaD*100,".csv",sep = ""))}
	reject  = round(rowMeans(results[1:K,]),2)
	fwer = round(ifelse(h < K - 1, sum(as.numeric(colSums(results[((1+h):(K)),],na.rm = T) >= 1))/n.sim, ifelse( h == K-1, sum(results[((1+h):(K)),],na.rm = T)/n.sim, NA ) ),2) # only one basket under the null
	obs = round(rowMeans(results[(2*K+1):(3*K),]),2)
	est = round(rowMeans(results[(3*K+1):(4*K),]),2)
	est.sd = round(apply(results[(3*K+1):(4*K),],1,sd),2)
	n = round(rowMeans(results[(4*K+1):(5*K),]))
	EN = round(mean(results[(5*K+1),]))
	mc.stuff = round(rowMeans(results[(5*K+2):(6*K+4),]),2)
	mc.sd.stuff = round(apply(results[(5*K+2):(6*K+4),],1,sd),2)
	keep = matrix(round(rowSums(results[(6*K+5):(6*K+4+no.interim*3*K),],na.rm = T)/n.sim,3),nrow = no.interim*3,ncol = K)
	n.interim = matrix(round(rowSums(results[(6*K+5+no.interim*3*K):(6*K+4+no.interim*6*K+K),],na.rm = T)/n.sim,3),nrow = no.interim*3+1,ncol = K)
	return(c(reject, fwer, obs, est, est.sd, n , EN, mc.stuff, mc.sd.stuff, keep, n.interim))
}

########################################################################################################
# specify simulation study parameters 

### number of simulated trials ###
n.sim = 1000

### number of baskets ###
K = 5

### scenarios ###
scenarios = 0:K

### select one: prior model & hyperparameters ###
# model = 1: IG(a,b) prior with a = hyp2/2, b = hyp1^2*hyp2/2
# model = 2: Uniform(hyp1,hyp2) prior

model = 1 # 2

if(model == 1){
	hyp1 = c(0.1,0.5,1,2,10) # prior mean for sigma2
	hyp2 = c(0.01,0.1,0.5,1,2,5,10) # prior weight for prior mean
}else if(model == 2){
	hyp1 = c(0.01,0.05,0.3,0.5,0.71) # lower bound
	hyp2 = c(1,2,3,10,100,10000) # upper bound
}

### trial design parameters ###
max.n = 20 # maximum number of patients/basket
h0 = rep(0.15,K)  # null/inactive response rate (null hypothesis) - can vary by basket
delta = 0.3 # treatment effect (must be common across all baskets)
gammaF = 0.05 # posterior threshold to stop for futility
gammaE = 0.9 # posterior threshold to stop for superior efficacy
int1 = 10 # avg no. of patients/basket req for first interim analysis
intk = 5 # avg no. of patients/basket between interim analyses
no.interim = round((max.n - int1)/intk) # number of interim analyses 
n.int.min = 10  # min. no. patients needed for interim analyses
n.final.min = 10  # min. no. patients needed to evaluated efficacy
AR = rep(2,K) # accrual rate/baskets (months)
write.csv.ind = 0 # write.csv of all raw results for every combination

### final decision rule for inference ###
# gammaD: posterior probability threshold that the RR is > null to declare drug works

### note: originally started with grid search here ###
#gammaD.seq = c(0.95, 0.99)
#combos = expand.grid(scenarios,m,w,gammaD.seq) 
#suggest using scenarios = 0 (or which ever scenario you want to calibrate under) to save computational time

combos = expand.grid(scenarios,hyp1,hyp2)
combos$gammaD = 0.95
if (model == 1){
# IG prior
	combos$gammaD[combos[,2] == 0.1] = 0.9
	combos$gammaD[combos[,2] == 0.1 & (combos[,3] == 0.5 | combos[,3] == 1)] = 0.89
	combos$gammaD[combos[,2] == 0.1 & combos[,3] == 2] = 0.85
	combos$gammaD[combos[,2] == 0.1 & combos[,3] == 5] = 0.84
	combos$gammaD[combos[,2] == 0.1 & combos[,3] == 10] = 0.845
	combos$gammaD[combos[,2] == 0.5] = 0.93
	combos$gammaD[combos[,2] == 0.5 & combos[,3] == 0.01] = 0.915
	combos$gammaD[combos[,2] == 0.5 & combos[,3] == 0.1] = 0.925
	combos$gammaD[combos[,2] == 0.5 & combos[,3] == 1] = 0.935
	combos$gammaD[combos[,2] == 0.5 & combos[,3] == 2] = 0.935
	combos$gammaD[combos[,2] == 1 & (combos[,3] == 0.01 | combos[,3] == 0.1)] = 0.93
	combos$gammaD[combos[,2] == 1 & (combos[,3] == 2 | combos[,3] == 10)] = 0.955
	combos$gammaD[combos[,2] == 2 & combos[,3] == 0.01] = 0.93
	combos$gammaD[combos[,2] == 2 & combos[,3] == 0.1] = 0.945
	combos$gammaD[combos[,2] == 2 & combos[,3] == 0.5] = 0.955
	combos$gammaD[combos[,2] == 2 & combos[,3] == 2] = 0.97
	combos$gammaD[combos[,2] == 2 & (combos[,3] == 1 | combos[,3] == 5 | combos[,3] == 10)] = 0.965
	combos$gammaD[combos[,2] == 10] = 0.975
	combos$gammaD[combos[,2] == 10 & combos[,3] == 0.01] = 0.955
	combos$gammaD[combos[,2] == 10 & combos[,3] == 0.1] = 0.97
	combos$gammaD[combos[,2] == 10 & combos[,3] == 1] = 0.97
	combos$gammaD[combos[,2] == 10 & combos[,3] == 10] = 0.98
}else if (model == 2){
# uniform prior
	combos$gammaD[combos[,2] == 0.01 | combos[,2] == 0.05] = 0.94
	combos$gammaD[combos[,2] == 0.01 & combos[,3] == 1] = 0.925
	combos$gammaD[combos[,2] == 0.01 & combos[,3] == 2] = 0.93
	combos$gammaD[combos[,2] == 0.01 & combos[,3] == 3] = 0.935
	combos$gammaD[combos[,2] == 0.05 & combos[,3] == 1] = 0.925
	combos$gammaD[combos[,2] == 0.05 & combos[,3] == 2] = 0.935
	combos$gammaD[combos[,2] == 0.3 & combos[,3] == 1] = 0.93
	combos$gammaD[combos[,2] == 0.3 & (combos[,3] == 100 | combos[,3] == 10000)] = 0.955
	combos$gammaD[combos[,2] == 0.5 & combos[,3] == 1] = 0.935
	combos$gammaD[combos[,2] == 0.5 & (combos[,3] == 3 | combos[,3] == 10 | combos[,3] == 10000)] = 0.955
	combos$gammaD[combos[,2] == 0.5 & combos[,3] == 100] = 0.96
	combos$gammaD[combos[,2] == 0.71] = 0.96
	combos$gammaD[combos[,2] == 0.71 & (combos[,3] == 1 | combos[,3] == 2)] = 0.955
	combos$gammaD[combos[,2] == 0.71 & combos[,3] == 3] = 0.96
	combos$gammaD[combos[,2] == 0.71 & (combos[,3] == 100 | combos[,3] == 10000)] = 0.965
}

########################################################################################################
# define function to organize results
output = function(ind){

	mat = matrix(1,nrow = K, ncol = K); mat[lower.tri(mat)] = 0
	avg.err.s = diag(t(all.results[1:K,(ind:(ind+K-1))])%*%t(mat))/seq(K,1)
	avg.err = round(mean(avg.err.s),2)
	rg.err = paste("(",round(range(avg.err.s)[1],2),",",round(range(avg.err.s)[2],2), ")", sep = "")

	avg.fwer = round(mean(as.numeric(all.results[(K+1),(ind:(ind+K-1))])),2)
	rg.fwer = paste("(",round(range(all.results[(K+1),(ind:(ind+K-1))])[1],2),",", round(range(all.results[(K+1),(ind:(ind+K-1))])[2],2),")",sep = "")

	avg.pow.s = diag(t(all.results[1:K,((ind+1):(ind+K))])%*%mat)/seq(1,K)
	avg.pow = round(mean(avg.pow.s),2)
	rg.pow = paste("(",round(range(avg.pow.s)[1],2),",",round(range(avg.pow.s)[2],2), ")", sep = "")

	avg.en = round(mean(as.numeric(all.results[(5*K+2),(ind:(ind+K))])))
	rg.en = paste("(",round(range(all.results[(5*K+2),(ind:(K+ind))])[1]),",", round(range(all.results[(5*K+2),(ind:(K+ind))])[2]),")",sep = "")

	return(c(avg.fwer, rg.fwer, avg.err,rg.err,avg.pow,rg.pow,avg.en,rg.en))
}

########################################################################################################
# run simulation study

start = Sys.time()
# using parallel processing
all.results.list = mclapply(1:nrow(combos),function(j, K, n.sim, max.n, h0, gammaF, gammaE, int1, intk, no.interim, AR, write.csv.ind, n.int.min, n.final.min, delta, model){sim(h = combos[j,1], n.sim, K, gammaD = combos[j,4], max.n, h0, delta, gammaF, gammaE, int1, intk, no.interim, AR, model, write.csv.ind, n.int.min, n.final.min, hyp1 = combos[j,2], hyp2 = combos[j,3])}, K = K, n.sim = n.sim, delta = delta, max.n = max.n, h0 = h0, gammaF = gammaF, gammaE = gammaE, int1 = int1, intk = intk, no.interim = no.interim, AR = AR, write.csv.ind = write.csv.ind, n.int.min = n.int.min, n.final.min = n.final.min, model = model, mc.silent = TRUE, mc.cores = no.cores, mc.preschedule = FALSE)

# run on desktop
#all.results.list = lapply(1:nrow(combos),function(j, K, n.sim, max.n, h0, gammaF, gammaE, int1, intk, no.interim, AR, write.csv.ind, n.int.min, n.final.min, delta, model){sim(h = combos[j,1], n.sim, K, gammaD = combos[j,4], max.n, h0, delta, gammaF, gammaE, int1, intk, no.interim, AR, model, write.csv.ind, n.int.min, n.final.min, hyp1 = combos[j,2], hyp2 = combos[j,3])}, K = K, n.sim = n.sim, delta = delta, max.n = max.n, h0 = h0, gammaF = gammaF, gammaE = gammaE, int1 = int1, intk = intk, no.interim = no.interim, AR = AR, write.csv.ind = write.csv.ind, n.int.min = n.int.min, n.final.min = n.final.min, model = model) #, mc.silent = TRUE, mc.cores = no.cores, mc.preschedule = FALSE)

# output
all.results = matrix(data = unlist(all.results.list),ncol = nrow(combos), nrow = 8 + (no.interim*6+8)*K) # col are combos
	write.csv(all.results,paste("AllResults_K",K,"baskets_Model",model,".csv", sep = ""),row.names = F)
print(Sys.time() - start)

ind = seq(1,nrow(combos),by = K+1)
gammaD.seq = combos[which(combos[,1] == 0),4]
out.combos = cbind(expand.grid(hyp1,hyp2),gammaD.seq)
out = cbind(out.combos,t(sapply(ind,output)))
colnames(out) = c("Hyp1","Hyp2","GammaD","AvgFWER","RgFWER","AvgErr","RgErr","AvgPow","RgPow","AvgEN","RgEN")
	write.csv(data.frame(out),paste("TabularResults_K",K,"baskets_Model",model,".csv", sep = ""),row.names = F)

# FWER under A = 0: check calibration
print(all.results[(K+1),ind])

############################################## JAGS CODE #####################################################
# JAGS code for Model 1: save in separate file as 'model1.txt'

model{
	for (k in 1:K){ # basket
		resp[k] ~ dbinom(p[k],n[k]) # data: resp and n
		logit(p[k]) <- theta[k] + logit(p.0[k])
		theta[k] ~ dnorm(mu,prec) # jags uses precision 
	}
		mu ~ dnorm(mu.h,0.1) 
		prec ~ dgamma(a,b) 
		sigma <- 1/sqrt(prec)
}

############################################## JAGS CODE #####################################################
# JAGS code for Model 2: save in separate file as 'model2.txt'

model{
	for (k in 1:K){ # basket
		resp[k] ~ dbinom(p[k],n[k]) # data: resp and n
		logit(p[k]) <- theta[k] + logit(p.0[k])
		theta[k] ~ dnorm(mu,prec) # jags uses precision 
	}
		mu ~ dnorm(mu.h,0.1) 
		prec <- 1/sigma^2
		sigma ~ dunif(a,b)
}

########################################################################################################
