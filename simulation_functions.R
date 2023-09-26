###simulation functions


prsSim = function(n, t, bprs, bst, beps, zStruct, zStructP){
  ### Function to simulate polygenic risk scores and phenotypes according to the Y_j = X beta_j + Z delta_j + epsilon_j model
  
  #n: number of individuals
  #t: number of polygenic risk scores
  #bprs: proportion of phenotypic variance attributable to polygenic risk score (X beta_j)
  #bst: proportion of phenotypic variance attributable to structural variation (Z delta_j)
  #beps: proportion of phenotypic variance attributable to random error (epsilon_j)
  #zStruct: distribution of Z term representing unobserved covariates. Three options:
    #1 -> unif(1,2)
    #2 -> unif(-2,2)
    #3 -> point-normal: Z~0 with probability zStructP, Z~N(1, 0.05) with probability 1-zStructP
  #zStructP: probability for point-normal distribution above
  
  ###create dataframe of PRS
  ###PRS (X beta_j) are simulated according to beta distribution, with shape hyperparameters each distributed unif(1, 10)
  prs = apply(X = data.frame(n = rep(n,t), shape1 = runif(t, 1, 10), shape2 = runif(t, 1, 10)),
             MARGIN=1,FUN=function(x){rbeta(x[1], x[2], x[3])})
  
  ###label columns
  prs = as.data.frame(prs)
  colnames(prs) = paste0("prs",1:t)
  
  ###scale PRS to have variance one
  prs = apply(prs, 2, FUN = function(x){x / sd(x)})
  
  ###simulate environmental factors
  if(zStruct==1){
    Z = runif(n, 1, 2)
    ###null placeholder
    chance = rep(NA, n)
  }else if(zStruct==2){
    Z = runif(n, -2, 2)
    ###null placeholder
    chance = rep(NA, n)
  }else if(zStruct==3){
    chance = rbinom(n, 1, zStructP)
    Z = rnorm(n,1, 0.05)
    Z[chance==0]=0
  }
  
  ###delta_j is drawn from N(2, .1) distribution
  deltas = rnorm(t, 2, .1)
  
  ###generate structural error
  zDelta = Z%o%deltas
  
  ###scale structural error to have variance one
  zDelta = apply(zDelta, 2, FUN = function(x){x / sd(x)})
  
  ###label columns
  zDelta = as.data.frame(zDelta)
  colnames(zDelta) = paste0("struct",1:t)
  
  ###generate random error
  epsilon = apply(X = data.frame(n = rep(n, t)), MARGIN=1, FUN=function(x){rnorm(x[1],0,1)})
  ###ensure variance is one for random error
  epsilon = apply(epsilon, 2, FUN = function(x){x / sd(x)})
  
  ###label columns
  epsilon = as.data.frame(epsilon)
  colnames(epsilon) = paste0("error",1:t)
  
  ###generate phenotypes
  phenos = matrix(NA,nrow = n, ncol = t)
  for(i in 1:t){
    phenos[,i] = sqrt(bprs)*prs[,i] + sqrt(bst)*zDelta[,i] + sqrt(beps)*epsilon[,i]
  }
  
  simData <- list("prs" = prs,
                  "struct" = zDelta,
                  "error" = epsilon,
                  "phenos" = phenos,
                  "clust" = chance,
                  "Z" = Z)
  
  return(simData)
}




getResiduals = function(phenos, prs){
  ###function to get residuals of univariate phenotype-on-prs regressions
  
  #phenos: n x t matrix of phenotypes
  #prs: n x t matrix of polygenic risk scores
  
  resids = matrix(NA, nrow=nrow(phenos),ncol = ncol(phenos))
  for(i in 1:ncol(phenos)){
    resids[,i] = residuals(lm(phenos[,i]~prs[,i]))
  }
  
  return(resids)
}




clustering = function(resids, clusters, method){
  ###function to perform clustering on the residuals and assess accuracy against the simulated bimodal structural error
  
  #resids: n x t matrix of the residuals of univariate phenotype-on-prs regressions
  #clusters: n x 1 vector of true cluster assignment (coded 0 or 1)
  #method: clustering method to apply. Currently, only k-means is implemented
  
  if(method == "kmeans"){
    estClust = kmeans(resids, centers = 2, nstart = 3)
    accuracy = adjustedRandIndex(estClust$cluster, clusters)
  }
  
  return(accuracy)
}

pcaModeling = function(resids, Z){
  ###function to perform multiple linear regression of the environmental error against the first five principal components
  
  #resids: n x t matrix of the residuals of univariate phenotype-on-prs regressions
  #Z: n x 1 vector of ground truth structural covariate 
  
  pcs = prcomp(resids)
  
  rsqs = apply(X = data.frame(pc = 1:5), MARGIN=1,
        FUN = function(y){summary(lm(Z~pcs$x[,y]))$r.squared})
  
  return(rsqs)
}