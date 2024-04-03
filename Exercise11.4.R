data <- as.matrix(read.table("mathstandard.dat"))
data <- data[-1,]

data <- as.data.frame(data)
data$V2 <- as.numeric(data$V2)
data$V3 <- as.numeric(data$V3)

ids<-sort(unique(data$V1)) 
m<-length(ids)
Y<-list() ; X<-list() ; N<-NULL
for(j in 1:m) 
{
  Y[[j]]<-data[data$V1==ids[j], 2] 
  N[j]<- sum(data$V1==ids[j])
  xj<-data[data$V1==ids[j], 3] 
  xj<-(xj-mean(xj))
  X[[j]]<-cbind( rep(1,N[j]), xj  )
}

#### OLS fits, same as ML for regression
S2.LS<-BETA.LS<-NULL
for(j in 1:m) {
  fit<-glm(Y[[j]]~-1+X[[j]],family = "binomial" )
  BETA.LS<-rbind(BETA.LS,c(fit$coef)) 
  S2.LS<-c(S2.LS, summary(fit)$sigma^2) 
} 

BETA.LS[7,2] = 0
BETA.LS[12,2] = 0
BETA.LS[30,2]=0
BETA.LS[35,2] = 0

#Logistic regression calculation

probabilities <- lapply(1:m, function(j) {
  exp_vals <- exp(BETA.LS[j, 1] + BETA.LS[j, 2] * X[[j]][, 2])
  probabilities <- exp_vals / (1 + exp_vals)
  return(probabilities)
})

pdf("LogRegPlots",family="Times",height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))

plot( range(data[,3]),range(data[,2]),type="n",xlab="s", 
      ylab="math score")
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,2],col="gray")  }
dev.off()


#Need to find the counties with 10 or more entries fo ML estimation

which(sapply(X, function(x) nrow(x) >= 10))

#Getting the beta values for these

TEN_BETA0 <- matrix(c(BETA.LS[c(6, 13, 14, 17, 18, 21, 27, 31, 32, 34, 37, 39), 1]), ncol = 1, nrow = 12)
TEN_BETA1 <- matrix(c(BETA.LS[c(6, 13, 14, 17, 18, 21, 27, 31, 32, 34, 37, 39), 2]), ncol = 1, nrow = 12)

MLS_BETA <- matrix(c(TEN_BETA0,TEN_BETA1),ncol=2,nrow=12)

BETA.MLS<-apply(MLS_BETA,2,mean)
BETA.MLS

#ad hoc estimates for sigma and theta

SIGMA0 <- cov(MLS_BETA)
SIGMA0

thetaMLS <- BETA.MLS
thetaMLS

#Functions for unit informative priors

## Wishart simulation
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}

#### Hierarchical regression model

## mvnormal simulation
rmvnorm<-function(n,mu,Sigma)
{ 
  E<-matrix(rnorm(n*length(mu)),n,length(mu))
  t(  t(E%*%chol(Sigma)) +c(mu))
}

mu0 <- apply(BETA, 2, mean)
S0 <- cov(BETA)
eta0 <- p + 2
iL0 <- iSigma <- solve(S0)
BETA <- BETA.LS

## MCMC
THETA.post <- NULL
set.seed(1)
for (s in 1:50000) {
  ## Update theta
  Lm <- solve(iL0 + m * iSigma)
  mum <- Lm %*% (iL0 %*% mu0 + iSigma %*% apply(BETA, 2, sum))
  theta <- t(rmvnorm(1, mum, Lm))
  ##
  ## Update Sigma
  mtheta <- matrix(theta, m, p, byrow = TRUE)
  iSigma <- rwish(1, eta0 + m,
                  solve(S0 + t(BETA - mtheta) %*% (BETA - mtheta)))
  ##
  ## Update beta
  Sigma <- solve(iSigma)
  dSigma <- det(Sigma)
  for (j in 1:m) {
    beta.p <- t(rmvnorm(1, BETA[j, ], 0.5 * Sigma))
    lr <- sum(log(dbinom(Y[j,, size = 1, prob = plogis(X[, , j] %*% beta.p))) -
               log(dbinom(Y[j, ], size = 1, prob = plogis(X[, , j] %*% BETA[j, ])))) +
          ldmvnorm(t(beta.p), theta, Sigma, iSigma = iSigma, dSigma = dSigma) -
          ldmvnorm(t(BETA[j, ]), theta, Sigma, iSigma = iSigma, dSigma = dSigma)
    if (log(runif(1)) < lr) { BETA[j, ] <- beta.p }
  }
  ##
  ## Store some output
  if (s %% 10 == 0) { THETA.post <- rbind(THETA.post, t(theta)) }
}