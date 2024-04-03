diabetes <- as.matrix(read.table("https://www2.stat.duke.edu/~pdh10/FCBS/Exercises/azdiabetes.dat"))
diabetes

# Extracting number of pregnancies, blood pressure, body mass index, diabetes pedigree, and age variables
X <- diabetes[, c(1, 3, 5, 6, 7)]

# Data conversion from character to number
X <- X[-1, ]
X <- apply(X, 2, function(x) as.numeric(as.character(x)))

# Removing missing values
X <- na.omit(X)

# Creating the response variable y
y <- (diabetes[, 8] == 'Yes') * 1

#Need to remove top entry as its a placeholder for the name
y <- y[-1]
y

#Standardizing the data

for (i in 1:dim(X)[2]){
  X[,i] <- (X[,i]-mean(X[,i]))/sd(X[,i])
}

#Sets a threshold for log values so that NA cannot be produced

fix = function(x) {
  ind = is.nan(x)|is.na(x)
  x[ind] = 0
  x
}

#Regression functions

logit = function(p) 
  {log(p/(1-p))}

expit = function(x) 
{exp(x)/(1+exp(x))}


update_beta0 = function(gamma, beta, beta0) {
  beta0_p =rnorm(1, beta0,sd =2)# propose beta0
  pi =expit(beta0 +as.matrix(X[,gamma==1]) %*% beta[gamma==1])
  pi_p =expit(beta0_p +as.matrix(X[,gamma==1]) %*% beta[gamma==1])
  log_a =sum(fix(y*log(pi_p))+fix((1-y)*log(1-pi_p))) +dnorm(beta0_p,mean=0,sd=4,log=TRUE)
  log_b =sum(fix(y*log(pi))+fix((1-y)*log(1-pi))) +dnorm(beta0,mean=0,sd=4,log=TRUE)
  log_r = log_a - log_b
  u =runif(1)
  beta0_new =ifelse(log(u)<log_r, beta0_p, beta0)
  return(beta0_new)
}

update_beta = function(gamma, beta, beta0) {
  p = length(beta)
  for(j in 1:p) {
    beta_p = beta
    beta_p[j] = rnorm(1, beta[j], sd = 1) # proposal beta
    pi = expit(beta0 + as.matrix(X[,gamma==1]) %*% beta[gamma==1]) # pi_original
    pi_p = expit(beta0 + as.matrix(X[,gamma==1]) %*% beta_p[gamma==1]) # pi_proposed
    log_a = sum(fix(y*log(pi_p))+fix((1-y)*log(1-pi_p))) + dnorm(beta_p[j], mean=0, sd=2, log=TRUE) # log numerator
    log_b = sum(fix(y*log(pi))+fix((1-y)*log(1-pi))) + dnorm(beta[j], mean=0, sd=2, log=TRUE) # log denominator
    log_r = log_a - log_b
    u = runif(1)
    beta[j] = ifelse(log(u)<log_r, beta_p[j], beta[j])
  }
  return(beta)
}

update_gamma = function(gamma, beta, beta0) {
  # randomly choose update order
  p = length(gamma)
  for(j in sample(p)) {
    gamma_a = gamma_b = gamma
    gamma_a[j] = 1 # for numerator
    gamma_b[j] = 0 # for denominator
    pi_a = expit(beta0 + as.matrix(X[,gamma_a==1]) %*% beta[gamma_a==1]) # pi_original
    log_a = sum(fix(y*log(pi_a))+fix((1-y)*log(1-pi_a)))# log numerator
    pi_b = expit(beta0 + as.matrix(X[,gamma_b==1]) %*% beta[gamma_b==1]) # pi_original
    log_b = sum(fix(y*log(pi_b))+fix((1-y)*log(1-pi_b)))# log numerator
    log_odds = log_a - log_b
    u = runif(1)
    gamma[j] = ifelse(u < expit(log_odds), 1, 0)
  }
  return(gamma)
}
# initial values
p = 5
gamma = rep(1, p)
beta = rep(0, p)
beta0 = 1
S = 10500
B = 500 # Burn-in
Gamma = matrix(NA, nrow = S, ncol = p)
Beta = matrix(NA, nrow = S, ncol = p)
Beta0 = rep(NA, S)
# update parameters
for(i in 1:S) {
  beta0 = update_beta0(gamma, beta, beta0)
  beta = update_beta(gamma, beta, beta0)
  gamma = update_gamma(gamma, beta, beta0)
  Beta0[i] = beta0
  Beta[i,] = beta
  Gamma[i,] = gamma
}

Beta0 = Beta0[-(1:B)]
Beta = Beta[-(1:B),]
Gamma = Gamma[-(1:B),]

# A plot for evaluating lack of convergence
stationarity.plot<-function(x,...){
  
  S<-length(x)
  scan<-1:S
  ng<-min( round(S/100),10)
  group<-S*ceiling(ng*scan/S)/ng
  
  boxplot(x~group,...)
}

#Autocorrelation plots to test convergence of mcmc

pdf("Stationarity105Beta.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,5),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))

stationarity.plot(Beta[,1],xlab="iteration",ylab=expression(beta[1]))
stationarity.plot(Beta[,2],xlab="iteration",ylab=expression(beta[2]))
stationarity.plot(Beta[,3],xlab="iteration",ylab=expression(beta[3]))
stationarity.plot(Beta[,4],xlab="iteration",ylab=expression(beta[4]))
stationarity.plot(Beta[,5],xlab="iteration",ylab=expression(beta[5]))
dev.off()



pdf("Stationarity105BetaGamma.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,5),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))

stationarity.plot(Beta[,1]*Gamma[,1],xlab="iteration",ylab=expression(paste(beta[1]*gamma[1])))
stationarity.plot(Beta[,2]*Gamma[,2],xlab="iteration",ylab=expression(beta[2]*gamma[2]))
stationarity.plot(Beta[,3]*Gamma[,3],xlab="iteration",ylab=expression(beta[3]*gamma[3]))
stationarity.plot(Beta[,4]*Gamma[,4],xlab="iteration",ylab=expression(beta[4]*gamma[4]))
stationarity.plot(Beta[,5]*Gamma[,5],xlab="iteration",ylab=expression(beta[5]*gamma[5]))
dev.off()

#Estimating posterior probability of top 5 gamma

library(tidyverse)

top5 = do.call(paste0, as.data.frame(Gamma)) %>%
  table() %>%
  sort(decreasing = TRUE) %>%
  .[1:5]
top5/10000

#Plotting posterior densities

pdf("Density105.pdf")
par(mfrow=c(3,2))

plot(density(Beta[,1]*Gamma[,1]))
plot(density(Beta[,2]*Gamma[,2]))
plot(density(Beta[,3]*Gamma[,3]))
plot(density(Beta[,4]*Gamma[,4]))
plot(density(Beta[,5]*Gamma[,5]))
dev.off()

#Posterior mean

mean_matrix <- matrix(c(mean(Beta[,1]*Gamma[,1]),mean(Beta[,2]*Gamma[,2]),mean(Beta[,3]*Gamma[,3]),mean(Beta[,4]*Gamma[,4]),mean(Beta[,5]*Gamma[,5])),ncol=5)
mean_matrix

prob_matrix <- matrix(c(mean(Gamma[,1]),mean(Gamma[,2]),mean(Gamma[,3]),mean(Gamma[,4]),mean(Gamma[,5])),ncol=5)
prob_matrix
