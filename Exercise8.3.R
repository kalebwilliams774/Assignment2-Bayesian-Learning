#Importing school math test data

school1 <- c(2.11, 9.75, 13.88, 11.3, 8.93, 15.66, 16.38, 4.54, 8.86, 11.94, 12.47, 11.11, 11.65, 14.53, 9.61, 7.38, 3.34, 9.06, 9.45, 5.98, 7.44, 8.5, 1.55, 11.45, 9.73)

school2 <- c(0.29, 1.13, 6.52, 11.72, 6.54, 5.63, 14.59, 11.74, 9.12, 9.43, 10.64, 12.28, 9.5, 0.63, 15.35, 5.31, 8.49, 3.04, 3.77, 6.22, 2.14, 6.58, 1.11)

school3 <- c(4.33, 7.77, 4.15, 5.64, 7.69, 5.04, 10.01, 13.43, 13.63, 9.9, 5.72, 5.16, 4.33, 12.9, 11.27, 6.05, 0.95, 6.02, 12.22, 12.85)

school4 <- c(12.46, 6.42, 5.96, 0.92, 11.43, 2.27, 1.54, 6.55, 2.3, 0.57, 7.4, 6.63, 7.02, 2.95, 4.44, 7.78, 8.36, 13.32, 8.81, 2.06, 14.17, 0.88, 10.36, 4.97)

school5 <- c(12.97, 13.6, 13.54, 5.49, 11.52, 8.23, 8.98, 6.42, 12.01, 15.08, 7.16, 10.84, 8.15, 4.27, 14.21, 15.93, 8.99, 10.12, 5.65, 14.94, 14.2, 8.43, 10.18, 17.47)

school6 <- c(2.5, 7.56, 5.79, 4.92, 3.32, 9.65, 2.58, 3.31, 5.47, 6.98, 9.74, 0.97, 6.2, 11.16, 13.45, 7.84, 10.43, 5.85, 5.56, 6.82, 5.23, 1.18)

school7 <- c(7.5, 11.15, 5.82, 0.39, 4.11, 4.82, 13.56, 3.11, 6.69, 7.33, 11.87, 9.14, 0.03, 1.76, 5.03, 3.72, 7.28, 7.15, 9.07, 8.59, 6.53, 0.27)

school8 <- c(6.41, 3.52, 7.65, 9.56, 9.49, 4.54, 14.72, 5.63, 4.24, 8.96, 8.59, 8.69, 6.18, 4.79, 11.67, 2.8, 7.03, 4.32, 11.51, 7.32)

# Combine all vectors into a single column
all_schools <- c(school1, school2, school3, school4, school5, school6, school7, school8)

# Display the combined vector
all_schools

#Creating vector for indices

indices <- c(rep(1,length(school1)), rep(2,length(school2)), rep(3,length(school3)), rep(4,length(school4)),rep(5,length(school5)), rep(6,length(school6)), rep(7,length(school7)),rep(8,length(school8)))

#Creating matrix 

Y <- matrix(c(indices,all_schools),ncol=2)
colnames(Y) <- c("school", "score")

### Weakly informative priors
nu0 <- 2 ; s20 <- 15
eta0 <- 2 ; t20 <- 10
mu0 <- 7 ; g20 <- 5

### Starting values
m <- length(unique(Y[, 1]))
n <- sv <- ybar <- rep(NA, m)

for (j in 1:m) {
  ybar[j] <- mean(Y[Y[, 1] == j, 2])
  sv[j] <- var(Y[Y[, 1] == j, 2])
  n[j] <- sum(Y[, 1] == j)
}

theta <- ybar ; sigma2 <- mean(sv)
mu <- mean(theta) ; tau2 <- var(theta)

### Setup MCMC
set.seed(1)
S <- 5000
THETA <- matrix(nrow = S, ncol = m)
MST <- matrix(nrow = S, ncol = 3)

### MCMC algorithm
for (s in 1:S) {
  # Sample new values of thetas
  for (j in 1:m) {
    v_theta <- 1 / (n[j] / sigma2 + 1 / tau2)
    e_theta <- v_theta * (ybar[j] * n[j] / sigma2 + mu / tau2)
    theta[j] <- rnorm(1, e_theta, sqrt(v_theta))
  }
  
  # Sample new value of sigma2
  nun <- nu0 + sum(n)
  ss <- nu0 * s20
  for (j in 1:m) {
    ss <- ss + sum((Y[Y[, 1] == j, 2] - theta[j])^2)
  }
  sigma2 <- 1 / rgamma(1, nun / 2, ss / 2)
  
  # Sample a new value of mu
  vmu <- 1 / (m / tau2 + 1 / g20)
  emu <- vmu * (m * mean(theta) / tau2 + mu0 / g20)
  mu <- rnorm(1, emu, sqrt(vmu))
  
  # Sample a new value of tau2
  etam <- eta0 + m
  ss <- eta0 * t20 + sum((theta - mu)^2)
  tau2 <- 1 / rgamma(1, etam / 2, ss / 2)
  
  # Store results
  THETA[s, ] <- theta
  MST[s, ] <- c(sigma2, mu, tau2)
}

# A plot for evaluating lack of convergence
stationarity.plot<-function(x,...){
  
  S<-length(x)
  scan<-1:S
  ng<-min( round(S/100),10)
  group<-S*ceiling(ng*scan/S)/ng
  
  boxplot(x~group,...)
}

#Autocorrelation plots to test convergence of mcmc

pdf("Stationarity83.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))

stationarity.plot(MST[,1],xlab="iteration",ylab=expression(mu))
stationarity.plot(MST[,2],xlab="iteration",ylab=expression(sigma^2))
stationarity.plot(MST[,3],xlab="iteration",ylab=expression(tau^2))
dev.off()


#Calculating effective sample size

library(coda)
effectiveSize(MST)

effectiveSize(THETA) -> esTHETA

MCERR<-  apply(MST,2,sd)/sqrt( effectiveSize(MST) )
apply(MST,2,mean)

#Posterior mean for sigma^2, mu and tau^2

mean((MST[,1]))
mean(sqrt(MST[,2]))
mean(sqrt(MST[,3]))

#95% confidence intervals for sigma^2, mu and tau^2

quantile(MST[,1],c(.025,.975))
quantile(MST[,2],c(.025,.975))
quantile(MST[,3],c(.025,.975))

# Priors for sigma^2, mu, and tau^2

library(invgamma)

seq1 = seq(0.1, 100, by = 0.1)

sigma_prior <- dinvgamma(seq1, nu0/2, nu0*s20/2)
mu_prior <- dnorm(seq1, mu0, sqrt(g20))
tau_prior <- dinvgamma(seq1, eta0/2, eta0*t20/2)


#Plots of posterior densities
pdf("PosteriorDensity83.pdf",family="Times",height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(density(MST[,1],adj=2),xlab=expression(sigma^2),main="",lwd=2,
     ylab=expression(paste(italic("p("),mu,"|",italic(y[1]),"...",italic(y[m]),")")),ylim = c(0, 0.35))
lines(seq1, sigma_prior,col='gray',lwd=2)
plot(density(MST[,2],adj=2),xlab=expression(mu),main="", lwd=2,
     ylab=expression(paste(italic("p("),sigma^2,"|",italic(y[1]),"...",italic(y[m]),")")))
lines(seq1,mu_prior,col='gray',lwd=2)
plot(density(MST[,3],adj=2),xlab=expression(tau^2),main="",lwd=2,
     ylab=expression(paste(italic("p("),tau^2,"|",italic(y[1]),"...",italic(y[m]),")")))
lines(seq1,tau_prior,col='gray',lwd=2)
legend('topright', lty = 1, col = c('black', 'gray'), legend = c('Posterior', 'Prior'))
dev.off()


#Plotting mean of posterior of theta and sample means

pdf('Means83.pdf')
plot(SA,theta.hat,xlab=expression(bar(italic(y))),ylab=expression(hat(theta)))
abline(0,1)
dev.off()

#Near linear relationship

#Calculating mean of all values and posterior mean of mu

mean(Y[,2])

mean(MST[,2])

#Posterior of R

R <- MST[,3]/(MST[,1]+MST[,3])

set.seed(123)
tau_p = 1/rgamma(S, eta0/2, eta0*t20/2)
s_p= 1/rgamma(S, nu0/2, nu0*s20/2)

R_prior <- tau_p/(s_p+tau_p)

pdf("Rpostpriorplot.pdf", family="Times", height=1.75, width=5)
par(mfrow=c(1,1), mar=c(2.75, 2.75, .5, .5), mgp=c(1.7, .7, 0))
plot(density(R), xlab=expression(tau^2/(sigma^2+tau^2)), ylab="Density",main='',ylim=c(0,10))
legend('topright', lty = 1, col = c('black', 'gray'), legend = c('Posterior', 'Prior'))
lines(density(R_prior), col='gray')
dev.off()

#Finding posterior probability that theta_7 is less than theta_6

mean(THETA[,7]<THETA[,6])

# Obtain the posterior probability that theta_7 is the smallest of all the theta's
mean(THETA[, 7] < THETA[,-7])


#Calculating sample averages for the data

SA <- c(mean(school1),mean(school2),mean(school3),mean(school4),
        mean(school5),mean(school6),mean(school7),mean(school8))

nj <- c(length(school1),length(school2),length(school3),length(school4),
        length(school5),length(school6),length(school7),length(school8))

theta.hat <- colMeans(THETA)



pdf("samplevsmeantheta.pdf", family="Times", height=1.75, width=5)
par(mfrow=c(1,2), mar=c(2.75, 2.75, .5, .5), mgp=c(1.7, .7, 0))
plot(density(SA), xlab=expression(bar(italic(y))~hat(theta)), ylab="Density",main='',ylim=c(0,0.35))
lines(density(theta.hat), col='gray')
plot(SA,theta.hat,xlab=expression(bar(y)),ylab =expression(hat(theta)))
abline(0,1)
dev.off()

#Calculating posterior mean of mu and of all data 

post_mu_mean <- mean(MST[,2])
tot_sample_mean <- mean(Y)

post_mu_mean
tot_sample_mean