# FIN251lab2.r			script file for econ 424 lab5 calculations
#
# author: Jane Yoo
# created: November 15, 2023
# revised: November 16, 2023
#
# comments:

options(digits=4, width=70)

# make sure packages are installed prior to loading them
library(PerformanceAnalytics)
library(zoo)
library(boot)
library(tseries)

# get monthly adjusted closing price data on VBLLX, FMAGX and "SBUX(YOUR STOCK)" from Yahoo
# using the tseries function get.hist.quote(). Set sample to Sept 2018 through
# Sep 2023. Note: if you are not careful with the start and end dates
# or if you set the retclass to "ts" then results might look weird

# get the last five years of monthly adjusted closing prices from Yahoo!
VBLLX.prices = get.hist.quote(instrument="vbllx", start="2018-09-01",
                             end="2023-09-30", quote="AdjClose",
                             provider="yahoo", origin="2018-01-01",
                             compression="m", retclass="zoo")
# change class of time index to yearmon which is appropriate for monthly data
# index() and as.yearmon() are functions in the zoo package 
#                             
index(VBLLX.prices) = as.yearmon(index(VBLLX.prices))
                             
class(VBLLX.prices)
colnames(VBLLX.prices)
start(VBLLX.prices)
end(VBLLX.prices)

FMAGX.prices = get.hist.quote(instrument="fmagx", start="2018-09-01",
                             end="2023-09-30", quote="AdjClose",
                             provider="yahoo", origin="2018-01-01",
                             compression="m", retclass="zoo")
index(FMAGX.prices) = as.yearmon(index(FMAGX.prices))

SBUX.prices = get.hist.quote(instrument="sbux", start="2018-09-01",
                             end="2023-09-30", quote="AdjClose",
                             provider="yahoo", origin="2018-01-01",
                             compression="m", retclass="zoo")
index(SBUX.prices) = as.yearmon(index(SBUX.prices))

# create merged price data
lab2Prices.z = merge(VBLLX.prices, FMAGX.prices, SBUX.prices)
# rename columns
colnames(lab2Prices.z) = c("VBLLX", "FMAGX", "SBUX")

# calculate cc returns as difference in log prices
lab2Returns.z = diff(log(lab2Prices.z))

#
# 3. Create timePlots of data
#

plot(lab2Returns.z, plot.type="single", lty=1:3, col=1:3, lwd=2)
legend(x="bottomleft", legend=colnames(lab2Returns.z), lty=1:3, col=1:3, lwd=2)
abline(h=0)
title("Monthly cc returns")


#
# 4. Create matrix of return data and compute pairwise scatterplots
#

ret.mat = coredata(lab2Returns.z)
colnames(ret.mat)
head(ret.mat)
VBLLX = ret.mat[,"VBLLX"]
FMAGX = ret.mat[,"FMAGX"]
SBUX = ret.mat[,"SBUX"]
pairs(ret.mat, col="blue")

#
# 5. Compute estimates of CER model parameters
#
muhat.vals = apply(ret.mat, 2, mean)
muhat.vals
sigma2hat.vals = apply(ret.mat, 2, var)
sigma2hat.vals
sigmahat.vals = apply(ret.mat, 2, sd)
sigmahat.vals
cov.mat = var(ret.mat)
cov.mat
cor.mat = cor(ret.mat)
cor.mat
covhat.vals = cov.mat[lower.tri(cov.mat)]
rhohat.vals = cor.mat[lower.tri(cor.mat)]
names(covhat.vals) <- names(rhohat.vals) <- 
c("VBLLX,FMAGX","VBLLX,SBUX","FMAGX,SBUX")
covhat.vals
rhohat.vals

# summarize the CER model estimates
cbind(muhat.vals,sigma2hat.vals,sigmahat.vals)
cbind(covhat.vals,rhohat.vals)

# plot mean vs. sd values
plot(sigmahat.vals, muhat.vals, pch=1:3, cex=2, col=1:3, 
     ylab = "mean", xlab="sd (risk)")
abline(h=0)     
legend(x="topright", legend=names(muhat.vals), pch=1:3, col=1:3, cex=1.5)     

#
# 6. Compute standard errors for estimated parameters
#

# compute estimated standard error for mean
nobs = nrow(ret.mat)
nobs
se.muhat = sigmahat.vals/sqrt(nobs)
se.muhat
# show estimates with SE values underneath
rbind(muhat.vals,se.muhat)

# compute estimated standard errors for variance and sd
se.sigma2hat = sigma2hat.vals/sqrt(nobs/2)
se.sigma2hat
se.sigmahat = sigmahat.vals/sqrt(2*nobs)
se.sigmahat

rbind(sigma2hat.vals,se.sigma2hat)
rbind(sigmahat.vals,se.sigmahat)

# compute estimated standard errors for correlation
se.rhohat = (1-rhohat.vals^2)/sqrt(nobs)
se.rhohat
rbind(rhohat.vals,se.rhohat)


#
# 7. Compute 5% and 1% Value at Risk
#

# function to compute Value-at-Risk
# note: default values are selected for 
# the probability level (p) and the initial
# wealth (w). These values can be changed
# when calling the function. Highlight the entire
# function, right click and select run line or selection
Value.at.Risk = function(x,p=0.05,w=100000) {
	x = as.matrix(x)
	q = apply(x, 2, mean) + apply(x, 2, sd)*qnorm(p)
	VaR = (exp(q) - 1)*w
	VaR
}

# 5% and 1% VaR estimates based on W0 = 100000

Value.at.Risk(ret.mat,p=0.05,w=100000)
Value.at.Risk(ret.mat,p=0.01,w=100000)



#
# 8. Evaluate bias and SE formulas using Monte Carlo
#

# generate 1000 samples from CER and compute sample statistics

mu = muhat.vals["FMAGX"]
sd = sigmahat.vals["FMAGX"]
n.obs = 60
set.seed(123)
n.sim = 1000

sim.means = rep(0,n.sim)
sim.vars = rep(0,n.sim)
sim.sds = rep(0,n.sim)
for (sim in 1:n.sim) {
	sim.ret = rnorm(n.obs,mean=mu,sd=sd)
	sim.means[sim] = mean(sim.ret)
	sim.vars[sim] = var(sim.ret)
	sim.sds[sim] = sqrt(sim.vars[sim])
}

par(mfrow=c(2,2))
hist(sim.means,xlab="mu hat", col="slateblue1")
abline(v=mu, col="white", lwd=2)
hist(sim.vars,xlab="sigma2 hat", col="slateblue1")
abline(v=sd^2, col="white", lwd=2)
hist(sim.sds,xlab="sigma hat", col="slateblue1")
abline(v=sd, col="white", lwd=2)
par(mfrow=c(1,1))

# 
# 9. compute MC estimates of bias and SE
#

c(mu, mean(sim.means))
mean(sim.means) - mu
c(sd^2, mean(sim.vars))
mean(sim.vars) - sd^2
c(sd, mean(sim.sds))
mean(sim.sds) - sd

# compute MC SE value and compare to SE calculated from simulated data

c(se.muhat["FMAGX"], sd(sim.means))
c(se.sigma2hat["FMAGX"], sd(sim.vars)) 
c(se.sigmahat["FMAGX"], sd(sim.sds))


#
# 10. bootstrapping SE for mean, variance, sd and correlation
#

?boot
# note: boot requires user-supplied functions that take
# two arguments: data and an index. The index is created
# by the boot function and represents random resampling with
# replacement

# function for bootstrapping sample mean
mean.boot = function(x, idx) {
# arguments:
# x 		data to be resampled
# idx		vector of scrambled indices created by boot() function
# value:
# ans		mean value computed using resampled data
     ans = mean(x[idx])
     ans
}

VBLLX.mean.boot = boot(VBLLX, statistic = mean.boot, R=999)
class(VBLLX.mean.boot)
names(VBLLX.mean.boot)

# print, plot and qqnorm methods
VBLLX.mean.boot
se.muhat["VBLLX"]

# plot bootstrap distribution and qq-plot against normal
plot(VBLLX.mean.boot)


# compute bootstrap confidence intervals from normal approximation
# basic bootstrap method and percentile intervals
boot.ci(VBLLX.mean.boot, conf = 0.95, type = c("norm","perc"))

#
# boostrap SD estimate
#
# function for bootstrapping sample standard deviation
sd.boot = function(x, idx) {
# arguments:
# x 		data to be resampled
# idx		vector of scrambled indices created by boot() function
# value:
# ans		sd value computed using resampled data
     ans = sd(x[idx])
     ans
}

VBLLX.sd.boot = boot(VBLLX, statistic = sd.boot, R=999)
VBLLX.sd.boot
se.sigmahat["VBLLX"]

# plot bootstrap distribution
plot(VBLLX.sd.boot)

# compute confidence intervals
boot.ci(VBLLX.sd.boot, conf=0.95, type=c("norm", "basic", "perc"))

# bootstrap correlation

# function to compute correlation between 1st 2 variables in matrix
rho.boot = function(x.mat, idx) {
# x.mat	n x 2 data matrix to be resampled
# idx		vector of scrambled indices created by boot() function
# value:
# ans		correlation value computed using resampled data

	ans = cor(x.mat[idx,])[1,2]
	ans
}
VBLLX.FMAGX.cor.boot = boot(ret.mat[,c("VBLLX","FMAGX")],
                           statistic=rho.boot, R = 999)
VBLLX.FMAGX.cor.boot
se.rhohat[1]

# plot bootstrap distribution
plot(VBLLX.FMAGX.cor.boot)

# bootstrap confidence intervals
boot.ci(VBLLX.FMAGX.cor.boot, conf=0.95, type=c("norm", "perc"))

#
# 11. Bootstrap VaR
#

# 5% Value-at-Risk
ValueAtRisk.boot = function(x, idx, p=0.05, w=100000) {
# x.mat	data to be resampled
# idx		vector of scrambled indices created by boot() function
# p		probability value for VaR calculation
# w		value of initial investment
# value:
# ans		Value-at-Risk computed using resampled data

	q = mean(x[idx]) + sd(x[idx])*qnorm(p)
	VaR = (exp(q) - 1)*w
	VaR
}

VBLLX.VaR.boot = boot(VBLLX, statistic = ValueAtRisk.boot, R=999)
VBLLX.VaR.boot
boot.ci(VBLLX.VaR.boot, conf=0.95, type=c("norm", "perc"))
plot(VBLLX.VaR.boot)




