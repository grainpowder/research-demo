# Demo code for semiparametric regression using data with measurement error
# Date   : 2020/09/30
# Author : Sunsik Kim

library(rstan)
source("./model_vb.R")


# Simulation data generation ----------------------------------------------
set.seed(10)
n_intknot = 5
N = 130
D = 6
RR = 0.9
xi2 = 0.8
sig2v = xi2/RR-xi2
mux = 1.5
beta = rnorm(D+1)
w = matrix(rnorm(N*D),N,D)
x = rnorm(N, mux, sqrt(xi2))
v = rnorm(N, x, sqrt(sig2v))
f = function(x) 2*x+2*sin(pi*x)
y = drop(cbind(1,w)%*%beta + f(x) + rnorm(N, sd=0.7))


# Create list of objects to pass into Stan --------------------------------
resolution = 200 # number of x's to estimate f(x)
order = 4 # To build cubic B-spline basis
knots = quantile(unique(v), seq(0,1,length=n_intknot+2)[-c(1,n_intknot+2)])
boundary = c(min(v)-sd(v)/2, max(v)+sd(v)/2)
xgrid = seq(boundary[1], boundary[2], length.out=resolution) 
bsgrid = bs(x=xgrid, knots=knots, intercept=TRUE, Boundary.knots=boundary)
data = list(
  n_obs = length(y),
  n_intknot = length(knots),
  n_regcoef = D,
  int_knots = knots,
  boundary = boundary,
  order = order,
  resolution = resolution,
  y = y,
  v = v,
  bsgrid = bsgrid,
  w = w)


# Estimate : Variational Bayes --------------------------------------------
start = Sys.time()
vb_result = model_vb(y,w,v,n_intknot=5)
Sys.time() - start


# Estimate : MCMC ---------------------------------------------------------
set.seed(100)
start = Sys.time()
mcmc_result = stan("./model_mcmc.stan", data=data, iter=1e3)
Sys.time() - start

write.csv(extract(mcmc_result, 'fxGrid')$fxGrid, "./samples_fxGrid.csv", row.names=FALSE)
write.csv(extract(mcmc_result, 'beta')$beta, "./samples_beta.csv", row.names=FALSE)
mcmc_curve = as.matrix(read.csv("./samples_fxGrid.csv", header=TRUE))
mcmc_beta = as.matrix(read.csv("./samples_beta.csv", header=TRUE))



# Result presentation : Parametric part -----------------------------------
# png('./result.png', width=900, height=250) # Uncomment when saving plots into file
# par(mfrow=c(1,4), mar=c(4,4,4,1)) # Uncomment when saving plots into file
deviance = 0.2
ord = order(beta[-1])
plot(1:length(beta[-1]), beta[-1][ord], pch=19, xlim=c(0.5,6.5), ylim=c(-2,1.2), xlab="", ylab="", main="Regression Coefficients")

vb_beta_mean = vb_result$mubeta.q[-1][ord]
vb_beta_var = diag(vb_result$sigbeta.q)[-1][ord]
mc_beta_mean = colMeans(mcmc_beta)[-1][ord]
mc_beta_upper = apply(mcmc_beta, 2, quantile, 0.975)[-1][ord]
mc_beta_lower = apply(mcmc_beta, 2, quantile, 0.025)[-1][ord]

for (idx in 1:D) {
  vb_upper = qnorm(0.975, vb_beta_mean[idx], sqrt(vb_beta_var[idx]))
  vb_lower = qnorm(0.025, vb_beta_mean[idx], sqrt(vb_beta_var[idx]))
  lines(c(idx-deviance, idx-deviance), c(vb_upper, vb_lower), lwd=3, col="darkgoldenrod4")
  lines(c(idx+deviance, idx+deviance), c(mc_beta_upper[idx], mc_beta_lower[idx]), lwd=3, col="deepskyblue4")
}
points((1:length(beta[-1]))-deviance, vb_result$mubeta.q[-1][ord], pch=19, col="darkgoldenrod1")
points((1:length(beta[-1]))+deviance, colMeans(mcmc_beta)[-1][ord], pch=19, col="deepskyblue1")

legend(0.5,1,col=c("darkgoldenrod4", "deepskyblue3"),lwd=c(4,4),legend=c("VB", "MCMC"),bty="n")



# Result presentation : Nonparametric part --------------------------------

plot(v, y, pch=19,cex=0.3, main="Contaminated Scatter Plot")

vb_x = c(vb_result$xgrid, rev(vb_result$xgrid))
vb_y = vb_result$mubeta.q[1] + c(vb_result$post_lower, rev(vb_result$post_upper))
mcmc_x = c(xgrid, rev(xgrid))
mcmc_y = colMeans(mcmc_beta)[1] + c(apply(mcmc_curve,2,quantile,0.025), rev(apply(mcmc_curve,2,quantile,0.975)))

plot(x,y-w%*%vb_result$mubeta.q[-1], main="Mean curve(VB)", ylab="Residual", cex=0.3, pch=19)
polygon(vb_x, vb_y, col="darkgoldenrod1", lty = "blank")
points(x,y-w%*%vb_result$mubeta.q[-1], cex=0.3, pch=19)
lines(vb_result$xgrid, f(vb_result$xgrid), lwd=3, col="grey60")
lines(vb_result$xgrid, vb_result$post_curve+vb_result$mubeta.q[1], col="darkgoldenrod4", lwd=3)

plot(x,y-w%*%colMeans(mcmc_beta)[-1], main="Mean Curve(MCMC)", ylab="Residual", cex=0.3, pch=19)
polygon(mcmc_x, mcmc_y, col="deepskyblue1", lty="blank")
points(x,y-w%*%colMeans(mcmc_beta)[-1],cex=0.3, pch=19)
lines(vb_result$xgrid, f(vb_result$xgrid), lwd=3, col="grey60")
lines(xgrid, colMeans(mcmc_curve)+colMeans(mcmc_beta)[1], col="deepskyblue4", lwd=4)
# dev.off() # Uncomment when saving plots into file