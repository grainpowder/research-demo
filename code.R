nonparametric_mer = function(y, v, n_intknot=10, prior=NULL, maxiter=500, tol=1e-4, n_grids=1e3, resolution=200)
{
  # Nonparametric Regression Model with Measurement Error
  # Model : y_i = f(x_i) + e_i
  # Input
  #   y       : response variable
  #   v       : contaminated explanatory variable
  
  library(splines)
  N = length(y)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sig2mu = 100
    anu = bnu = axi = bxi = aups = bups = asig = bsig = 1e-3
  }
  else
  {
    sig2mu = prior$sig2mu
    anu = prior$anu
    bnu = prior$bnu
    axi = prior$axi
    bxi = prior$bxi
    aups = prior$aups
    bups = prior$bups
    asig = prior$asig
    bsig = prior$bsig
  }
  
  # Compose interior knots of the basis
  intKnots = quantile(unique(v), seq(0,1,length=n_intknot+2)[-c(1,n_intknot+2)])
  boundary = c(min(v)-sd(v)/2, max(v)+sd(v)/2)
  
  # Compose spline basis of grids
  grids = seq(min(v)-sd(v)/2, max(v)+sd(v)/2, length.out=n_grids)
  vphig = bs(x=grids, knots=intKnots, intercept=TRUE, Boundary.knots=boundary)
  n_knot = ncol(vphig)
  
  # Initialize variational parameters
  asigtl = asig + (N/2)
  anutl = anu + (N/2)
  axitl = axi + (N/2)
  aupstl = aups + (n_knot/2)
  sig.ratio = asig/bsig
  nu.ratio = anu/bnu
  xi.ratio = axi/bxi
  ups.ratio = aups/bups
  muu.q = rep(0, n_knot)
  sigu.q = diag(rep(ups.ratio, n_knot))
  mutl = mean(v)
  
  # Update as
  for (iter in 1:maxiter)
  {
    sig.ratio.old = sig.ratio
    nu.ratio.old = nu.ratio
    xi.ratio.old = xi.ratio
    ups.ratio.old = ups.ratio
    
    # denoised values
    common = -0.5*(sig.ratio*diag(vphig%*%(outer(muu.q,muu.q)+sigu.q)%*%t(vphig)) + (xi.ratio+nu.ratio)*(grids^2) - 2*xi.ratio*mutl*grids)
    lnpgrids = common + nu.ratio*outer(grids, v) + sig.ratio*outer(drop(vphig%*%muu.q), drop(y))
    pgrids = exp(t(lnpgrids))
    normalizers = apply(pgrids, 1, sum)
    ex = drop(pgrids%*%grids)/normalizers
    ex2 = drop(pgrids%*%(grids^2))/normalizers
    varx = ex2 - (ex)^2
    vphiq = pgrids%*%vphig/outer(normalizers,rep(1,n_knot))
    vphiqtvphiq = t(vphig)%*%diag(apply(pgrids/normalizers,2,sum))%*%vphig
    
    # nu
    bnutl = bnu + 0.5*sum((v-ex)^2 + varx)
    nu.ratio = anutl/bnutl
    
    # mu
    sig2mutl = 1/(1/sig2mu + N*xi.ratio)
    mutl = xi.ratio*sig2mutl*sum(ex)
    
    # xi
    bxitl = bxi + 0.5*(sum((ex-mutl)^2) + sum(varx) + N*sig2mutl)
    xi.ratio = axitl/bxitl
    
    # spline coefficient
    sigu.q = solve(diag(rep(ups.ratio,n_knot)) + sig.ratio*vphiqtvphiq)
    muu.q = drop(sig.ratio*sigu.q%*%t(vphiq)%*%y)
    
    # upsilon
    bupstl = bups + 0.5*sum(muu.q^2+diag(sigu.q))
    ups.ratio = aupstl/bupstl
    
    # sigma
    cpterm = sum(y^2) - 2*sum(y*(vphiq%*%muu.q)) + sum(diag(vphiqtvphiq%*%(sigu.q+outer(muu.q,muu.q))))
    bsigtl = bsig + 0.5*cpterm
    sig.ratio = asigtl/bsigtl
    
    # Convergence check
    bool1 = abs(nu.ratio.old-nu.ratio) < tol
    bool2 = abs(xi.ratio.old-xi.ratio) < tol
    bool3 = abs(ups.ratio.old-ups.ratio) < tol
    bool4 = abs(sig.ratio.old-sig.ratio) < tol
    if (bool1&bool2&bool3&bool4) break
  }
  print(paste('Number of Iteration :', iter))
  if (iter < maxiter) print('Algorithm Converged')
  
  xgrid = seq(min(v)-sd(v)/2, max(v)+sd(v)/2, length.out=resolution) 
  vphi = bs(x=xgrid, knots=intKnots, intercept=TRUE, Boundary.knots=boundary)
  post_curve=drop(vphi%*%muu.q)
  return(list(
    xgrid=xgrid, ex=ex, varx=varx, muu.q=muu.q, sigu.q=sigu.q,
    sig.ratio=sig.ratio, nu.ratio=nu.ratio, xi.ratio=xi.ratio, ups.ratio=ups.ratio,
    post_curve=post_curve,
    post_lower=vphi%*%qnorm(0.025,muu.q,sqrt(diag(sigu.q))),
    post_upper=vphi%*%qnorm(0.975,muu.q,sqrt(diag(sigu.q)))
  ))
}

# Data generation
set.seed(10)
N = 130
RR = 0.7
xi2 = 0.8
sig2v = xi2/RR-xi2
mux = 1.5
x = rnorm(N, mux, sqrt(xi2))
v = rnorm(N, x, sqrt(sig2v))
f = function(x) 2*x+2*sin(pi*x)
y = f(x) + rnorm(N, sd=0.7)

# Apply Algorithm
vb_result = nonparametric_mer(y, v)

# Result Presentation
coverage = 0
for (idx in 1:N)
{
  lb = qnorm(0.025, vb_result$ex[idx], sqrt(vb_result$varx[idx]))
  ub = qnorm(0.975, vb_result$ex[idx], sqrt(vb_result$varx[idx]))
  if ((x[idx] > lb) & (x[idx] < ub)) coverage = coverage + 1
}
coverage = round(coverage/N*100, 2)
vb_x = c(vb_result$xgrid, rev(vb_result$xgrid))
vb_y = c(vb_result$post_lower, rev(vb_result$post_upper))

png('./result.png', width=900, height=320)
par(mfrow=c(1,3), mar=c(5,4.2,4,0.5))
plot(v, y, main='Contaminated Scatter Plot', pch=19, cex=0.5, xlab='v (Contaminated Data)', ylab='y (Response)')
plot(x, vb_result$ex, main='Denoised Values', pch=19, cex=0.5, xlab='x (True Data)', ylab=expression(paste(E[q](x), '(Estimated)')), sub=paste0('Coverage : ',coverage,'%'))
lines(-10:10,-10:10, lwd=2, col=2)
points(x, vb_result$ex, pch=19, cex=0.5)
plot(x, y, pch=19, cex=0.5, main='Mean curve and 95% credible set', xlab='x (True Data)', ylab='y (Response)')
polygon(vb_x, vb_y, col="darkgoldenrod1", lty = "blank")
points(x, y, pch=19, cex=0.5)
lines(vb_result$xgrid, f(vb_result$xgrid), col="deepskyblue3", lwd=4)
lines(vb_result$xgrid, vb_result$post_curve, col="darkgoldenrod4", lwd=4)
legend(-0.5,8,col=c("darkgoldenrod4", "deepskyblue3"),lwd=c(4,4),legend=c("Estimated", "True"),bty="n", seg.len=1)
dev.off()