model_vb = function(y, w, v, n_intknot=10, prior=NULL, maxiter=500, tol=1e-3, n_grids=1e3, resolution=200)
{
  # Semiparametric measurement error regression using B-spline basis
  # Author : Sunsik Kim
  #
  # y : response
  # w : ordinary explanatory variable(s)
  # v : variate containing measurement error
  
  library(splines)
  N = length(y)
  D = ncol(w)
  W = cbind(1, w)
  WtW = crossprod(W)
  
  # Hyperparameters
  if (is.null(prior))
  {
    sig2mu = sig2beta = 100
    anu = axi = bxi = aups = bups = asig = bsig = 1e-3
    bnu = 1e-4
  }
  else
  {
    anu = prior$anu
    bnu = prior$bnu
    sig2mu = prior$sig2mu
    axi = prior$axi
    bxi = prior$bxi
    sig2beta = prior$sig2beta
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
  anutl = anu + (N/2)
  axitl = axi + (N/2)
  aupstl = aups + (n_knot/2)
  asigtl = asig + (N/2)
  nu.ratio = anu/bnu
  xi.ratio = axi/bxi
  ups.ratio = aups/bups
  sig.ratio = asig/bsig
  muu.q = rep(0, n_knot)
  sigu.q = diag(rep(ups.ratio, n_knot))
  mutl = mean(v)
  vphiq = matrix(0, N, n_knot)
  
  # Update as
  lb = rep(0, maxiter)
  lbold = -Inf
  for (iter in 1:maxiter)
  {
    # beta
    sigbeta.q = solve(diag(rep(1/sig2beta,D+1)) + sig.ratio*WtW)
    mubeta.q = drop(sig.ratio*sigbeta.q%*%t(W)%*%(y-vphiq%*%muu.q))
    
    # denoised values
    common = -0.5*(sig.ratio*diag(vphig%*%(outer(muu.q,muu.q)+sigu.q)%*%t(vphig)) + (xi.ratio+nu.ratio)*(grids^2) - 2*xi.ratio*mutl*grids)
    lnpgrids = common + nu.ratio*outer(grids, v) + sig.ratio*outer(drop(vphig%*%muu.q), drop(y-W%*%mubeta.q)); pgrids = exp(t(lnpgrids))
    normalizers = apply(pgrids, 1, sum)
    ex = drop(pgrids%*%grids)/normalizers
    ex2 = drop(pgrids%*%(grids^2))/normalizers
    varx = ex2 - (ex)^2
    vphiq = pgrids%*%vphig/outer(normalizers,rep(1,n_knot))
    vphiqtvphiq = t(vphig)%*%diag(apply(pgrids/normalizers,2,sum))%*%vphig
    
    # nu
    bnutl = bnu + 0.5*(sum((v-ex)^2) + sum(varx))
    nu.ratio = anutl/bnutl
    
    # mu
    sig2mutl = 1/(1/sig2mu + xi.ratio*N)
    mutl = xi.ratio*sig2mutl*sum(ex)
    
    # xi
    bxitl = bxi + 0.5*(sum((ex-mutl)^2) + sum(varx) + N*sig2mutl)
    xi.ratio = axitl/bxitl
    
    # spline coefficients
    sigu.q = solve(diag(rep(ups.ratio,n_knot)) + sig.ratio*vphiqtvphiq)
    muu.q = drop(sig.ratio*sigu.q%*%t(vphiq)%*%(y-W%*%mubeta.q))
    
    # upsilon
    bupstl = bups + 0.5*sum(muu.q^2+diag(sigu.q))
    ups.ratio = aupstl/bupstl
    
    # sigma
    cpterm = sum((y-W%*%mubeta.q)^2) + sum(diag(WtW%*%sigbeta.q))
    cpterm = cpterm - 2*sum((y-W%*%mubeta.q)*(vphiq%*%muu.q)) + sum(diag(vphiqtvphiq%*%(sigu.q+outer(muu.q,muu.q))))
    bsigtl = bsig + 0.5*cpterm
    sig.ratio = asigtl/bsigtl
    
    # ELBO
    lbnew = -0.5*N*log(2*pi) - 0.5*N*(log(bsigtl)-digamma(asigtl)) - 0.5*sig.ratio*cpterm
    lbnew = lbnew - 0.5*N*log(2*pi) - 0.5*N*(log(bnutl)-digamma(anutl)) - 0.5*nu.ratio*(sum((v-ex)^2)+sum(varx))
    lbnew = lbnew - 0.5*N*(log(bxitl)-digamma(axitl)) - 0.5*xi.ratio*(sum((ex-mutl)^2)+sum(varx)+N*sig2mutl)
    lbnew = lbnew + 0.5*sum(log(varx)) + 0.5*N
    lbnew = lbnew - 0.5*(D+1)*log(sig2beta) - 0.5*sum(mubeta.q^2+diag(sigbeta.q))/sig2beta
    lbnew = lbnew + 0.5*determinant(sigbeta.q,logarithm=TRUE)$modulus[1] + 0.5*(D+1)
    lbnew = lbnew - 0.5*n_knot*(log(bupstl)-digamma(aupstl)) - 0.5*ups.ratio*sum(muu.q^2+diag(sigu.q))
    lbnew = lbnew + 0.5*determinant(sigu.q,logarithm=TRUE)$modulus[1] + 0.5*n_knot
    lbnew = lbnew - 0.5*log(sig2mu) - 0.5*(mutl^2+sig2mutl)/sig2mu
    lbnew = lbnew + 0.5*log(sig2mutl) + 0.5
    lbnew = lbnew - lgamma(axi) + axi*log(bxi) - (axi+1)*(log(bxitl)-digamma(axitl)) - xi.ratio*bxi
    lbnew = lbnew + lgamma(axitl) - axitl*log(bxitl) + (axitl+1)*(log(bxitl)-digamma(axitl)) + xi.ratio*bxitl
    lbnew = lbnew - lgamma(anu) + anu*log(bnu) - (anu+1)*(log(bnutl)-digamma(anutl)) - nu.ratio*bnu
    lbnew = lbnew + lgamma(anutl) - anutl*log(bnutl) + (anutl+1)*(log(bnutl)-digamma(anutl)) + nu.ratio*bnutl
    lbnew = lbnew - lgamma(asig) + asig*log(bsig) - (asig+1)*(log(bsigtl)-digamma(asigtl)) - sig.ratio*bsig
    lbnew = lbnew + lgamma(asigtl) - asigtl*log(bsigtl) + (asigtl+1)*(log(bsigtl)-digamma(asigtl)) + sig.ratio*bsigtl
    lbnew = lbnew - lgamma(aups) + aups*log(bups) - (aups+1)*(log(bupstl)-digamma(aupstl)) - ups.ratio*bups
    lbnew = lbnew + lgamma(aupstl) - aupstl*log(bupstl) + (aupstl+1)*(log(bupstl)-digamma(aupstl)) + ups.ratio*bupstl
    lb[iter] = lbnew
    if (abs(lbnew-lbold)<tol) break
    lbold = lbnew
  }
  lb = lb[1:iter]
  xgrid = seq(min(v)-sd(v)/2, max(v)+sd(v)/2, length.out=resolution) 
  vphi = bs(x=xgrid, knots=intKnots, intercept=TRUE, Boundary.knots=boundary)
  post_curve=drop(vphi%*%muu.q)
  return(list(
    lb=lb, ex=ex, varx=varx, mubeta.q=mubeta.q, sigbeta.q=sigbeta.q, muu.q=muu.q, sigu.q=sigu.q,
    sig.ratio=sig.ratio, nu.ratio=nu.ratio, xi.ratio=xi.ratio, ups.ratio=ups.ratio,
    xgrid=xgrid,
    post_curve=post_curve,
    post_lower=vphi%*%qnorm(0.025,muu.q,sqrt(diag(sigu.q))),
    post_upper=vphi%*%qnorm(0.975,muu.q,sqrt(diag(sigu.q)))
  ))
}