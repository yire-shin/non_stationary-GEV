library(lmomco)
library(nleqslv)
library(robustbase)

#' Estimate the nonstationary GEV parameters
#'
#' @description
#' This method combines L-moments, a goodness-of-fit measure, and robust regression,
#' using standardized residuals obtained through a transformation to a stationary sequence.
#' The function fits a nonstationary Generalized Extreme Value (GEV) distribution to the given data,
#' which is particularly useful in hydrology and climatology for modeling extreme events.
#' For more details, refer to the paper "L-moment-based algorithm for nonstationary generalized extreme value model (2024)".
#'
#' @param xdat A numeric vector of data to be fitted. This is typically a series of annual maximums,
#' such as peak river flows or maximum daily temperatures.
#' @param ntry The number of attempts to perform in the optimization process. This parameter allows
#' the function to explore multiple starting points to avoid local minima and find a more global solution.
#' The default value is 20.
#' @param ftol The function tolerance for the optimization process. It determines the precision of the
#' optimization algorithm. A smaller value results in a more precise but potentially more time-consuming
#' optimization process. The default value is 1e-6.
#'
#' @return This function returns a list containing the following elements:
#' \itemize{
#'   \item \code{para.prop}: The proposed parameters for the nonstationary GEV model, including location,
#'   scale, and shape parameters, along with any nonstationarity coefficients.
#'   \item \code{precis}: The precision of the optimization process, indicating how closely the optimization
#'   criteria were met. A smaller value indicates a better fit.
#' }
#'
#' @examples
#' data("Trehafod")
#' result1 <- nsgev(Trehafod$r1)
#' print(result1)
#'
#' @author Jeong-Soo Park
#'
#' @export
nsgev <- function(xdat, ntry=20, ftol=1e-6){
  model='gev11'
  ns=length(xdat)
  year=seq(1,length(xdat))

  reg.dat=data.frame( cbind(year, xdat) )

  mu.init= lmrob(xdat~year, reg.dat)$coefficients  # robust regression
  m0.rob= mu.init[1]
  m1.rob= mu.init[2]

  orig.para=c(m0.rob, m1.rob, 1.0,-0.001,0)

  qlist= make.qmax_11(xdat, orig.para=orig.para, rob=T)

  orig.para=c(m0.rob, m1.rob, qlist$sig0, qlist$sig1, 0)

  gado.rob = time.m_11(qmax=qlist$qmax, orig.para=orig.para)

  z = multi.m0s0_11(xdat, ntry=ntry, ftol=ftol,
                    pretheta=gado.rob$para.org, model=model)  #proposed method
  # z$para.prop = proposed est
  if(z$precis > ftol) { z$para.prop = gado.rob$para.org
  cat("no optim for proposed","\n") }

  return(z)
}

#' Estimate the nonstationary GEV parameters
#'
#' @description
#' This function estimates the parameters of a nonstationary Generalized Extreme Value (GEV) distribution
#' for a given dataset. It applies a combination of techniques including Weighted Least Squares (WLS),
#' robust regression, and the GN16 method, to handle both stationary and nonstationary data. The function
#' is designed to provide a robust estimation of GEV parameters, particularly useful in fields like hydrology
#' and climatology where extreme value analysis is critical. For more details on the methods mentioned above, refer to
#' the paper "L-moment-based algorithm for nonstationary generalized extreme value model (2024)".
#'
#' @param xdat A numeric vector of data to be fitted. This data typically represents annual maximum values
#' of environmental variables like precipitation or river flow.
#' @param ntry The number of attempts for the optimization routine to avoid local minima and ensure a more
#' reliable estimation. The default is set to 20.
#' @param ftol The tolerance level for the optimization function, controlling the precision of the estimation.
#' The default value is set to 1e-6.
#'
#' @return A list containing the following elements:
#' \itemize{
#'  \item \code{para.prop}: Parameters estimated by the proposed method. For more details on the proposed methods, refer to
#' the paper "L-moment-based algorithm for nonstationary generalized extreme value model (2024)".
#'  \item \code{para.gado}: Parameters estimated by the GN16 method.
#'  \item \code{strup.sta}: Stationary parameters estimated by the Strup WLS method.
#'  \item \code{strup.org}: Original nonstationary parameters estimated by the Strup WLS method.
#'  \item \code{strup.final}: Final nonstationary parameters adjusted by the specified WLS method.
#'  \item \code{lme.sta}: Stationary parameters estimated by the L-moments method.
#' }
#' Each of these elements is a named vector of the GEV parameters, providing a comprehensive set of estimations
#' for further analysis.
#'
#'
#' @examples
#' data("Trehafod")
#' result2 <- gado.prop_11(Trehafod$r1)
#' print(result2)
#'
#' @author Jeong-Soo Park
#'
#' @export

gado.prop_11= function (xdat, ntry=20, ftol=1e-6){

  z <- list()
  ns=length(xdat)
  year=seq(1,length(xdat))

  model='gev11'
  name_gev11_ns  =c("mu0","mu1","sigma0","sigma1","xi")
  name_gev00_sta =c("mu_sta","sigma_sta","xi_sta")

# --------------Strup WLS ---------------------------------------------
  reg.dat=data.frame( cbind(year, xdat) )

  mu.init= lm(xdat~year, reg.dat)$coefficients      # regression
  m0= mu.init[1]
  m1= mu.init[2]

  orig.para=c(m0,m1,1.0,-0.001,0)

  strup = strup_11(xdat, orig.para=orig.para)        # strup wls

#---------------- GN16 ------------------------------------------------
  qlist= make.qmax_11(xdat, orig.para=orig.para, rob=F)

  orig.para=c(m0,m1, qlist$sig0, qlist$sig1, 0)

  gado = time.m_11(qmax=qlist$qmax, orig.para=orig.para)  # GN16

#------ proposed method --------------------------------------------
  mu.init= lmrob(xdat~year, reg.dat)$coefficients  # robust regression
  m0.rob= mu.init[1]
  m1.rob= mu.init[2]

  orig.para=c(m0.rob, m1.rob, 1.0,-0.001,0)

  qlist= make.qmax_11(xdat, orig.para=orig.para, rob=T)

  orig.para=c(m0.rob, m1.rob, qlist$sig0, qlist$sig1, 0)

  gado.rob = time.m_11(qmax=qlist$qmax, orig.para=orig.para)

  z = multi.m0s0_11(xdat, ntry=ntry, ftol=ftol,
                    pretheta=gado.rob$para.org, model=model)  #proposed method
                                              # z$para.prop = proposed est
  if(z$precis > ftol) { z$para.prop = gado.rob$para.org
  cat("no optim for proposed","\n") }

# --------------------------------------------------------------------
  z$para.gado   = gado$para.org            # GN16 orginal est
  z$strup.sta   = strup$strup.sta          # stationary wlse
  z$strup.org   = strup$strup.para         # wlse by strup
  z$strup.final = strup$strup.final        # specified WLSE

  z$lme.sta = pargev(lmoms(xdat,nmom=5))$para   # stationary L-ME

  names(z$para.prop)     <-name_gev11_ns
  names(z$para.gado)     <-name_gev11_ns
  names(z$strup.org)     <-name_gev11_ns
  names(z$strup.final)   <-name_gev11_ns
  names(z$strup.sta)     <-name_gev00_sta
  names(z$lme.sta)       <-name_gev00_sta

  return(z)
}

#-----strup wls -----------------------------------
strup_11 =function(xdat, orig.para=orig.para){

  w=list()
  ns=length(xdat)
  year=seq(1,ns)

  m0=orig.para[1]
  m1=orig.para[2]
  res=xdat -(m0+m1*year)

    stand=wls.park_11(xdat,res)              # regression and Steps 3, 4, and 5

    new.para=c(stand$m, stand$sig, 0)
    ares=stand$res
    sigt =  exp(new.para[3]+new.para[4]*year)

  w$strup.sta = pargev(lmoms(ares))$para
  w$strup.para = c(new.para[1:4], w$strup.sta[3])

#------to specify final parameter values -----------------------------
  yt=rep(0,ns+1)
  mu_st=w$strup.sta[1]
  year2=seq(0,ns)

  for (k in 1:(ns+1) ) {
     t=k-1
     yt[k]=  w$strup.para[1] + w$strup.para[2] * t
     yt[k]= yt[k] + mu_st * exp(w$strup.para[3]+ w$strup.para[4]* t)
  }

  nh=round((ns/2))
  yt[nh-1] = yt[nh-1] + 0.03;   yt[nh-2] = yt[nh-2] - 0.02
  yt[nh+1] = yt[nh+1] - 0.03;   yt[nh+2] = yt[nh+2] + 0.02

  reg.dat=data.frame( cbind(year2, yt) )

  mu.init= lm(yt~year2, reg.dat)$coefficients

  sigmaf_0 = w$strup.para[3] + log(w$strup.sta[2])  #w$strup.sta[2] = sig_st
  sigmaf_1 = w$strup.para[4]
  xif = w$strup.para[5]

  w$strup.final= c(mu.init, sigmaf_0, sigmaf_1, xif)
  return(w)
}
#----------------------------------------------------------
wls.park_11 = function(xdat,res ){

  z=list()
  ns=length(res)
  year=seq(1,ns)

  lres.pr=log(abs(res))
  sig.dat=data.frame( cbind(year, lres.pr) )
  z$sig= lm(lres.pr~year, sig.dat)$coefficients             # Step 3

  sigt =  exp(z$sig[1] + z$sig[2]* year)

  res.n = xdat/sigt
  ytran0 = rep(1,ns)/sigt
  ytran1 = year/sigt

  new.data= data.frame( cbind(ytran0, ytran1, res.n) )

  z$m= lm(res.n~0+ytran0+ytran1, new.data)$coefficients      # Step 4

  z$res = res.n -(z$m[1]*ytran0+z$m[2]*ytran1)               # Step 5
  return(z)
}
#------------------------------------------------------
#' @export False
time.m_11 = function(qmax=NULL, orig.para=NULL){

  z=list()
  para.gado=rep(NA,5)
  ns=length(qmax)
  year=seq(1,ns)

  m0=orig.para[1]
  m1=orig.para[2]
  sig0=orig.para[3]
  sig1=orig.para[4]

    lmom_q = lmoms(qmax)
    q.sta = pargev(lmom_q)$para
    xi = q.sta[3]
    if( xi <= -0.5) xi = -0.4999

  cd = sqrt( (xi^2) /( gamma(1+2*xi) - gamma(1+xi)^2 ) )
  alpha_t = exp(sig0 +sig1*year) * cd
  z$mu_t=  - (1-gamma(1+ xi) )*alpha_t/xi  + m0+ m1*year

  #--------------------------------------------------------
  nh=round((ns/2))
  mu.gado=z$mu_t
  mu.gado[nh-1] = mu.gado[nh-1] + 0.02
  mu.gado[nh+1] = mu.gado[nh+1] - 0.02

  mu.data= data.frame( cbind(year, mu.gado) )
  loc.gado =lm(mu.gado~year, mu.data)$coefficients

  alpha0 = log(cd) + sig0
  alpha1 = sig1

  z$para.org= c(loc.gado, alpha0, alpha1, xi )

# modify GN16 ---------------------------------------------
#  alpha0_up = log(q.sta[2])- sig1*nh
#  z$para.up = c(loc.gado, alpha0_up, alpha1, xi )

#  alpha_t = exp(alpha0_up +sig1*year)
#  z$mu_t.up=  - (1-gamma(1+ xi) )*alpha_t/xi  + m0+ m1*year
# ---------------------------------------------------------
  return(z)
}

#-------------------------------------------------
gev.Ldist.m0s0_11 <- function(a, xdat=xdat, pretheta=pretheta)
{
  zz=rep(100,3)

  mu0 <- a[1]
  mu1 <- pretheta[2]
  sig0 <- a[2]
  sig1 = pretheta[4]
  xi= a[3]

  ns=length(xdat)
  year=seq(1,ns)
  gum01=rep(NA, ns)
  gum.dat=rep(NA, ns)
  newg2=rep(NA,ns)

  gum.dat[1:ns]= xdat[1:ns]-(mu0 + mu1*year[1:ns])
  gum.dat[1:ns]= gum.dat[1:ns]/exp(sig0 + sig1*year[1:ns])

  gum01[1:ns]= 1-xi*gum.dat[1:ns]   #^(1/(-xi))

  for (it in 1:ns ) {
    if( is.na(gum01[it]) ){
      newg2[it]=NA

    }else if( gum01[it] <= 0 ) {
      newg2[it]= NA
    }else if( gum01[it] > 0) {
      newg2[it]= log(gum01[it])/(-xi)
    }
  }

  newg=newg2[!is.na(newg2)]
  newg= newg*20 +300

  if( length(newg) < ns/2 ) {
    #      cat("Too many NA in newg","\n")
    zz[1:3]=1000
    return(zz)
  }

  lam=list()
#  if(mtrim==F){
    lam= lmomgum(vec2par(c(0,1),'gum'))
#  }else if(mtrim==T){
#    lam= theoTLmoms(vec2par(c(0,1),'gum'), leftrim=5)
#  }

  lgum=list()
  lgum=lmoms.md.park(newg, mtrim=F, no.stop=T)

  if(lgum$ifail == 1) {
    zz[1:3]=1000
    return(zz)
  }

  pen= max(abs(xi)- 1.0, 0)

  zz[1] = lam$lambdas[1]*20 + 300 - lgum$lambdas[1]
  zz[2] = lam$lambdas[2]*20  - lgum$lambdas[2]
  zz[3] = lam$ratios[3]  - lgum$ratios[3]

  zz[3] = zz[3] + sign(zz[3])*pen
  return(zz)
}
#---------------------------------------------------
multi.m0s0_11= function(xdat, ntry=20, ftol=1e-6,
                        pretheta=pretheta, model=model)
{
  z=list()
  value=list()
  k=list()

  init = matrix(0, nrow=ntry, ncol=3)
  init = ginit.m0s0(xdat, ntry, pretheta)

  if(model=='gev10') npar=4
  if(model=='gev20') npar=5
  if(model=='gev11') npar=5

  precis=rep(1000, ntry)
  para.sel=matrix(NA,ntry+1,ncol=npar)

  tryCatch({
    for(i in 1:ntry) {
      value =  tryCatch( nleqslv( x=as.vector(init[i,1:3]),
                                  fn= gev.Ldist.m0s0_11,
                                  method="Broyden",
                                  xdat=xdat, pretheta=pretheta ) )
      k[[i]] <- value

      if(is(value)[1]=="try-error"){
        k[[i]]$fvec <- 10^6
        k[[i]]$termcd = 5
      }else{
        precis[i]=  mean(abs(k[[i]]$fvec) )

        if( precis[i] < ftol) {
          k[[i]]$root = value$x
          para.sel[i,1:5]=c( k[[i]]$root[1], pretheta[2],
                             k[[i]]$root[2], pretheta[4],  k[[i]]$root[3])
        }
      }                           # ks[i]= KS_gev110(xdat, para)$KS_ST[1]

      precis[is.na(precis[i])]=1000
      if( abs( k[[i]]$termcd ) > 3 ) {
        precis[i]=1000
        para.sel[i,]=NA
      }
    } #end for
  }) # trycathch

  z$para.prop =sel.para_all(xdat, para.sel, model)$para

  z$precis =precis[which.min(precis)]
  return(z)
}
#-------------------------------------------------
sel.para_all =function(xdat, para.sel, model=NULL){

  z=list()
  upara.sel = para.sel

  gof=rep(NA,nrow(upara.sel))
  ns=length(xdat)
  npar=ncol(upara.sel)         # if(model=='gev10') npar=4
                               # if(model=='gev20') npar=5
                               # if(model=='gev11') npar=5
  vecT=c(5,10,20,40,60)
  if(ns >= 100) vecT=c(5,10,20,40,80,120)
  if(ns <= 30) vecT=c(5,10,20,40)

    for(i in 1:nrow(upara.sel) ){
      par.vec = as.vector(upara.sel[i,1:npar])
      gof[i]  = gof.ene_all(xdat, vecT, par.vec, model)
    }

  if(length(unique(gof))==1) {
    z$para =upara.sel[length(gof),]
  }else{
    z$para =upara.sel[which.min(gof),]
  }

  z$gof=gof
  return(z)
}
#-------------------------------------------------
gof.ene_all = function(xdat, vecT=c(5,10,20,40,80), para, model=NULL){

  ns=length(xdat)
  nT = length(vecT)
  year=seq(1,ns)
  chi=rep(NA,nT)

  for(i in 1:nT){
    T =  vecT[i]
    qt = qns.gev_all(T, para, year, model)
    ene = ns/T
    sne = sum(xdat >= qt)
    chi[i] = abs(ene-sne) /ene
  }
  chi2 = sum(chi)
  return(chi2)
}

#--------------------------------------------------
lmoms.md.park =
  function (x, nmom = 5, mtrim=F, no.stop = FALSE, vecit = FALSE)
  {
    z=list()
    ifail=0

    n <- length(x)
    if (nmom > n) {
      if (no.stop) {
        ifail=1
        z$ifail=ifail
        return(z)
      }else{
        stop("More L-moments requested by parameter 'nmom' than data points available in 'x'")
      }
    }
    if (length(unique(x)) == 1) {
      if (no.stop) {
        ifail=1
        z$ifail=ifail
        return(z)
      }else{stop("all values are equal--Lmoments can not be computed")
      }
    }

    if(ifail == 0) {

      if(mtrim==F){
        z <- TLmoms(x, nmom = nmom)
      }else if(mtrim==T){
        z <- TLmoms(x, nmom = nmom, leftrim=5)
      }
      z$source <- "lmoms"
      if (!vecit)
        z$ifail=ifail
      return(z)
      if (nmom == 1) {
        z <- z$lambdas[1]
      }
      else if (nmom == 2) {
        z <- c(z$lambdas[1], z$lambdas[2])
      }
      else {
        z <- c(z$lambdas[1], z$lambdas[2], z$ratios[3:nmom])
      }
      attr(z, which = "trim") <- NULL
      attr(z, which = "rightrim") <- NULL
      attr(z, which = "leftrim") <- NULL
      attr(z, which = "source") <- "lmoms"

    }
    z$ifail=ifail
    return(z)
  }
#-----------------------------------------------------
make.qmax_11 =function(xdat=NULL, orig.para=NULL, rob=NULL)
{

  z=list()
  m0= orig.para[1]
  m1= orig.para[2]
  ns = length(xdat)
  year= seq(1,ns)

  res = xdat - (m0 + m1* year)
  mres= mean(res)
  res.pr = abs(res - mres)

  lres.pr=log(res.pr)
  sig.dat=data.frame( cbind(year, lres.pr) )

  if(rob==F){
     sig.lm= lm(lres.pr~year, sig.dat)$coefficients
  }else if(rob==T){
     sig.lm= lmrob(lres.pr~year, sig.dat)$coefficients
  }

  sig0 = sig.lm[1]
  sig1 = sig.lm[2]
  sigt = exp(sig0 + sig1* year)

  qmax= rep(NA, ns)

  for(i in 1:ns){
    if(sig1 >= 0){
      if(res[i] >= mres ) {
        qmax[i] = res[i] - sigt[i]
      }else{
        qmax[i] = res[i] + sigt[i]
      }
    }else{
      if(res[i] >= mres ) {
        qmax[i] = res[i] + sigt[i]
      }else{
        qmax[i] = res[i] - sigt[i]
      }
    }
  }
  z$qmax=qmax
  z$sig0=sig0
  z$sig1=sig1
  return(z)
}
#------------------------------------------------------
ginit.m0s0 <-function(data, ntry=ntry, pretheta){

  init <-matrix(0, nrow=ntry, ncol=3)
  if(abs(pretheta[5]) > 0.5) pretheta[5] = sign(pretheta[5])*0.48

  lmom_init = lmoms(data,nmom=5)
  lmom_est <- pargev(lmom_init)

  init[1,1] <- lmom_est$para[1]
  init[1,2] <- log(lmom_est$para[2])
  init[1,3] =  lmom_est$para[3]

  if( abs(lmom_est$para[3]) > 0.5) {
    init[1,3] = sign(init[1,3])*0.48 }

  maxm1=ntry-2; maxm2=ntry-3
  init[2:maxm1,1] <- init[1,1]+rnorm(n=maxm2,mean=0,sd = 20)
  init[2:maxm1,2] <- log(lmom_est$para[2])+rnorm(n=maxm2,mean=0,sd = 1)
  init[2:maxm1,3] <- runif(n=maxm2,min= -0.49, max=0.49)

  mx = mean(data)
  sx= log(sqrt(var(data)))
  init[ntry-1,1:3] = c(mx, sx, pretheta[5])
  init[ntry,1:3] =   c(pretheta[1], pretheta[3], pretheta[5]+.05)
  return(init)
}
#----------------------------------------------------
qns.gev_all= function(T=NULL, para, year, model=NULL){

  nsample=length(year)
  ns=nsample
  zpT=rep(NA, nsample)
  year2=year^2

  sp=set.para.model(para,model)

  xi= sp$xi
  zpc= (1- ( -log(1-(1/T) ) )^xi ) /xi
  zpT[1:ns] = sp$mu[1] + sp$mu[2]*year[1:ns] +sp$mu[3]*year2[1:ns]
  zpT[1:ns] = zpT[1:ns]+ zpc* exp( sp$sig[1] +sp$sig[2]*year[1:ns] )

  return(zpT)
}

#--------------------------------------------
set.para.model = function(para, model=NULL){

  z=list()
  mu0=para[1]
  mu1=para[2]

  if(model=='gev10'){
    sig0=para[3]
    xi=para[4]
    mu2=0
    sig1=0
  }else if(model=='gev20'){
    mu2=para[3]
    sig0=para[4]
    xi=para[5]
    sig1=0
  }else if(model=='gev11'){
    sig0=para[3]
    sig1=para[4]
    xi=para[5]
    mu2=0
  }

  z$mu=c(mu0,mu1,mu2)
  z$sig=c(sig0,sig1)
  z$xi=xi
  return(z)
}
#---------------------------------------------------------
