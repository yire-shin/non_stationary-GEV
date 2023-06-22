library(lmomco)
library(nleqslv)
library(robustbase)

# R code for NS GEV11 model by Yire Shin, Sanghoo Yoon, and JS Park (2023)
#----------------------------------------------------------------------
gado110.m0s0= function (xdat, ntry=20, ftol=1e-6){
  
  z <- list()
  ns=length(xdat)
  year=seq(1,length(xdat))
  para.gado=rep(NA,5)
  
  name_gev11_ns  =c("mu0","mu1","sigma0","sigma1","xi")
  name_gev00_sta =c("mu","sigma","xi")
  
  reg.dat=data.frame( cbind(year, xdat) )
  
  mu.init= lmrob(xdat~year, reg.dat)$coefficients  # robust regression
  m0= mu.init[1]
  m1= mu.init[2]
  
  orig.para=c(m0,m1,1.0,-0.001,0)
  
  strup.tF = strup(xdat, orig.para=orig.para, mtrim=F, 
                   rewt=T)                            #WLS by strup
  
                               # rewt=T: iterative re-weighted LSE
                               # mtrim=T: trimmed L-ME
  
  qlist= make.qmax(xdat, orig.para=orig.para) 
  
  orig.para=c(m0,m1, qlist$sig0, qlist$sig1, 0)
  
  gado.rob = time.m(qmax=qlist$qmax, orig.para=orig.para, 
                    mtrim=F)                          # method of GN16
  
  z = multi.m0s0(xdat, ntry=ntry, ftol=ftol,
                 pretheta=gado.rob$para.org, mtrim=F)  #proposed method
  
                                            # z$para.prop = proposed est
  
  if(z$precis > ftol) z$para.prop = gado.rob$para.up  
  
  z$para.gado.org  =gado.rob$para.org        # GN16 orginal est
  z$para.gado.up   =gado.rob$para.up         # mdfd GN16 est
  #  z$gado.mu_t = gado.rob$mu_t
  #  z$gado.mu_t.up = gado.rob$mu_t.up
  
  z$strup.para.sta =strup.tF$strup.sta      # stationary wlse
  z$strup.org      =strup.tF$strup.para     # wlse by strup
  z$strup.final    =strup.tF$wls.par        # specified WLSE
  
  z$lme.sta= pargev(lmoms(xdat,nmom=5))$para   # stationary L-ME
  
  names(z$para.prop)     <-name_gev11_ns
  names(z$para.gado.org) <-name_gev11_ns
  names(z$para.gado.up)  <-name_gev11_ns
  names(z$strup.org)     <-name_gev11_ns
  names(z$strup.final)   <-name_gev11_ns
  
  names(z$strup.para.sta) <-name_gev00_sta
  names(z$lme.sta)        <-name_gev00_sta
  
  return(z)
}

#-----strup wls -----------------------------------
strup =function(xdat, orig.para=orig.para, mtrim=F, rewt=T){
  
  w=list()
  ns=length(xdat)
  year=seq(1,ns)
  
  ftol=1e-4
  iter=0
  dist=1000
  
  m0=orig.para[1]
  m1=orig.para[2]
  orig.para[3:4]= 0
  res=xdat -(m0+m1*year)
  new.para=orig.para
  
  while( dist > ftol) { 
    iter=iter+1
    para=new.para  
    
    stand=wls.park(res)
    
    new.para=para[1:4]+c(stand$m, stand$sig)
    res=stand$res
    dist= sum(abs(new.para[1:4]-para[1:4]))
    if(iter > 20) dist=0
    if(rewt==F) dist=0      # rewt=T: iterative re-weighted LSE
  }
  
  ares= xdat-(new.para[1]+new.para[2]*year)
  ares= ares/exp(new.para[3]+new.para[4]*year)   # standardization
  
  w$strup.sta = pargev(lmoms(ares))$para
  if(mtrim==T){                            # mtrim=T: trimmed L-ME
    w$strup.sta = TLpargev(ares,leftrim=3)
  }
  
  # modified wls by park: standardization first ----------
  
  res.std = (ares -w$strup.sta[1])/w$strup.sta[2]
  w$wls.sta = pargev(lmoms(res.std))$para       
  
  res.std2 = (res.std -w$wls.sta[1])/w$wls.sta[2]
  w$wls.sta2 = pargev(lmoms(res.std2))$para      # standardization-2
  
  w$strup.para = c(new.para[1:4], w$wls.sta2[3])
  
  #------ optim to specify parameter values for return period T=Trp
  
  Trp=50
  init4=matrix(0,20,5)
  init4[1,]= w$strup.para[1:5]
  init4[2,]= c(w$strup.sta[1],0,w$strup.sta[2],0,w$strup.sta[3])
  init4[3,]= (init4[1,]+init4[2,])/2
  init4[4,]= init4[3,]+ 0.1
  mx = mean(xdat)
  sx= log(sqrt(var(xdat)))
  init4[5,] = c(mx,0.2, sx, 1.0, -0.3)
  init4[6,] = c(mx,0.1, 2.0, -0.1, -0.1)
  init4[7,] = c(w$strup.para[1:2],2.0, 0.03, 0.4)
  init4[8,] = c(mx,0.1, sx, -0.1, w$wls.sta2[3])
  init4[9,] = c(mx,-0.1, sx, 0.1, -0.4)
  init4[10,] = c(mx,0.01,w$strup.para[3:5])
  
  ming=c(w$strup.sta, w$strup.para, w$wls.sta,  
         new.para[4], w$wls.sta2[3] ) 
  wls=list()
  sq=rep(NA,20)
  
  # ------ number of trials -----------------  
  for (i in 1:10){
    init=as.vector(init4[i,1:3])      # optim for first 3 paras
    wls[[i]] = tryCatch( optim(par=init, fn=wls.para, method=c("BFGS"),
                               ming=ming, year=year,T=Trp) )
    sq[i] = wls[[i]]$value
  }
  
  w$wls.par= c(wls[[which.min(sq)]]$par, new.para[4], w$wls.sta2[3])
  # sigma_1 and xi are fixed
  return(w)
}
#---------------------------------------------------------- 
#----------------------------------------------------------

wls.para=function(init4, ming=ming, year=year,T=100){
  
  strup.para.sta=ming[1:3]
  strup.para=ming[4:8]
  wls.sta =ming[9:11]
  
  paraT=c(as.vector(init4), ming[12:13])
  Trp=T
  
  zpT.stF= quagev(1-(1/Trp), vec2par(strup.para.sta, 'gev'))
  zpT.lme= zpT.stF*exp(strup.para[3]+strup.para[4]*year)
  zpT.lme= zpT.lme+ strup.para[1] +strup.para[2]*year
  zpT.lme= zpT.lme*wls.sta[2]+wls.sta[1]
  
  qt.wls = qns.gev110(T=Trp, paraT, year)
  
  q= sum((qt.wls-zpT.lme)^2)
  
  return(q)
}
#------------------------------------------------------
wls.park = function(res ){
  
  z=list()
  ns=length(res)
  year=seq(1,ns)
  # m0=orig.para[1]
  # m1=orig.para[2]
  # res = xdat -(m0+m1*year)
  
  lres.pr=log(abs(res))
  sig.dat=data.frame( cbind(year, lres.pr) )
  z$sig= lm(lres.pr~year, sig.dat)$coefficients
  
  sigt = exp(z$sig[1] + z$sig[2]* year)
  
  res.n= res/sigt
  new.data= data.frame( cbind(year, res.n) )
  z$m= lmrob(res.n~year, new.data)$coefficients
  z$res = res.n -(z$m[1]+z$m[2]*year) 
  
  return(z)
}
#------------------------------------------------------
#---------------------------------------------------
TLpargev=function(xdat, leftrim=5){
  
  value=list()
  ftol=1e-5
  
  q.sta = pargev(lmoms(xdat))$para
  lmom_q = TLmoms(xdat, leftrim=leftrim)
  
  value =  nleqslv( x=q.sta, 
                    fn= run.TLgev, lmom=lmom_q,
                    leftrim=leftrim) 
  
  precis=  mean(abs(value$fvec))
  if( precis < ftol) { 
    return(value$x)
  }else{
    stop("no result at TLpargev","\n")
  }
}
#---------------------------------------------------    
run.TLgev = function(x, lmom=NULL, leftrim=5){
  
  z=rep(NA,3)  
  theo=theoTLmoms(vec2par(x,'gev'), leftrim=leftrim)
  
  z[1]=lmom$lambdas[1]-theo$lambdas[1]
  z[2]=lmom$lambdas[2]-theo$lambdas[2]
  z[3]=lmom$ratios[3]-theo$ratios[3]
  return(z)
}
#---------------------------------------------------------------  
time.m = function(qmax=NULL, orig.para=NULL, mtrim=T){  
  
  z=list()
  para.gado=rep(NA,5)
  ns=length(qmax)
  year=seq(1,ns)
  
  m0=orig.para[1]
  m1=orig.para[2]
  sig0=orig.para[3]
  sig1=orig.para[4]
  
  if(mtrim==F){
    lmom_q = lmoms(qmax)
    q.sta = pargev(lmom_q)$para
    
  }else if(mtrim==T){
    #    leftrim= round(ns*0.1)
    
    leftrim=5
    q.sta = TLpargev(qmax, leftrim=leftrim)
  }
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
  #--------------------------------------------------------
  
  mu.data= data.frame( cbind(year, mu.gado) )
  loc.gado =lm(mu.gado~year, mu.data)$coefficients
  
  alpha0 = log(cd) + sig0  #+0.5
  alpha1 = sig1
  
  alpha0_up = log(q.sta[2])- sig1*nh
  
  z$para.org= c(loc.gado, alpha0, alpha1, xi )
  z$para.up = c(loc.gado, alpha0_up, alpha1, xi )
  
  alpha_t = exp(alpha0_up +sig1*year)
  z$mu_t.up=  - (1-gamma(1+ xi) )*alpha_t/xi  + m0+ m1*year
  
  return(z)
}

#-------------------------------------------------  
gev.Ldist.m0s0 <- function(a, xdat=xdat, pretheta=pretheta,
                           mtrim=F) 
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
  if(mtrim==F){
    lam= lmomgum(vec2par(c(0,1),'gum'))
  }else if(mtrim==T){
    lam= theoTLmoms(vec2par(c(0,1),'gum'), leftrim=5)
  }
  
  lgum=list()
  lgum=lmoms.md.park(newg, mtrim=mtrim, no.stop=T)
  
  if(lgum$ifail == 1) {
    zz[1:3]=1000
    #      cat("Only one value in newg","\n")
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
multi.m0s0= function(xdat, ntry=20, ftol=1e-4, 
                     pretheta=pretheta, mtrim=F)
{
  z=list()
  value=list()
  k=list()
  
  init= matrix(0, nrow=ntry, ncol=3)
  init = ginit.m0s0(xdat, ntry, pretheta)
  
  precis=rep(1000, ntry)
  para.sel=matrix(NA,ntry+1,ncol=5)
  
  tryCatch({
    for(i in 1:ntry) {
      
      value =  tryCatch( nleqslv( x=as.vector(init[i,1:3]), 
                                  fn= gev.Ldist.m0s0,
                                  method="Broyden",
                                  xdat=xdat, pretheta=pretheta,
                                  mtrim=mtrim) )
      
      k[[i]] <- value
      
      if(is(value)[1]=="try-error"){
        k[[i]]$fvec <- 10^6
        k[[i]]$termcd = 5
      }else{
        
        precis[i]=  mean(abs(k[[i]]$fvec) )
        
        if( precis[i] < ftol) {
          
          k[[i]]$root = value$x
          
          # if(abs(k[[i]]$root[3]) >0.5) {
          #   k[[i]]$root[3] =sign(k[[i]]$root[3])*0.4999
          # }
          
          para.sel[i,1:5]=c( k[[i]]$root[1], pretheta[2],
                             k[[i]]$root[2], pretheta[4],
                             k[[i]]$root[3])
          
          #            ks[i]= KS_gev110(xdat, para)$KS_ST[1]
        }
      }
      
      precis[is.na(precis[i])]=1000
      if( abs( k[[i]]$termcd ) > 3 ) {
        precis[i]=1000
        # k[[i]]$root = c(pretheta[1], pretheta[3],
        #                 pretheta[5])
        para.sel[i,]=NA
      }
    } #end for
  }) # trycathch
  
  z$para.prop =sel.para(xdat,para.sel)$para
  
  #   ns=length(xdat)
  #   lmom_init = lmoms(xdat,nmom=5)
  #   lme.sta <- pargev(lmom_init)$para
  # 
  #   mu0=  lme.sta[1]-pretheta[2]*ns
  #   sig0= log(lme.sta[2])-pretheta[4]*ns
  # 
  #   para.sel[ntry+1,1:5] =c(mu0, pretheta[2], sig0,
  #                           pretheta[4], lme.sta[3])
  # 
  # z$para.incy =sel.para(xdat,para.sel)
  z$precis =precis[which.min(precis)]
  
  return(z)
}
#-------------------------------------------------
sel.para =function(xdat, para.sel){
  
  z=list()
  upara.sel = para.sel
  gof=rep(NA,nrow(upara.sel))
  
  ns=length(xdat)
  #  Tmax =min(40,ns)
  vecT=c(5,10,20,40)
  if(ns >= 100) vecT=c(5,10,20,40,80)
  
  for(i in 1:nrow(upara.sel) ){
    gof[i] = gof.ene(xdat, vecT, as.vector(upara.sel[i,1:5]))
    #sum(abs(upara.sel[i,]))
  }
  
  if(length(unique(gof))==1) {
    z$para =upara.sel[length(gof),]
  }else{
    z$para =upara.sel[which.min(gof),]
  }
  z$gof=gof
  return(z)
}
#-----------------------------------------
#-------------------------------------------------
gof.ene = function(xdat, vecT=c(5,10,20,40,80,120), para){
  
  ns=length(xdat)
  nT = length(vecT)
  year=seq(1,ns)
  chi=rep(NA,nT)
  
  for(i in 1:nT){
    T=vecT[i]
    qt=qns.gev110(T, para, year)
    ene = ns/T
    sne = sum(xdat >= qt)
    #    chi[i] = ((ene-sne)^2 )/ene
    
    chi[i] = ( abs(ene-sne) )/ene
  }
  chi2 =sum(chi)
  return(chi2)
}

#------------------------------------------------  
eq =function(q){
  eq=quagum(q,vec2par(c(0,1),'gum') )*20 +300
  return(eq)
}
sq = function(q, newg=newg){
  sq=quantile(newg, q)
  return(sq)
}

#----------------------------------------------------
qns.gev110= function(T=50, para, year){
  
  nsample=length(year)
  ns=nsample
  zpT=rep(NA, nsample)
  
  mu0=para[1]
  mu1=para[2]
  sigma0=para[3]
  sigma1=para[4]
  xi=para[5]
  
  zpc= (1- ( -log(1-(1/T) ) )^xi ) /xi
  zpT[1:ns]= mu0+mu1*year[1:ns] + zpc*exp(sigma0 + sigma1*year[1:ns])
  
  return(zpT)
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
make.qmax =function(xdat=NULL, orig.para=NULL)
  # iter=F, new.sig=NULL)
{
  
  z=list()
  m0= orig.para[1]                     
  m1= orig.para[2]     
  ns = length(xdat)
  year= seq(1,ns)
  
  res = xdat - (m0 + m1* year)
  mres= mean(res)
  res.pr = abs(res - mres)
  
  #  if(iter==F){
  lres.pr=log(res.pr)
  sig.dat=data.frame( cbind(year, lres.pr) )
  sig.lm= lm(lres.pr~year, sig.dat)$coefficients
  sig0 = sig.lm[1]
  sig1 = sig.lm[2]
  # }else if(iter==T){
  #   sig0= new.sig[3]
  #   sig1= new.sig[4]
  # }
  
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
  
  #     qmax2 = res/sigt

  z$qmax=qmax
  #  z$qmax2=qmax2
  z$sig0=sig0
  z$sig1=sig1
  z$sigt=sigt
  
  return(z)
}
#----------------------------------------------
#------------------------------------------------------
ginit.m0s0 <-function(data,ntry=ntry, pretheta){
  
  init <-matrix(0, nrow=ntry, ncol=3)
  if(abs(pretheta[5]) > 0.5) pretheta[5] = sign(pretheta[5])*0.48
  
  lmom_init = lmoms(data,nmom=5)
  lmom_est <- pargev(lmom_init)
  
  init[1,1] <- lmom_est$para[1]
  init[1,2] <- log(lmom_est$para[2])
  init[1,3] = lmom_est$para[3]
  
  if( abs(lmom_est$para[3]) > 0.5) {
    init[1,3] = sign(init[1,3])*0.48
  }
  
  maxm1=ntry-2; maxm2=ntry-3
  init[2:maxm1,1] <- init[1,1]+rnorm(n=maxm2,mean=0,sd = 20)
  init[2:maxm1,2] <- log(lmom_est$para[2])+rnorm(n=maxm2,mean=0,sd = 1)
  init[2:maxm1,3] <- runif(n=maxm2,min= -0.49, max=0.49)
  
  mx = mean(xdat)
  sx= log(sqrt(var(xdat)))
  init[ntry-1,1:3] = c(mx, sx, pretheta[5])
  init[ntry,1:3]= c(pretheta[1], pretheta[3], pretheta[5]+.05)
  return(init)
}
# #--------------------------------------------------------



