# lemma.R
# 2010-03-19
# version 1.3
# Haim Bar, Elizabeth Schifano {hyb2,eds27}@cornell.edu
#
# Allow for two nonnull groups
# Allow topgenes to be all, or nonnull (to print all genes, or
#   all the nonnulls, repectively)
# Implmented a function to print the top genes separately from
#   the estimation phase
# Incorporated the plotting function into the main R file
# 
# ver. 1.3:
# Added an option to get Method-of-Moments estimates for the
#   hyperparameters of prior distribution of the error variance

if (getRversion() >= "2.15.1") globalVariables(c("G","dg","n2","n1","p0","tau","sig2eb","p1","psi","sig2psi","p2","mg","f","alpha_hat","beta_hat","pBH0","RRfdr0","RRfdr1","RRfdr2","genename"))

##############################
#
# The EM algorithm functions - estimating the parameters of the
# mixture model.  We are assuming that the variable dg is a mixture of
# 3 Normal distributions (one null and two nonnull groups):
# dg ~ p0*N(tau,s2) + p1*N(tau+psi,s2+sig2psi) + p2*N(tau-psi,s2+sig2psi)
# where p0+p1+p2 = 1
# If the user wants to fit a mixture of two normal distributions
# by setting modes=2, we set p2=0
# 
get.EM.est3 <- function(dg, sig2, n1, n2, logf, tol=0.0001, maxIts=20000,
     modes=3){
  out = vector("list", 5)
  names(out) = c("sig2psi", "psi", "p1", "p2", "tau")

  # set initial values:
  tau.old     = 0
  psi.old     = 20
  p1.old      = 0.1
  sig2psi.old = 1

  if (modes == 3) {
    p2.old = 0.1
  }
  else { p2.old = 0 }

  if (n2 == 1) {   # within-group (paired) analysis
    sd0 = sqrt((1/n1)*sig2)
    sd1.old = sqrt((1/n1)*sig2 + sig2psi.old)
  }
  else {           # two-group comparison:
    sd0 = sqrt((1/n1+1/n2)*sig2)
    sd1.old = sqrt((1/n1+1/n2)*sig2 + sig2psi.old)
  }

  # w1 correspond to the indicator variables that are 1 if gene g is non-null #1
  w1.old = posteriorProbability("nonnull1",(1-p1.old-p2.old),p1.old,p2.old,
     dg,tau.old,psi.old,sd0,sd1.old)

  # w2 correspond to the indicator variables that are 1 if gene g is non-null #2
  if (modes == 3) {
    w2.old = posteriorProbability("nonnull2",(1-p1.old-p2.old),p1.old,p2.old,
          dg,tau.old,psi.old,sd0,sd1.old)
  }
  else { w2.old = 0 }
  aerr  = 1 # aerr = the cumulative squared difference between the estimates in
            # the k-th and (k+1)-th iterations
  count = 0
  while (aerr>tol) {  # check convergence of the EM algorithm
    out.sig2psi = try(uniroot(f=Mstep.sig2psi, lower = -200,
        upper = 200, X=dg-tau.old,sig2=sig2, w1=w1.old,
        w2=w2.old,n1=n1, n2= n2, psi=(psi.old)),silent=T)
# load("projects/LEMMA-paper/simulations/SIM/LEMMA/5.00/0.01/1.00/4.00/S04.RData")

    if (!is.character(out.sig2psi))
    { sig2psi.cur = exp(out.sig2psi$root) }
    else {
      sig2psi.cur = 0
      p1.old=0
      p2.old=0
    }
    
    psi.cur = Mstep.psi(dg-tau.old, sig2, w1.old, w2.old,  n1, n2,
          sig2psi.cur)

    if (n2 == 1) {
      v0 = (1-w1.old-w2.old)/((1/n1)*sig2)
      v1 = w1.old/(((1/n1)*sig2)+sig2psi.cur)
      v2 = w2.old/(((1/n1)*sig2)+sig2psi.cur)
    }
    else {
      v0 = (1-w1.old-w2.old)/((1/n1+1/n2)*sig2)
      v1 = w1.old/(((1/n1+1/n2)*sig2)+sig2psi.cur)
      v2 = w2.old/(((1/n1+1/n2)*sig2)+sig2psi.cur)
    }
    if(sum(v1+v2) == 0) {
      tau.cur = sum(dg/((1/n1+1/n2)*sig2))/sum(1/((1/n1+1/n2)*sig2))
    } else {
      u1 = sum((v1-v2)*dg)/sum(v1+v2)
      u2 = sum(v1-v2)/sum(v1+v2)
      tau.cur = ( ( sum((v0+v1+v2)*dg) - (sum(v1-v2))*u1 ) /
                  ( sum(v0+v1+v2) - (sum(v1-v2))*u2 ))
    }

    if (n2 == 1) {
      sd1.cur= sqrt(sig2psi.cur + (1/n1)*sig2)
    }
    else {
      sd1.cur= sqrt(sig2psi.cur + (1/n1+1/n2)*sig2)
    }

    w1.cur = posteriorProbability("nonnull1",(1-p1.old-p2.old),p1.old,p2.old,
          dg,tau.cur,psi.cur,sd0,sd1.cur)
    if (modes == 3) {
      w2.cur = posteriorProbability("nonnull2",(1-p1.old-p2.old),p1.old,p2.old,
          dg,tau.cur,psi.cur,sd0,sd1.cur)
    }
    else {
      w2.cur = 0
    }
    p1.cur  = mean(w1.cur)
    p2.cur  = mean(w2.cur)

    aerr = sqrt((p1.cur - p1.old)^2 + (p2.cur - p2.old)^2 + 
               (sig2psi.cur-sig2psi.old)^2 +
               (psi.cur-psi.old)^2 + (tau.cur-tau.old)^2 )
    
    sd1.old = sd1.cur
    sig2psi.old = sig2psi.cur
    w1.old = w1.cur  
    w2.old = w2.cur  
    p1.old = p1.cur
    p2.old = p2.cur
    psi.old = psi.cur
    tau.old = tau.cur
    count = count+1
    if (count > maxIts) {
      cat("WARNING! Convergence not achieved after ",maxIts,
         " iterations (norm_2(phi.new - phi.old) = ", aerr," > ", tol,")\n",
        file=logf);
      break
    }
  }

  cat("  Number of EM iterations: ",count,"\n",file=logf)
  out$sig2psi = sig2psi.cur
  out$psi = psi.cur
  out$p1 = p1.cur
  out$p2 = p2.cur
  out$tau = tau.cur
  return(out)
}

Mstep.sig2psi <- function(lsig2psi, X, sig2, w1,w2, n1, n2, psi){
  sig2psi = exp(lsig2psi)
  if (sum(w1+w2) == 0) {
    0
  }
  else {
    if (n2 == 1) {
      V = (1/n1)*sig2 + sig2psi
    }
    else {
      V= (1/n1+1/n2)*sig2 + sig2psi
    }
    sum( (w1+w2)/V - (w1*(psi-X)^2 + w2*(psi+X)^2)/(V^2) )
  }
}

Mstep.psi <- function(X, sig2, w1, w2, n1, n2, sig2psi){
  if (sum(w1+w2) == 0) {
    0
  } else {
    if (n2 == 1) {
      V = (1/n1)*sig2 + sig2psi
    }
    else {
      V = (1/n1+1/n2)*sig2 + sig2psi
    }
    sum(((w1-w2)/V)*X) / sum((w1+w2)/V)
  }
}

mg.marginal.nloglik <- function(pars, MSE, f){
  alpha = pars[1]
  beta = pars[2]
  -sum(log ((gamma(alpha+f/2)*(MSE^(f/2-1)) /
         ((1/beta+MSE*f/2)^(f/2+alpha)*(beta^alpha)*gamma(alpha))) ))
}


###################
#
# Computing confidence intervals for the estimates
#
varAlphaBeta <- function(alpha,beta,f,mg) {
  M = matrix(0,nrow=2,ncol=2)
  M[1,1] = length(mg)*(trigamma(f/2+alpha)-trigamma(alpha))
  M[1,2] = sum(-(mg*f)/(beta*mg*f+2))
  M[2,1] = M[1,2]
  M[2,2] = sum((-2*beta*mg*f^2-2*f+alpha*(beta*mg*f)^2)/(mg*f*beta^2+2*beta)^2)
  -solve(M)
}

varPhi <- function(sig2eb,n1,n2,sig2psi,p0,dg,tau,p1,psi,p2) {
  if (n2 == 1) {   # within-group (paired) analysis
    v0 = (1/n1)*sig2eb
  } else {
    v0 = sig2eb*(1/n1+1/n2)
  }
  sd0 = sqrt(v0)
  v1 = v0 + sig2psi
  sd1 = sqrt(v1)
  w0 = posteriorProbability("null",(1-p1-p2),p1,p2, dg,tau,psi,sd0,sd1)
  w1 = posteriorProbability("nonnull1",(1-p1-p2),p1,p2, dg,tau,psi,sd0,sd1)
  w2 = posteriorProbability("nonnull2",(1-p1-p2),p1,p2, dg,tau,psi,sd0,sd1)

  # using Oakes' method:
  B = matrix(0,nrow=5, ncol=5)
  B[1,1] = -sum((w1)/(p1^2) + w0/(p0^2))
  B[1,2] = -sum(w0/(p0^2))
  B[2,1] = B[1,2]
  B[2,2] = -sum((w2)/(p2^2) + w0/(p0^2))
  B[3,3] = -sum(w0/v0 + (w1+w2)/v1)
  B[3,4] = -sum((w1-w2)/v1)
  B[3,5] = -sum((w1*(dg-tau-psi)+w2*(dg-tau+psi))/(v1^2))
  B[4,3] = B[3,4]
  B[5,3] = B[3,5]
  B[4,4] = -sum((w1+w2)/v1)
  B[4,5] = -sum((w1*(dg-tau-psi) - w2*(dg-tau+psi))/(v1^2))
  B[5,4] = B[4,5]
  B[5,5] = -sum( (w1*(dg-tau-psi)^2 + w2*(dg-tau+psi)^2)/(v1^3) -
                 (w1+w2)/(2*v1^2) )
  # B is the matrix of second derivatives

  f0 = dnorm(dg-tau,  0,   sd0)
  f1 = dnorm(dg-tau,  psi, sd1)
  f2 = dnorm(dg-tau, -psi, sd1)
  M = p1*f1 + p2*f2 + p0*f0
  Mp1 = f1 - f0
  Mp2 = f2 - f0
  Mtau = p0*f0*(dg-tau)/v0 + p1*f1*(dg-tau-psi)/v1 + p2*f2*(dg-tau+psi)/v1
  Mpsi = p1*f1*(dg-tau-psi)/v1 - p2*f2*(dg-tau+psi)/v1
  Msig2psi = (p1*f1*((dg-tau-psi)^2-v1)/(2*v1^2) +
              p2*f2*((dg-tau+psi)^2-v1)/(2*v1^2) )
  
  C = matrix(0,nrow=5, ncol=5)
  C[1,1] = sum( ((f1*M-p1*f1*Mp1)/p1 + ((f0*M + p0*f0*Mp1)/p0))/(M^2))
  C[1,2] = sum((-f1*Mp2 + (f0*M + p0*f0*Mp2)/p0)/(M^2))
  C[2,1] = C[1,2]
  C[1,3] = sum((f1*((dg-tau-psi)*M/v1 -Mtau) - f0*((dg-tau)*M/v0 -Mtau))/(M^2))
  C[3,1] = C[1,3]
  C[1,4] = sum((f1*((dg-tau-psi)*M/v1 - Mpsi) + f0*Mpsi)/(M^2))
  C[4,1] = C[1,4]
  C[1,5] = sum((f1*( ((dg-tau-psi)^2 - v1)*M/(2*v1^2) -Msig2psi)
              + f0*Msig2psi)/(M^2))
  C[5,1] = C[1,5]
  
  C[2,2] = sum(((f2*M-p2*f2*Mp2)/p2 + (f0*M+p0*f0*Mp2)/p0)/(M^2))
  C[2,3] = sum((f2*((dg-tau+psi)*M/v1-Mtau) - f0*((dg-tau)*M/v0-Mtau))/(M^2))
  C[3,2] = C[2,3]
  C[2,4] = sum( (f2*(-(dg-tau+psi)*M/v1-Mpsi) + f0*Mpsi)/(M^2))
  C[4,2] = C[2,4]
  C[2,5] = sum( (f2*(((dg-tau+psi)^2-v1)*M/(2*v1^2)-Msig2psi)
                + f0*Msig2psi)/(M^2))
  C[5,2] = C[2,5]
  
  C[3,3] = sum( (p0*f0*((dg-tau)*M/v0-Mtau)*(dg-tau)/v0
           + (p1*f1*((dg-tau-psi)*M/v1-Mtau)*(dg-tau-psi) 
           +  p2*f2*((dg-tau+psi)*M/v1-Mtau)*(dg-tau+psi)) /v1)  /(M^2) )
  C[3,4] = sum(( -p0*f0*Mpsi*(dg-tau)/v0
           + (p1*f1*((dg-tau-psi)*M/v1-Mpsi)*(dg-tau-psi) 
           -  p2*f2*((dg-tau+psi)*M/v1+Mpsi)*(dg-tau+psi))/v1 )/(M^2))
  C[4,3] = C[3,4]
  C[3,5] = sum(( -p0*f0*Msig2psi*(dg-tau)/v0 
            + (p1*f1*( ((dg-tau-psi)^2-v1)*M/(2*v1^2)-Msig2psi )*(dg-tau-psi)
            +  p2*f2*( ((dg-tau+psi)^2-v1)*M/(2*v1^2)-Msig2psi )*(dg-tau+psi))
              /v1 )/(M^2))
  C[5,3] = C[3,5]
  
  C[4,4] = sum(((p1*f1*((dg-tau-psi)*M/v1-Mpsi)*(dg-tau-psi)
              + p2*f2*((dg-tau+psi)*M/v1+Mpsi)*(dg-tau+psi))/v1)/(M^2))
  C[4,5] = sum(((p1*f1*( ((dg-tau-psi)^2-v1)*M/(2*v1^2)-Msig2psi )*(dg-tau-psi)
          - p2*f2*( ((dg-tau+psi)^2-v1)*M/(2*v1^2)-Msig2psi )*(dg-tau+psi))/v1
           )/(M^2))
  C[5,4] = C[4,5]
  
  C[5,5] = sum( (
      ( p1*f1*( ((dg-tau-psi)^2-v1)*M/(2*v1^2)-Msig2psi ) *(dg-tau-psi)^2
       +p2*f2*( ((dg-tau+psi)^2-v1)*M/(2*v1^2)-Msig2psi ) *(dg-tau+psi)^2 )/
      (2*v1^2)
     -( p1*f1*( ((dg-tau-psi)^2-v1)*M/(2*v1^2)-Msig2psi )
       +p2*f2*( ((dg-tau+psi)^2-v1)*M/(2*v1^2)-Msig2psi ) )/
      (2*v1) ) /
     (M^2) )
  INM = -(B+C)
}

posteriorProbability <- function(group,p0,p1,p2,y,tau,psi,sd0,sdA) {
  m = p0*dnorm(y-tau, 0, sd0) + p1*dnorm(y-tau, psi, sdA) +
        p2*dnorm(y-tau, -psi, sdA)
  switch(group,
     null = p0*dnorm(y-tau, 0, sd0)/m,
     nonnull1 = p1*dnorm(y-tau, psi, sdA)/m,
     nonnull2 = p2*dnorm(y-tau, -psi, sdA)/m)
}

###################
#
# Plotting functions:
#
#   note: if you run these functions separately from the estimation,
#         load the AllData.RData file, first

###################
# Fitting the M.S.E.'s
#
posteriorIG <- function(m, f, alpha, beta) {
  res = lgamma(f/2 + alpha)+(f/2-1)*log(m)+(f/2)*log(f/2) -
        (f/2+alpha)*log(1/beta+m*f/2) - lgamma(f/2)- lgamma(alpha) -
        alpha*log(beta)
  return(exp(res))
}

lemmaPlots <- function(outdir, mgq=0.99, titletext, modes=3) {
  # load the data:
  lemmaDataFile = sprintf("%s/AllData.RData",outdir)
  load(lemmaDataFile)

###################
# Plot dg and the fitted mixture:
#
  plotname1 = sprintf("%s/LEMMAplots1.ps", outdir)
  postscript(file=plotname1, width=10, height=7.5, paper="letter",horizontal=T,
   family="Times")
     ## allocate figure 1 all of row 1
     ## allocate figure 2 the intersection of column 2 and row 2
  layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4), 2, 6, byrow = TRUE))
  par(oma=c(0,0,4,0),mar=c(5,4.1,4.1,2.1))

#  par(oma=c(0,0,4,0),mfrow=c(1,2),mar=c(5,4.1,4.1,2.1))
  brks = min(100,max(floor(G*0.1),20))
  x = seq(min(dg),max(dg),by=((max(dg)-min(dg))/100))
  if (n2 == 1) {
    dn = 1/n1
  }
  else {
    dn = 1/n1 + 1/n2
  }
  ynull = p0*dnorm(sort(dg),tau,sqrt((dn*sig2eb[order(dg)])))
  yalt1 = p1*dnorm(sort(dg),tau+psi,
            sqrt(sig2psi+(dn*sig2eb[order(dg)])))
  yalt2 = p2*dnorm(sort(dg),tau-psi,
            sqrt(sig2psi+(dn*sig2eb[order(dg)])))
  ymix = (p1*dnorm(x,tau+psi,sqrt(sig2psi + mean((dn*sig2eb))))
           +p2*dnorm(x,tau-psi,sqrt(sig2psi + mean((dn*sig2eb))))
           +p0*dnorm(x,tau,sqrt(mean((dn*sig2eb)))))

  ymax = max(ynull,yalt1,yalt2,ymix,max(hist(dg,breaks=40,plot=F)$density))
  ylim = c(0,ymax*1.2)
  maxsig = max(sig2psi,sig2psi)+ mean(dn*sig2eb)
  xmax = tau+abs(psi)+2.5*sqrt(maxsig )
  xmin = tau-abs(psi)-2.5*sqrt(maxsig )
  plot(sort(dg), ynull,type='l',col='lightskyblue1', lwd=1,
    xlab=expression(d[g]),
    main=bquote(bold(paste("(a) mean(",d[g],")=",.(round(mean(dg),2)),
    "  s.d.(",d[g],")=",.(round(sd(dg),2))))), ylab="Density", ylim=ylim,
    xlim=c(xmin,xmax), axes=F)
  axis(1,floor(min(dg)):ceiling(max(dg)))
  axis(2)
  lines(sort(dg), yalt1, col='pink3', lwd=1)
  lines(sort(dg), yalt2, col='yellow', lwd=1)

  hist(dg,breaks=brks,probability=TRUE,border="gray66",add=T)
  lines(x,p0*dnorm(x,tau,sqrt(mean((dn*sig2eb)))),col="navyblue",
     lwd=3,lty=2)
  lines(x,p1*dnorm(x,tau+psi,sqrt(sig2psi + mean((dn*sig2eb)))),
     col="red", lwd=2, lty=4)
  if (modes == 3) {
    lines(x,p2*dnorm(x,tau-psi,sqrt(sig2psi + mean((dn*sig2eb)))),
       col="orange", lwd=2, lty=4)
  }
  lines(x,ymix,col="black", lwd=1,lty=1)

  dig = 4
  if (abs(min(p1,tau,psi,sig2psi)) >= 0.01) {
    dig = 2
  } else if (abs(min(p1,tau,psi,sig2psi)) >= 0.001) {
    dig = 3
  }

  mixtext = bquote(paste(psi[g]," ~ N(",psi,",",sigma[psi]^2,")"))
  text(xmin,0.95*ylim[2],mixtext, pos=4, cex=0.8)
  text(xmin,0.85*ylim[2],bquote(paste(hat(psi)," = ",
       .(round(psi,dig)))), pos=4, cex=0.8)
  text(xmin,0.75*ylim[2],bquote(paste(hat(sigma)[psi]^2," = ",
       .(round(sig2psi,dig)))), pos=4, cex=0.8)
  text(xmin,0.65*ylim[2],bquote(paste(hat(tau)," = ",
       .(round(tau,dig)))), pos=4, cex=0.8)
  text(xmin,0.55*ylim[2],bquote(paste(hat(p)[1]," = ",
       .(round(p1,dig)))), pos=4, cex=0.8)
  if (modes == 3) {
    text(xmin,0.45*ylim[2],bquote(paste(hat(p)[2]," = ",
       .(round(p2,dig)))), pos=4, cex=0.8)
  }

  if (modes == 3) {
    legend(0.65*xmax+0.35*xmin,1.02*ylim[2],
      col=c("navyblue","red","orange","black"),lwd=c(2,2,2,1), lty=c(2,4,4,1),
       legend=c("Null","Nonnull 1","Nonnull 2","Mixture"), bty="n",cex=1)
  }
  else {
    legend(0.65*xmax+0.35*xmin,1.02*ylim[2],
      col=c("navyblue","red","black"),lwd=c(2,2,1), lty=c(2,4,1),
       legend=c("Null","Nonnull","Mixture"), bty="n",cex=1)
  }

##################
# Plot mg and the fitted posterior distribution:
#
  m = seq(min(mg),max(mg), length.out=length(mg))
  pst = posteriorIG(m,f,alpha_hat, beta_hat)
  brks = min(100,max(floor(G*0.1),20))
  histmax = max(hist(mg,breaks=40,plot=F)$density)
  hist(mg, breaks=brks, probability=TRUE,
    main=bquote(bold(paste("(b) mean(",m[g],")=",
    .(round(mean(mg),2)), " , s.d.(",m[g],")=",.(round(sd(mg),2))))),
    xlab=expression(m[g]),ylim=c(0,max(pst,histmax)),
    xlim=c(0,quantile(mg,mgq)))
  lines(m, pst, col=2, lwd=2)

  mtext(paste(titletext),
      NORTH<-3, line=0, adj=0.5, cex=2, col="red", outer=TRUE)

  text(0.7*quantile(mg,mgq),0.7*max(pst),bquote(
    paste(hat(alpha),"=",.(round(alpha_hat,2)),", ",hat(beta),
     "=",.(round(beta_hat,2)))))

##################
# create a volcano plot:
  plot(dg,pBH0,axes=F,xlab=expression(d[g]),ylab="Adjusted p-value",
    pch=20,main=bquote(bold(paste("(c) volcano plot (p-values)"))),col="grey")
  axis(1)
  axis(2)
  plot(dg,RRfdr0,axes=F,xlab=expression(d[g]),ylab=expression(p[0]),
    pch=20,main=bquote(bold(paste("(d) volcano plot (posterior probability)"))),
    col="grey")
  axis(1)
  axis(2)

  dev.off()

##################
  plotname2 = sprintf("%s/LEMMAplots2.ps", outdir)
  postscript(file=plotname2, width=10, height=7, paper="letter",horizontal=T,
   family="Times")
  par(oma=c(0,0,4,0),mfrow=c(1,1),mar=c(5,4.1,4.1,2.1))
  cols = seq(261,361,by=10)
  ts = c(0.00001,0.0001,0.001,0.01,0.02,0.05,0.1,0.2,0.5,1)
  ths = kronecker((rep(1,G)),t(ts))
  shades = 1+rowSums(pBH0-ths>0)
  L = length(dg)
  T = length(cols)
  ym = min(dg) - 0.15*(max(dg)-min(dg))
  plot(dg,pch=20,axes=F,col=colors()[cols[shades]],ylab=expression(d[g]),
   xlab="Gene #",ylim=c(ym,1.1*max(dg)),main=titletext,cex=(1.5-shades/T))
  axis(1)
  axis(2)
  abline(h=tau,col=2,lwd=2)
  text(1,1.1*tau,expression(hat(tau)),cex=1.3)

  for (i in 1:(T-1)) {
    rect((i-1)*L/T,ym,i*L/T,(ym+0.05*(max(dg)-min(dg))),col=colors()[cols[i]],
      border="white")
    text(0.5*((i-1)*L/T + i*L/T),(ym+0.08*(max(dg)-min(dg))),
      sprintf("p<=%s",ts[i]),cex=0.7)
  }
  minp = which.min(pBH0)
  text(minp,1.05*dg[minp],genename[minp],cex=0.7,col=3)
  points(minp,dg[minp],col=3,cex=1.5)
  dev.off()

  pstopdf1 = sprintf("ps2pdf %s",plotname1)
  pstopdf2 = sprintf("ps2pdf %s",plotname2)
  if (.Platform$OS == "windows") {
    shell(pstopdf1)
    shell(pstopdf2)
  } else { # assuming linux
    pdfplotname = sprintf("%s/LEMMAplots1.pdf", outdir)
    pstopdf = sprintf("ps2pdf %s %s",plotname1,pdfplotname)
    system(pstopdf)
    pdfplotname = sprintf("%s/LEMMAplots2.pdf", outdir)
    pstopdf = sprintf("ps2pdf %s %s",plotname2,pdfplotname)
    system(pstopdf)
  }
}

###################
#
# Output functions:
#
#   note: if you run these functions separately from the estimation,
#         load the AllData.RData file, first

# print the top genes using the local FDR criteria.
# The function expects the following parameters:
# type - either "RR" or "FDR"
# outdir - the output directory
# data - either the RR (RRfdr) or FDR (pBH) test statistics
# geneid - the gene IDs
# genename - the gene names
# topgenes - can be a positive integer less than or equal 
#   to the number of genes, or "all", or "nonnull")
# geneid - the gene IDs
# titletext - the title associated with the dataset
printTopGenes <- function(type,outdir,data0, data1, data2,
                    geneid,genename,topgenes="nonnull",titletext,
                    cutoff,modes=3) {
  resultsf <- file(sprintf("%s/results%s.txt",outdir,type),"w")
  cat("  Dataset: ",titletext,"\n\n", file=resultsf)
  # trim spaces in the beginning and end of gene names:
  genename = sub('^ +', '', genename)
  genename = sub(' +$', '', genename)
  shortgname = substr(genename,1,40)
  if (type == "RR") {
    if(topgenes == "nonnull") { topgenes = sum(data0<=cutoff) }
  }
  if (type == "FDR") {
    if(topgenes == "nonnull") { topgenes = sum(data0<=cutoff) }
  }
  if(topgenes == "all") { topgenes = length(data0) }
  geneorder = order(data0)
  orderedtable = data0[geneorder]
  cat(sprintf("Top %d genes, using the %s test statistics\n\n",
     topgenes,type), file=resultsf)
  if (type == "RR") {
    cat(sprintf("%62sPosterior probability\n"," "), file=resultsf)
    if (modes == 3) {
      cat(sprintf("%5sId            Gene%39sNull    Group1  Group2  Gene Effect\n"," "," "), file=resultsf)
    }
    else {
      cat(sprintf("%5sId            Gene%39sNull    Nonnull         Gene Effect\n"," "," "), file=resultsf)
    }
  }
  if (type == "FDR") {
    cat(sprintf("%5sId            Gene%39sGene effect  p-value(BH)\n"," "," "), file=resultsf)
  }
  if (topgenes == 0) {
    cat("No genes were found to be in the nonnull group for the
        given threshold\n",file=resultsf)
  } else {
    for (i in 1:topgenes) {
      if (type == "RR") {
        if (modes == 3) {
         cat(sprintf("%-5d%-14s%-43s%1.4f  %1.4f  %1.4f  %+1.4f\n",
          i, geneid[geneorder[i]], shortgname[geneorder[i]],
          data0[geneorder[i]],data1[geneorder[i]],
          1-(data0[geneorder[i]]+data1[geneorder[i]]),
          data2[geneorder[i]]), file=resultsf)
        }
        else {
         cat(sprintf("%-5d%-14s%-43s%1.4f  %1.4f          %+1.4f\n",
          i, geneid[geneorder[i]], shortgname[geneorder[i]],
          data0[geneorder[i]], 1-data0[geneorder[i]],
          data2[geneorder[i]]), file=resultsf)
        }
      }
      if (type  == "FDR") {
       cat(sprintf("%-5d%-14s%-43s%+1.4f      %1.4f\n",
        i, geneid[geneorder[i]], shortgname[geneorder[i]],
        data1[geneorder[i]],data0[geneorder[i]]), file=resultsf)
      }
    }
  }
  close(resultsf)
}

saveAsCsv <- function(outdir,loadData=TRUE) {

  if (loadData) {
    varsfile = sprintf("%s/AllData.RData",outdir)
    load(varsfile)
  }
  csvf <- file(sprintf("%s/results.csv",outdir),"w")
  cat("Gene,ID,dg,BH,p_0,p_1,p_2\n", file=csvf)
  genename = gsub('^ +', '', genename)
  genename = gsub(' +$', '', genename)
  genename = gsub(',', ';', genename)
  geneid = gsub('^ +', '', geneid)
  geneid = gsub(' +$', '', geneid)
  geneid = gsub(',', ';', geneid)
  cat(paste(geneid,genename,dg,pBH0,RRfdr0,RRfdr1,RRfdr2,sep=","),
      sep="\n",file=csvf)
  close(csvf)
}

printEst <- function(name, est, sd) {
  if (sd == "*") {
    s = sprintf("%s: %2.5f(%s)",name,est,sd)
  } else {
    s = sprintf("%s: %2.5f(%2.5f)",name,est,as.numeric(sd))
  }
  s
}

###################
#
# The main program:
#

lemma <- function(dataframe, locfdrcutoff=0.2, fdrcutoff=0.2, mgq=1,
    titletext="",outdir=tempdir(),topgenes="nonnull", tol=1e-6, maxIts=50000,
    modes=3, plots=TRUE, saveascsv=TRUE, ErrVarEst="MLE") {
  if (!file.exists(outdir)) {
    if (.Platform$OS == "windows") {
      tmpoutdir = sub('/', '\\\\', outdir)
      mkoutdir = sprintf("mkdir %s",tmpoutdir)
      shell(mkoutdir)
    } else { # assuming linux
      mkoutdir = sprintf("mkdir -p %s",outdir)
      system(mkoutdir)
    }
  } else { cat("Directory ",outdir," already exists.\n") }
  logf <- file(sprintf("%s/log.txt",outdir),"w")
  cat("  Dataset: ",titletext,"\n\n", file=logf)
  cat(" ",date(),"\n", file=logf)

  geneid = dataframe[,1]
  genename = dataframe[,2]
  Y1 = dataframe[,grep("Y1_",names(dataframe))]
  Y2 = dataframe[,grep("Y2_",names(dataframe))]
  if (dim(Y2)[2] == 0) {    # one treatment group / paired data
    Y2 = matrix(0,nrow=dim(Y2)[1],ncol=1)
  }

  cat("Fitting the model...\n")
  cat("Fitting the model...\n", file=logf)
  G  = dim(Y1)[1]
  n1 = dim(Y1)[2]
  n2 = dim(Y2)[2]
  f  = (n1-1) + (n2-1)
  dg = rowMeans(Y1) - rowMeans(Y2)
  mg = (rowSums((Y1 - rowMeans(Y1))^2) +
        rowSums((Y2 - rowMeans(Y2))^2)) / f  # mg = the MSEs

  # remove genes with extreme values of mg (mgq=1 means - use all genes):
  if (mgq < 1) {
    idx = setdiff(seq(1:G),which(mg>quantile(mg,mgq)))
    Y1 = Y1[idx,]
    Y2 = as.matrix(Y2[idx,])
    genename = genename[idx]
    geneid = geneid[idx]
    G  = dim(Y1)[1]
    dg = rowMeans(Y1) - rowMeans(Y2)
    mg = (rowSums((Y1 - rowMeans(Y1))^2) +
          rowSums((Y2 - rowMeans(Y2))^2)) / f  # mg = the MSEs
  }
  mgbar = mean(mg) 
  mgvar = var(mg)
  cat("  Number of genes (G)=",G,"\n", file=logf)
  cat("  Subjects in group 1 (n1)=",n1,"\n", file=logf)
  if (n2 == 1) { cat("    (one-treatment only)\n", file=logf) }
  else { cat("  Subjects in group 2 (n2)=",n2,"\n", file=logf) }
  cat("  Mean(d_g)=",mean(dg)," s.d.(d_g)=",sd(dg),"\n", file=logf)
  cat("  Mean(m_g)=",mgbar," s.d.(m_g)=",sqrt(mgvar),"\n", file=logf)

  cat("Estimating alpha, beta...\n", file=logf)

  if (ErrVarEst == "MM") { # Method of Moments (slightly faster, but always
                           # gives unique estimates, with alpha>2)
    cat("  using method of moments estimation\n", file=logf)
    M1 = mgbar
    M2 = mean(mg^2)
    alpha_hat = (2*M2-(1+2/f)*M1^2)/(M2-(1+2/f)*M1^2)
    beta_hat  = 1/(M1*(alpha_hat-1))
  }
  else { # MLE is the default
    cat("  using maximum likelihood estimation\n", file=logf)
    alphainit = 2+mgbar^2/mgvar
    betainit = 1/(mgbar*(alphainit-1))
    out = nlminb(start=c(alphainit,betainit), objective= mg.marginal.nloglik,
      lower = c(0.5,1e-6), upper = c(Inf,Inf), MSE=mg, f=f)
    if(out$convergence!=0) {
      cat("\n*** convergence problems in finding alpha_hat and beta_hat [",
        out$convergence,"]\n   ",out$message,"\n",file=logf)
    }
    alpha_hat = out$par[1]
    beta_hat = out$par[2]
  }
  M = varAlphaBeta(alpha_hat,beta_hat,f,mg)
  sig2eb = (f/2)*mg/( f/2 + alpha_hat +1 ) +
            alpha_hat*mgbar/( f/2 + alpha_hat +1 )

  cat("  Empirical Bayes estimates for alpha,beta:\n", file=logf)
  alphaest = printEst("alpha", alpha_hat, sqrt(M[1,1]))
  betaest = printEst("beta", beta_hat, sqrt(M[2,2]))
  cat("    ",alphaest,"\n     ",betaest,"\n",file=logf)
  cat(sprintf("  Fitted error variance: mean = %2.6f,  s.d. = %2.6f\n",
      mean(sig2eb),sd(sig2eb)), file=logf)
  if (modes == 3) {
    cat("Estimating tau, psi, sig2psi, p1, p2...\n", file=logf)
  }
  else { cat("Estimating tau, psi, sig2psi, p1...\n", file=logf) }

  em.out = get.EM.est3(dg, sig2eb, n1, n2,logf,tol, maxIts,modes)
  tau = em.out$tau
  psi = em.out$psi
  sig2psi = em.out$sig2psi
  p1 = em.out$p1
  p2 = em.out$p2
  p0 = 1-p1-p2
  
  INM = varPhi(sig2eb,n1,n2,sig2psi,p0,dg,tau,p1,psi,p2)
  if (p1+p2 < 1/(10*G)) {
    sds = c("*","*",1/INM[3,3],"*","*")
  }
  else {
    if (p1 < 1/(10*G)) {
      INM = INM[2:5,2:5]
      if (min(eigen(INM,only.values=TRUE)$values) > 0) {
        sds = c("*",sqrt(diag(solve(INM))))
      }
      else {
        sds = c("*","*","*","*","*")
        cat("  Information matrix not positive-definite.
                Cannot obtains standard errors for p1,p2,tau,psi,sig2psi\n",
            file=logf)
      }
    }
    else {
      if (p2 < 1/(10*G)) {
        INM = INM[c(1,3,4,5),c(1,3,4,5)]
        if (min(eigen(INM,only.values=TRUE)$values) > 0) {
          ests = sqrt(diag(solve(INM)))
          sds = c(ests[1],"*",ests[2:4])
        }
        else {
          sds = c("*","*","*","*","*")
          cat("  Information matrix not positive-definite.
                  Cannot obtains standard errors for p1,p2,tau,psi,sig2psi\n",
              file=logf)
        }
      }
      else {
        if (min(eigen(INM,only.values=TRUE)$values) > 0) {
          sds = sqrt(diag(solve(INM)))
        }
        else {
          sds = c("*","*","*","*","*")
          cat("  Information matrix not positive-definite.
                  Cannot obtains standard errors for p1,p2,tau,psi,sig2psi\n",
              file=logf)
        }
      }
    } 
  }
  
  tauest = printEst("tau", tau, sds[3])
  cat("       ",tauest,"\n",file=logf)
  psiest = printEst("psi", psi, sds[4])
  cat("       ",psiest,"\n",file=logf)
  sig2psiest = printEst("sig2psi", sig2psi, sds[5])
  cat("   ",sig2psiest,"\n",file=logf)
  p1est = printEst("p1", p1, sds[1])
  cat("        ",p1est,"\n",file=logf)
  if (modes == 3) {
    p2est = printEst("p2", p2, sds[2])
    cat("        ",p2est,"\n",file=logf)
  }

  cat("Calculating the posterior probabilities...\n", file=logf)
  if (n2 == 1)  {
    sig2g = sig2eb*(1/n1)
  }
  else {
    sig2g = sig2eb*(1/n1+1/n2)
  }
  sd0 = sqrt(sig2g)
  sd1 = sqrt(sig2g+sig2psi)
  RRfdr0 = posteriorProbability("null",(1-p1-p2),p1,p2, dg,tau,psi,sd0,sd1)
  RRfdr1 = posteriorProbability("nonnull1",(1-p1-p2),p1,p2, dg,tau,psi,sd0,sd1)
  RRfdr2 = 1 - (RRfdr0+RRfdr1)

  # Using FDR:
  # Using the Benjamini-Hochberg procedure (get all the p-values,
  # and declare as nonnull the ones with p-value less than or
  # equal to a user-defined threshold.
  pBH0 = p.adjust(2*pnorm(-abs(dg),mean=tau,sd=sqrt(sig2g)), method="BH")

  # print the top genes using the local FDR criteria:
  printTopGenes("RR", outdir, RRfdr0,RRfdr1,(dg-tau), geneid, genename,
         topgenes,titletext,locfdrcutoff,modes)
  printTopGenes("FDR", outdir, pBH0,(dg-tau),c(), geneid, genename,
         topgenes,titletext,fdrcutoff)

  cat("  Using locfdr cutoff of <=", locfdrcutoff," there are ",
       sum(RRfdr0<=locfdrcutoff)," non-null genes.\n", file=logf)
  cat("  Using FDR threshold of <=", fdrcutoff," there are ",
       sum(pBH0<=fdrcutoff)," non-null genes.\n", file=logf)
     
  varsfile = sprintf("%s/AllData.RData",outdir)
  save(geneid,genename,dg,mg,n1,n2,f,G,RRfdr0,RRfdr1,RRfdr2,
       alpha_hat,beta_hat,sig2eb,tau,psi,sig2psi,p1,p2,p0,pBH0,modes,
       file=varsfile)
  cat("Saved the following variables in ",varsfile," :
    dg,mg,n1,n2,f,G,RRfdr0,RRfdr1,RRfdr2, alpha_hat,beta_hat,sig2eb, 
    tau,psi,sig2psi,p1,p2,p0,pBH0,modes\n", file=logf)     
  if (saveascsv) {
    saveAsCsv(outdir)
  }
  if (plots) {
    cat("Creating diagnostics plots...\n")
    cat("Creating diagnostics plots...\n", file=logf)
    lemmaPlots(outdir, titletext, mgq=0.99, modes=modes)
  }
  cat("Done.\n", date(),"\n", file=logf)
  cat("\n\nDone.\n   Check output in ",outdir,"\n\n")
  close(logf)
}
