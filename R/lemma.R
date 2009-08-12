# lemma.R
# 2009-08-12
# version 1.2
# Haim Bar, Elizabeth Schifano {hyb2,eds27}@cornell.edu
#
# Allow for two nonnull groups
# Allow topgenes to be all, or nonnull (to print all genes, or
#   all the nonnulls, repectively
# Implmented a function to print the top genes separately from
#   the estimation phase
# Incorporated the plotting function into the main R file
# 

##############################
#
# the EM algorithm functions:
#

# assuming 3 modes (one null, two nonnull with the same variance, and
#   with opposite sign means)
get.EM.est3 <- function(Y1, Y2, dg, sig2, n1, n2,logf,tol,maxIts,modes=3){
  out = vector("list", 7)
  names(out) = c("sig2psi", "psi", "p1","p2","tau")

  # set initial values:
  tau.old = 0
  psi.old = mean(dg-tau.old)/2
  p1.old = min(.95,mean(var(dg-tau.old)>=((1/n1+1/n2)*sig2)))
  p1.old = max(0.05,p1.old)
  p1.old = 0.1
  sig2psi.old =1

  if (modes == 3) {
    p2.old = min(.95,mean(var(dg-tau.old)>=((1/n1+1/n2)*sig2)))
    p2.old = max(0.05,p2.old)
    p2.old = 0.1
  }
  else { p2.old = 0 }

  sd0 = sqrt((1/n1+1/n2)*sig2)
  sd1.old = sqrt((1/n1+1/n2)*sig2 + sig2psi.old)

  # w1 correspond to the indicators that are 1 if gene g is non-null #1
  w1.old = p1.old*dnorm(dg-tau.old, psi.old, sd1.old)/
     ( p1.old*dnorm(dg-tau.old, psi.old, sd1.old) +
       p2.old*dnorm(dg-tau.old, -psi.old, sd1.old) +
       (1-p1.old-p2.old)*dnorm(dg-tau.old, 0, sd0) )
  # w2 correspond to the indicators that are 1 if gene g is non-null #2
  if (modes == 3) {
    w2.old = p2.old*dnorm(dg-tau.old, -psi.old, sd1.old)/
     ( p1.old*dnorm(dg-tau.old, psi.old, sd1.old) +
       p2.old*dnorm(dg-tau.old, -psi.old, sd1.old) +
       (1-p1.old-p2.old)*dnorm(dg-tau.old, 0, sd0) )
  }
  else { w2.old = 0 }
  aerr  = 1
  count = 0
  
  while (aerr>tol){
    out.sig2psi = try(uniroot(f=Mstep.sig2psi, lower = -200,
        upper = 200, X=dg-tau.old,sig2=sig2, w1=w1.old,
        w2=w2.old,n1=n1, n2= n2, psi=(psi.old)),silent=T)
    if (!is.character(out.sig2psi))
    { sig2psi.cur = exp(out.sig2psi$root) }
    else { sig2psi.cur = min(sig2psi.old,.00001) }
    
    psi.cur = Mstep.psi(dg-tau.old, sig2, w1.old, w2.old,  n1, n2,
          sig2psi.cur)

    v0 = (1-w1.old-w2.old)/((1/n1+1/n2)*sig2)
    v1 = w1.old/(((1/n1+1/n2)*sig2)+sig2psi.cur)
    v2 = w2.old/(((1/n1+1/n2)*sig2)+sig2psi.cur)
    u1 = sum((v1-v2)*dg)/sum(v1+v2)
    u2 = sum(v1-v2)/sum(v1+v2)
    tau.cur = ( ( sum((v0+v1+v2)*dg) - (sum(v1-v2))*u1 ) /
                ( sum(v0+v1+v2) - (sum(v1-v2))*u2 ))

    sd1.cur= sqrt(sig2psi.cur + (1/n1+1/n2)*sig2)

    w1.cur = p1.old*dnorm(dg-tau.cur, psi.cur, sd1.cur)/
       ( p1.old*dnorm(dg-tau.cur, psi.cur, sd1.cur) +
         p2.old*dnorm(dg-tau.cur, -psi.cur, sd1.cur) +
         (1-p1.old-p2.old)*dnorm(dg-tau.cur, 0, sd0) )
    if (modes == 3) {
      w2.cur = p2.old*dnorm(dg-tau.cur, -psi.cur, sd1.cur)/
       ( p1.old*dnorm(dg-tau.cur, psi.cur, sd1.cur) +
         p2.old*dnorm(dg-tau.cur, -psi.cur, sd1.cur) +
         (1-p1.old-p2.old)*dnorm(dg-tau.cur, 0, sd0) )
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
  V= (1/n1+1/n2)*sig2 + sig2psi
  sum( (w1+w2)/V - (w1*(psi-X)^2 + w2*(psi+X)^2)/(V^2) )
}

Mstep.psi <- function(X, sig2, w1, w2, n1, n2, sig2psi){
  V = (1/n1+1/n2)*sig2 + sig2psi
  sum( ((w1-w2)/V)*X) / sum( (w1+w2)/V )
}

mg.marginal.nloglik <- function(pars, MSE, f){
  alpha = pars[1]
  beta = pars[2]
  -sum(log ((gamma(alpha+f/2)*(MSE^(f/2-1)) /
         ((1/beta+MSE*f/2)^(f/2+alpha)*(beta^alpha)*gamma(alpha))) ))
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
  plotname = sprintf("%s/LEMMAplots.ps", outdir)
  postscript(file=plotname, width=10, height=6, paper="letter",horizontal=T,
   family="Times")
  par(oma=c(0,0,4,0),mfrow=c(1,2),mar=c(5,4.1,4.1,2.1))
  brks = min(100,max(floor(G*0.1),20))
  x = seq(min(dg),max(dg),by=((max(dg)-min(dg))/100))
  ynull = p0*dnorm(sort(dg),tau,sqrt(((1/n1 + 1/n2)*sig2eb[order(dg)])))
  yalt1 = p1*dnorm(sort(dg),tau+psi,
          sqrt(sig2psi+((1/n1 + 1/n2)*sig2eb[order(dg)])))
  yalt2 = p2*dnorm(sort(dg),tau-psi,
          sqrt(sig2psi+((1/n1 + 1/n2)*sig2eb[order(dg)])))
  ymix = (p1*dnorm(x,tau+psi,sqrt(sig2psi + mean(((1/n1 + 1/n2)*sig2eb))))
         +p2*dnorm(x,tau-psi,sqrt(sig2psi + mean(((1/n1 + 1/n2)*sig2eb))))
         +p0*dnorm(x,tau,sqrt(mean(((1/n1 + 1/n2)*sig2eb)))))

  ymax = max(ynull,yalt1,yalt2,ymix,max(hist(dg,breaks=40,plot=F)$density))
  ylim = c(0,ymax*1.2)
  maxsig = max(sig2psi,sig2psi)+ mean((1/n1 + 1/n2)*sig2eb)
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
  lines(x,p0*dnorm(x,tau,sqrt(mean(((1/n1 + 1/n2)*sig2eb)))),col="navyblue",
     lwd=3,lty=2)
  lines(x,p1*dnorm(x,tau+psi,sqrt(sig2psi + mean(((1/n1 + 1/n2)*sig2eb)))),
     col="red", lwd=2, lty=4)
  if (modes == 3) {
    lines(x,p2*dnorm(x,tau-psi,sqrt(sig2psi + mean(((1/n1 + 1/n2)*sig2eb)))),
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

  dev.off()
  pstopdf = sprintf("ps2pdf %s",plotname)
  if (.Platform$OS == "windows") {
    shell(pstopdf)
  } else { # assuming linux
    pdfplotname = sprintf("%s/LEMMAplots.pdf", outdir)
    pstopdf = sprintf("ps2pdf %s %s",plotname,pdfplotname)
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

###################
#
# The main program:
#

lemma <- function(dataframe, locfdrcutoff=0.2, fdrcutoff=0.2, mgq=1,
    titletext="",outdir="OUT",topgenes="nonnull", tol=1e-6, maxIts=50000,
    modes=3, plots=TRUE) {
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
  idx = setdiff(seq(1:G),which(mg>quantile(mg,mgq)))
  Y1 = Y1[idx,]
  Y2 = Y2[idx,]
  genename = genename[idx]
  geneid = geneid[idx]
  G  = dim(Y1)[1]
  n1 = dim(Y1)[2]
  n2 = dim(Y2)[2]
  f  = (n1-1) + (n2-1)
  dg = rowMeans(Y1) - rowMeans(Y2)
  mg = (rowSums((Y1 - rowMeans(Y1))^2) +
        rowSums((Y2 - rowMeans(Y2))^2)) / f  # mg = the MSEs
  mgbar = mean(mg) 
  mgvar = var(mg)
  cat("  Number of genes (G)=",G,"\n", file=logf)
  cat("  Subjects in group 1 (n1)=",n1,"\n", file=logf)
  cat("  Subjects in group 2 (n2)=",n2,"\n", file=logf)
  cat("  Mean(d_g)=",mean(dg)," s.d.(d_g)=",sd(dg),"\n", file=logf)
  cat("  Mean(m_g)=",mgbar," s.d.(m_g)=",sqrt(mgvar),"\n", file=logf)

  cat("Estimating alpha, beta...\n", file=logf)

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
  sig2eb = (f/2)*mg/( f/2 + alpha_hat +1 ) +
            alpha_hat*mgbar/( f/2 + alpha_hat +1 )

  cat(sprintf("  Empirical Bayes estimates for alpha,beta: = %2.5f,%2.5f\n",
      alpha_hat,beta_hat), file=logf)
  cat(sprintf("  Fitted error variance: mean = %2.6f,  s.d. = %2.6f\n",
      mean(sig2eb),sd(sig2eb)), file=logf)
  if (modes == 3) {
    cat("Estimating tau, psi, sig2psi, p1, p2...\n", file=logf)
  }
  else { cat("Estimating tau, psi, sig2psi, p1...\n", file=logf) }

  em.out = get.EM.est3(Y1, Y2, dg, sig2eb, n1, n2,logf,tol, maxIts,modes)
  tau = em.out$tau
  psi = em.out$psi
  sig2psi = em.out$sig2psi
  p1 = em.out$p1
  p2 = em.out$p2
  p0 = 1-p1-p2
  cat("  hat(tau) = ",tau,"\n", file=logf)
  cat("  hat(psi) = ",psi,"\n", file=logf)
  cat("  hat(sig2psi) = ",sig2psi,"\n", file=logf)
  cat("  hat(p1) = ",p1,"\n", file=logf)
  if (modes == 3) {
    cat("  hat(p2) = ",p2,"\n", file=logf)
  }
  
  cat("Calculating the posterior probabilities...\n", file=logf)
  sig2g = sig2eb*(1/n1+1/n2)
  RRfdr0 = p0*dnorm(dg-tau, mean=0, sd=sqrt(sig2g))/
            (p1*dnorm(dg-tau, mean=psi, sd=sqrt(sig2g+sig2psi)) +
             p2*dnorm(dg-tau, mean=-psi, sd=sqrt(sig2g+sig2psi)) +
             p0*dnorm(dg-tau, mean=0, sd=sqrt(sig2g)) )
  RRfdr1 = p1*dnorm(dg-tau, mean=psi, sd=sqrt(sig2g+sig2psi))/
            (p1*dnorm(dg-tau, mean=psi, sd=sqrt(sig2g+sig2psi)) +
             p2*dnorm(dg-tau, mean=-psi, sd=sqrt(sig2g+sig2psi)) +
             p0*dnorm(dg-tau, mean=0, sd=sqrt(sig2g)) )
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
  save(dg,mg,n1,n2,f,G,RRfdr0,RRfdr1,RRfdr2, alpha_hat,beta_hat,sig2eb, 
       tau,psi,sig2psi,p1,p2,p0,pBH0,modes, file=varsfile)
  cat("Saved the following variables in ",varsfile," :
    dg,mg,n1,n2,f,G,RRfdr0,RRfdr1,RRfdr2, alpha_hat,beta_hat,sig2eb, 
    tau,psi,sig2psi,p1,p2,p0,pBH0,modes\n", file=logf)     
  if (plots) {
    cat("Creating diagnostics plots...\n")
    cat("Creating diagnostics plots...\n", file=logf)
    lemmaPlots(outdir, titletext, mgq=0.99, modes=modes)
  }
  cat("Done.\n", date(),"\n", file=logf)
  cat("\n\nDone.\n   Check output in ",outdir,"\n\n")
  close(logf)
}
