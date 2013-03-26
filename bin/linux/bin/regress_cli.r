#!/usr/bin/env Rscript

##############################################################################
## Functions

stdAnalysis = function(obs.infile,
        sim.infile,
        out.file="posterior-summary.txt",
        adjusted_path="adjusted-posterior-samples.txt",
        tol=1,
        rejmethod=F,
        stat_prefixes=c("pi","wattTheta","pi.net","tajD.denom"),
        continuous_prefixes=c("PRI.E.t", "PRI.omega"),
        discrete_prefixes=c("PRI.Psi", "PRI.model"),
        stat_indices = NULL,
        discrete_indices = NULL,
        continuous_indices = NULL) {
    header = parse_header(sim.infile)
    simDat <- getData(sim.infile)

    nPairs <- simDat[["numTaxonPairs"]]
    params.from.priorDistn <- simDat[["prior.names"]]
    summary.stat.names <- simDat[["summary.stats"]]
    simDat <- simDat[["dat"]]

    models = NULL
    model_idx = grep('PRI.model', names(simDat), ignore.case=TRUE)
    if (length(model_idx) == 1) {
        models = unique(simDat[,model_idx])
    }
    if (length(models) < 2) {
        discrete_prefixes = discrete_prefixes[grep(
                'model',
                discrete_prefixes,
                invert=T,
                ignore.case=T)]
        discrete_indices = discrete_indices[discrete_indices != model_idx]
    }

    # if tol is NA, set default tol to get 1000 best matches.
    if (is.na(tol)) {
      tol <- 1000/nrow(simDat)
    }

    # construct the column names
    usedColNames = NULL
    if (is.null(stat_indices)) {
        usedColNames <- as.vector(sapply(stat_prefixes, paste,
                                     1:nPairs, sep="."))
    } else {
        usedColNames = names(simDat)[stat_indices]
    }

    #load OBSERVED summary stat vector
    obsDat <-getData(obs.infile, header=header)
    if (obsDat[["numTaxonPairs"]] != nPairs) {
      cat("ERROR: The number of taxon pairs are not same between the\n      ",
          "observed data,", obsDat$numTaxonPairs, "pairs, and simulated data,",
          nPairs, "pairs.\n")
      return(NA)
    }
    obsDat <- obsDat[["dat"]]

    # acceptance/regression, .... ie  the meat
    # The 1st column is PRI.numTauClass, which should be removed from analysis
    result = list(prior.names=params.from.priorDistn)
    dummy_idx = grep('PRI.numTauClass', params.from.priorDistn, ignore.case=TRUE)
    if (length(dummy_idx) > 0) {
        stopifnot(dummy_idx == c(1))
        result <- list(prior.names=params.from.priorDistn[-dummy_idx])
        if (!is.null(continuous_indices)) {
            continuous_indices = continuous_indices[continuous_indices != 1]
            continuous_indices = continuous_indices - 1
        }
        if (!is.null(discrete_indices)) {
            discrete_indices = discrete_indices[discrete_indices != 1]
            discrete_indices = discrete_indices - 1
        }
    }

    noPsiAnalysis <- F
    constrained = F
    prior.names = result$prior.names
    if (is.null(continuous_indices)) {
        continuous_patterns = gsub('[.]', '[.]', continuous_prefixes)
        continuous_indices = get_indices_of_patterns(
                target_vector = prior.names,
                patterns = continuous_patterns,
                ignore_case = TRUE,
                sort_indices = TRUE)
    }
    if (is.null(discrete_indices)) {
        discrete_patterns = gsub('[.]', '[.]', discrete_prefixes)
        discrete_indices = get_indices_of_patterns(
                target_vector = prior.names,
                patterns = discrete_patterns,
                ignore_case = TRUE,
                sort_indices = TRUE)
    }
    prior.names.cont = prior.names[continuous_indices]
    prior.names.discrete = prior.names[discrete_indices]

    # prior.names.pretty = list(PRI.Psi="psi",
    #   		    PRI.omega="omega",
    #   		    PRI.E.t="tau_mean",
    #   		    PRI.var.t="tau_var")
    # set min, max boundary of prior parameters, and verbose print message
    # PRI.omega >= 0, PRI.E.t >= 0, PRI.Psi > 0
    min.val <- list(PRI.Psi = 1, PRI.var.t = 0, PRI.E.t = 0, PRI.omega = 0 )
    max.val <- list(PRI.Psi = nPairs)
    verbose.print <- list(PRI.omega = "(=Var(t)/E(t))",
                          PRI.Psi="(= number of possible divtimes)",
                          PRI.E.t="(= E(t))")
    # PRI.Tau* >= 0
    tauNames <- params.from.priorDistn[grep("^PRI[.]Tau", params.from.priorDistn)]
    if(length(tauNames) > 0) {
      temp.val <- sapply(rep(0, length(tauNames)), list)
      names(temp.val) <- tauNames
      min.val <- c(min.val, temp.val)
    }
    # 1 <= PRI.Psi.* <= nPairs - (numTauClasss - 1)
    if(constrained) {
      psiNames <- params.from.priorDistn[grep("^PRI[.]Psi[.]", params.from.priorDistn)]
      temp.val <- sapply(rep(1, length(psiNames)), list)
      names(temp.val) <- psiNames
      min.val <- c(min.val, temp.val)
      # I could use length of psiNames below, instead of simDat[1,1] to get numTauClass
      temp.val <- sapply(rep(nPairs - (simDat[1,"PRI.numTauClass"]-1), length(psiNames)), list)
      names(temp.val) <- psiNames
      max.val <- c(max.val, temp.val)
      
      temp.val <- sapply(rep("(= number of taxon pairs that divergence at corresponding tau)", length(psiNames)), list)
      names(temp.val) <- psiNames
      verbose.print <- c(verbose.print, temp.val)
    }
    
    # run makepdANY for each para for continuous
    for (i in 1:length(prior.names.cont)) {
      thisPriorName <- prior.names.cont[i]
      # might need to work on constrained vs unconstrained here
      temp <- list(makepdANY(as.vector(obsDat[1,usedColNames],mode="numeric"),
                             simDat[,thisPriorName], simDat[,usedColNames], tol,
                             rep(T,len=nrow(simDat)),rejmethod=rejmethod))
      names(temp) <- thisPriorName
      
      # absorbing boundary
      if ( !  is.null(min.val[[thisPriorName]])) {
        temp[[thisPriorName]]$x[which(temp[[thisPriorName]]$x<min.val[[thisPriorName]])] <- min.val[[thisPriorName]]
      }
      if ( !  is.null(max.val[[thisPriorName]])) {
        temp[[thisPriorName]]$x[which(temp[[thisPriorName]]$x>max.val[[thisPriorName]])] <- max.val[[thisPriorName]]
      }

      result <- c(result, temp)
    }

    # deal with the discrete priors
    calmod.fail <- c()
    if(length(prior.names.discrete) > 0) {
      for (i in 1:length(prior.names.discrete)) {
        thisPriorName <- prior.names.discrete[i]
        calmod.res <- try(calmod(as.vector(obsDat[1,usedColNames],mode="numeric"),
                                 simDat[,thisPriorName], simDat[,usedColNames], tol,
                                 rep(T,len=nrow(simDat)),rejmethod=rejmethod))
        if(class(calmod.res) == "try-error") {
          calmod.res <- calmod(as.vector(obsDat[1,usedColNames],mode="numeric"),
                               simDat[,thisPriorName], simDat[,usedColNames], tol,
                               rep(T,len=nrow(simDat)),rejmethod=T)
          calmod.fail <- c(calmod.fail, thisPriorName)
          this.failed <- T
        } else {
          this.failed <- F
        }
        temp <- list(calmod.res)
        names(temp) <- thisPriorName
        
        if (rejmethod || this.failed) {
          # with simple rejection, $x contains the accepted values
          # Need to copy to $vals to make the later handling easier.
          temp[[thisPriorName]]$vals <- temp[[thisPriorName]]$x
        }
        result <- c(result, temp)
      }
    }
    
    post.modes = list()
    for (p in prior.names.cont) {
        m = try(loc1stats(result[[p]]$x,prob=0.95),silent=T)
        if (class(m) == "try-error") {
            post.modes[p] = "None"
        } else {
            post.modes[p] = m[1]
        }
    }
    adjusted_samples = list()
    sink(out.file)
    for (p in prior.names.discrete) {
        pname = sub("[.]", "_", sub("PRI[.]", "", p))
        values = c(1:100)
        if (regexpr('psi', p, ignore.case=TRUE)[1] != -1) {
            values = c(1:nPairs)
        }
        if (regexpr('model', p, ignore.case=TRUE)[1] != -1) {
            values = models
        }
        cat("[", pname, "]\n", sep="")
        if (p %in% calmod.fail) {
            cat("failed = True\n")
        } else {
            cat("failed = False\n")
            post.prob.matrix = result[[p]]$x2
            post.prob.vector = as.vector(post.prob.matrix)
            names(post.prob.vector) = colnames(post.prob.matrix)
            m = as.numeric(names(post.prob.vector)[which.max(post.prob.vector)])
            cat("\t[[adjusted_results]]\n")
            cat("\tmode = ", m, "\n", sep="")
            cat("\t\t[[[post_probs]]]\n")
            write_probabilites(post.prob.vector, values)
        }
        # m = as.numeric(colnames(post.probs)[which.max(post.probs)])
        # cat("mode = ", m, "\n", sep="")
        raw.post.counts = table(simDat[,p])
        raw.post.prob.table = raw.post.counts / sum(raw.post.counts)
        m = as.numeric(names(raw.post.prob.table)[which.max(raw.post.prob.table)])
        cat("\t[[unadjusted_results]]\n")
        cat("\tmode = ", m, "\n", sep="")
        cat("\t\t[[[post_probs]]]\n")
        write_probabilites(raw.post.prob.table, values)
    }
    for (p in prior.names.cont) {
        pname = sub("[.]", "_", sub("PRI[.]", "", p))
        cat("[", pname, "]\n", sep="")
        cat("mode = ", post.modes[[p]], "\n", sep="")
        cat("mean = ", mean(result[[p]]$x), "\n", sep="")
        cat("median = ", median(result[[p]]$x), "\n", sep="")
        quants = quantile(result[[p]]$x, prob=c(0.025,0.975))
        cat("\t[[quantiles]]\n")
        cat("\t'0.025' = ", quants[[1]], "\n", sep="")
        cat("\t'0.975' = ", quants[[2]], "\n", sep="")
        if (p == "PRI.omega") {
            omega.post = result[[p]]$x
            omega.prob = length(omega.post[omega.post<0.01]) / length(omega.post)
            cat("post_prob_zero = ", omega.prob, "\n", sep="")
        }
        adjusted_samples[[p]] = result[[p]]$x
    }
    cat("[settings]\n")
    cat("stats_used = ")
    cat(usedColNames, sep=', ')
    cat("\n")
    cat("continuous_parameters = ")
    cat(prior.names.cont, sep=', ')
    cat("\n")
    cat("discrete_parameters = ")
    cat(prior.names.discrete, sep=', ')
    cat("\n")
    sink()
    adj_samples = data.frame(adjusted_samples)
    write.table(adj_samples,
            adjusted_path,
            sep='\t',
            row.names=F,
            col.names=T,
            quote=F,
            fileEncoding='UTF-8')
}

write_probabilites = function(probs, n) {
    for (i in n) {
        val.str = as.character(i)
        if (val.str %in% names(probs)) {
            cat("\t\t", i, " = ", probs[val.str], "\n", sep="")
        } else {
        cat(i, " = 0.0\n", sep="")
        }
    }
}
# This function takes an file name as an argument, read in the data
# from the file, and assign appropriate names.
# Returns a list(dat, numTaxonPairs)
# dat is the data.frame, numTaxonPairs is an integer indicating the number
# of taxon pairs.

parse_header = function(infile) {
    header = scan(infile, what="character", nlines=1, quiet=T)
    return(header)
}

getData <- function (infile, header=NULL) {
    first.line = header
    skip = 0
    if (is.null(header)) {
        first.line <- scan(infile, what="character", nlines=1, quiet=T) #header
        skip = 1
    }
    dat <- scan(infile, skip=skip, quiet=T)
    dat <- data.frame(matrix(dat, ncol=length(first.line), byrow=T))
    names(dat) <- first.line  # assign the column names to the data.frame

    prior.names <- first.line[grep("^PRI[.]", first.line)]
    num.prior <- length(prior.names)
    # sum stats column-names (header) have the following form
    #   c("pi.b.1", "pi.b.2", "pi.b.3", "pi.w.1", "pi.w.2", "pi.w.3", ...)
    # Here, I'm getting rid of .digits part and taking unique names.
    sum.stat.names <- unique(sub("[.][0-9]+$", "", first.line[(num.prior+1):length(first.line)],fixed=F))  

    # number of taxon pairs can be calculated from
    nTaxPairs = (ncol(dat) - num.prior) / (length(sum.stat.names))

    return (list(dat=dat, numTaxonPairs=nTaxPairs, prior.names=prior.names, summary.stats=sum.stat.names))
}


# -----------------------------------------------------------------------
# 2D Kernel density estimates: q1 X q2
# -----------------------------------------------------------------------
# horizontal view
plotKernDensity <- function (res1, res2, title="q1 and q2", xlab="q1", ylab="q2") {
    bwq1 <- try(dpik(res1$x),silent=T)
    bwq2 <- try(dpik(res2$x),silent=T)
    # I think bandwidth choice by dpik may fail if there aren't enough unique
    # values (e.g. mostly 0), I'm not sure following is ok or not, but
    # bw.nrd0 seems to be more robust
    if(class(bwq1) == "try-error") {
      cat("INFO: In plotKernDensity(), simpler bandwidth selection used for ",
              xlab, "\n", file=stderr())
      bwq1 <- bw.nrd0(res1$x)
    }
    if(class(bwq2) == "try-error") {
      cat("INFO: In plotKernDensity(), simpler bandwidth selection used for ",
              ylab, "\n", file=stderr())
      bwq2 <- bw.nrd0(res1$x)
    }  
    x <- cbind(res1$x, res2$x)
    est <- bkde2D(x, bandwidth=c(bwq1,bwq2))
    par(ask = FALSE)
    persp(est$x1, est$x2, est$fhat, theta = 145, phi = 25, col = "lightblue1",
          xlab=xlab, ylab =ylab,
          zlab = paste("Pr(", xlab, ", ", ylab, "| X)",sep=""),
          main = paste("Joint Density of", title), axes = TRUE, nticks = 5,
          ticktype = "detailed", ltheta = -135, lphi = 145, shade = TRUE)
}


make.hist <-function(vect, res.makepd, title="", xlim, ...) {
      #old.mfcol <- par()$mfcol
    #  par(mfcol=c(3,1))
    bw_vect<-max(vect)/100
    bw_res.makepd<-max(res.makepd$x)/50
      hist.col = "white"
      if(missing(xlim)) {
        hist(vect,col=hist.col,border ="white",xlim=c(0,max(vect)),ylim=c(0,max(density(vect,bw=bw_vect)$y,density(res.makepd$x,bw=bw_res.makepd)$y)),prob=TRUE,main=paste(title),xlab=title, ...)
      } else {
        hist(vect,col=hist.col,freq=F,xlim=xlim,prob=TRUE,main=paste(title),
           xlab=title, ...)
      }
      lines(density(vect,bw=bw_vect),lty=2,col="red",pch=3)
    #  if(missing(xlim)) {
    #    hist(res.makepd$x,col="blue",prob=TRUE,freq=F,xlim=c(0,max(vect)),
    #         main=paste(title, ": Posterior Dist'n "), xlab=title, ...)
    #  } else {
    #    hist(res.makepd$x,col=hist.col,xlim=xlim,freq=F,prob=TRUE,
    #         main=paste(title, ": Posterior Dist'n "), xlab=title, ...)
    #  
      lines(density(res.makepd$x,bw=bw_res.makepd),lty=1,col="blue",pch=3)
      #par(mfcol=old.mfcol)
    legend("topright",c("Prior","Posterior"),text.col=c("red","blue"),lty=c(2,1),col=c("red","blue"))
}

plot.bf <- function(prior,posterior,...) {
    bf.out <- make.bf.vect(prior,posterior)
    if (is.null(bf.out)) {
      return(NULL);
    }
    # finding some reasonable y.lim for the plot, all values here are arbitrary
    y.max <- 40
    if (max(bf.out$bf) * 1.1 < 40) {
      y.max <- max(bf.out$bf)* 1.1
    } else if (mean(bf.out$bf > 40) > 0.6) {
      y.max <- sort(bf.out$bf)[ceiling(length(bf.out$bf) * 0.4)] * 1.1
    }
    plot(bf.out$crit.vals, bf.out$bf, type="l",
         ylim=c(0,y.max), 
         xlab="Hyper Parameter Thresholds",ylab="Bayes Factor",...)
    line.type = 3
    abline (1, 0, lty=2)
    abline(3, 0, lty=line.type)
    abline(10, 0, lty=line.type)
    abline(1/3, 0, lty=line.type)
    abline(1/10, 0, lty=line.type)
}

# Compares two Model of x < crit.val and x >= crit.val, and
# calculate the bayes factor.
# It returns a data.frame with two columns
make.bf.vect <- function(prior, posterior, crit.vals = NULL, num.points=100) {

  prior.len <- length(prior)
  posterior.len <- length(posterior)
  if (prior.len * posterior.len == 0) {  # some error
    return(NULL);
  }
  
  if (is.null(crit.vals)) {

    post.range <- range(posterior)

    # A single value of posterior is observed (e.g. at the boundary)
    if(post.range[1] == post.range[2]) {
      if(post.range[1] == 0) {
        post.range <- c(-0.01, 0.01)
      } else  {
        post.range <- c(0.99, 1.01) * post.range[1]
      }
    }

    # following is probably not required
    pri.min <- min(prior)
    if (post.range[1] <= pri.min) {
      if(pri.min >= post.range[2]) {  # impossible, but better check
        cat("WARN: In make.bf.vect(), weird prior/posterior range encountered\n",
            file=stderr())
        return(NULL)
      }
      post.range[1] <- pri.min
    }
    crit.vals <- seq(post.range[1], post.range[2], length.out=num.points+2)
    # get rid of the ends to avoid bf = Inf
    crit.vals <- crit.vals[2:(length(crit.vals) - 1)]
  }
  
  if (is.null(crit.vals)) {
    temp.pos <- posterior[posterior != max(posterior) &
                          posterior != min(posterior)]
     post.range <- range(temp.pos)
     crit.vals <- seq(post.range[1], post.range[2], length.out=num.points)
  }
  
  bf.vect <- numeric (0)
  for(i in 1:length(crit.vals)) {
    prior.below <- length(which(prior < crit.vals[i]))
    posterior.below <- length(which(posterior < crit.vals[i]))
    if (prior.below == 0 || posterior.below == posterior.len) {
      this.bf <- Inf
    } else {
      this.bf <- posterior.below * (prior.len-prior.below)/
        ((posterior.len - posterior.below) * prior.below)
    }
    bf.vect <- c(bf.vect, this.bf)
  }

  return (data.frame(crit.vals=crit.vals, bf=bf.vect))
}

# Marging two tables created by table()
# this is not generic, only for 1 dimensional table
merge.2tbl.byName <- function(arr1, arr2) {
  #this.name <- dimnames(arr1)[[length(dimnames(arr1))]]
  this.name <- names(arr1)
  d1 <- data.frame(cbind(this.name, arr1))
  names(d1)[1] <- 'name'
  
  #this.name <- dimnames(arr2)[[length(dimnames(arr2))]]
  this.name <- names(arr2)
  d2 <- data.frame(cbind(this.name, arr2))
  names(d2)[1] <- 'name'
  
  m1 <- merge(d1,d2,by='name',all=T)
  m1 <- m1[order(as.numeric(levels(m1$name)[m1$name])),]  #sorting rows
  m1 <- as.matrix(m1)
  m1[is.na(m1)] <- 0  # replacing NA with 0
  result <- as.table(t(m1[,2:3]))
  colnames(result) <- m1[,1]
  # convert character table to numeric
  return(type.convert(result))
}

##### from make_pd2005.r #####
makepdANY <- function(target,x,sumstat,tol,gwt,rejmethod=T)
{
# target is the set of target summary stats
# x is the parameter vector (long vector of numbers from the simulations) and is the dependent variable for the regression
# sumstat is an array of simulated summary stats (i.e. independent variables).
# NBB this function assumes 4 summary stats - I edited by hand for other numbers
# tol is the required proportion of points nearest the target values
# gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
# if rejmethod=T it doesn't bother with the regression, and just does rejection.


# If rejmethod=F it returns a list with the following components:-

# $x regression adjusted values
# $vals - unadjusted values in rejection region (i.e. normal rejection)
# $wt - the regression weight (i.e. the Epanechnikov weight)
# $ss - the sumstats corresponding to these points
# $predmean - estimate of the posterior mean
# $fv - the fitted value from the regression

if(missing(gwt))gwt <- rep(T,length(sumstat[,1]))

nss <- length(sumstat[1,])


# scale everything 

    scaled.sumstat <- sumstat

    for(j in 1:nss){
    	scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
    }
    target.s <- target

    for(j in 1:nss){

    	target.s[j] <- normalise(target[j],sumstat[,j][gwt])
    }

# calc euclidean distance

    sum1 <- 0
    for(j in 1:nss){
    	sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
   }
   dst <- sqrt(sum1)
# includes the effect of gwt in the tolerance
    dst[!gwt] <- floor(max(dst[gwt])+10)


# wt1 defines the region we're interested in 
    abstol <- quantile(dst,tol)
# making sure all simulated results are included when tol = 1.
    if(tol == 1) {abstol <- abstol * 1.1}
    wt1 <- dst < abstol

    if(rejmethod){
        l1 <- list(x=x[wt1],wt=0)
    }
    else{
        # this is proportional to k_delta(t) expression (5) of
        # Beaumont et al (2002)
        regwt <- 1-dst[wt1]^2/abstol^2

        # Note that scaled.sumstat is used instead of difference between
        # scaled.sumstat and target.s.  Since target.s (constant) is
        # subtracted from each simulated sumstat in Beaumont et al (2002),
        # we do not need to include it in the lsfit below.  The estimated
        # slopes are the same.  However intercept needs to be adjusted
        # below (see predmean).  predmean = hat(alpha) of the paper.
        fit1 <- lsfit(scaled.sumstat[wt1,],x[wt1],wt=regwt)
        predmean <- fit1$coeff %*% c(1,target.s)

        # x, correspond to phi*_i (between expressions (3) and (4) of the paper)
        l1 <- list(x=fit1$residuals+predmean,vals=x[wt1],wt=regwt,ss=sumstat[wt1,],predmean=predmean,fv = x[wt1]-fit1$residuals)

    }
    l1
}


normalise <- function(x,y){

if(var(y) == 0)return (x - mean(y))
(x-(mean(y)))/sqrt(var(y))
}

##### from loc2plot.r #####
loc2plot <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,...)
{
        sc1 <- sqrt(var(x))
        sc2 <- sqrt(var(y))
        if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        maxk=100,mint=100,cut=0.8,maxit=100)
        else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        xlim=xlim,maxk=100,mint=100,cut=0.8,maxit=100)
        lev <- sort(fitted(fit))[floor(cprob*length(x))]
        plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
        xlim=gxlim,...)
	
}

loc2plotw <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,wt,...)
{
        sc1 <- sqrt(var(x))
        sc2 <- sqrt(var(y))
        if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        maxk=100,mint=100,cut=0.8,maxit=100,weight=wt)
        else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        xlim=xlim,maxk=100,mint=100,cut=0.8,maxit=100,weight=wt)
        lev <- sort(fitted(fit))[floor(cprob*length(x))]
        plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
        xlim=gxlim,...)
}

gethpdprob2 <- function(x,y,px,py,alpha=0.5,xlim,gxlim,...)
{
        sc1 <- sqrt(var(x))
        sc2 <- sqrt(var(y))
        if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        maxk=100,mint=100,cut=0.8,maxit=100)
        else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        xlim=xlim,maxk=100,mint=100,cut=0.8,maxit=100)
#       d1 <- (x-px)^2+(y-py)^2
#       best <- d1 == min(d1)
#       lev <- mean(fitted(fit)[best])
        lev <- predict.locfit(fit,list(px,py))
        slev <- sort(fitted(fit))
        indic <- slev <= lev
        sum(indic)/length(x)
}
        
        


loc2mode <- function(x,y,alpha=0.5,xlim,...)
{
        sc1 <- sqrt(var(x))
        sc2 <- sqrt(var(y))
        if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        maxk=100,mint=100,cut=0.8,maxit=100)
        else fit <- locfit(~x+y,alpha=alpha,scale=c(sc1,sc2),
        xlim=xlim,maxk=100,mint=100,cut=0.8,maxit=100)
        tt <- max(fitted(fit))
        wt <- fitted(fit) == tt
        c(x[wt][1],y[wt][1])
}


loc1stats <- function(x,prob=0.95,alpha=0.5,xlim,...)
{
        if(missing(xlim)) fit <- locfit(~x,alpha=alpha)
        else fit <- locfit(~x,alpha=alpha,xlim=xlim)
        fx <- fitted(fit)
        x.modef <- max(fx)
        x.mode <- x[fx == x.modef]
        if(length(x.mode)>1)x.mode <- x.mode[1]
        lev <- sort(fx)[floor(prob*length(x))]
       
        cat("INFO: loc1stats: difference in log(prob. density) between the ",
            "max (mode) and the ends of confidence interval (at p=", prob,
            ") is ", log(x.modef)-log(lev), ".\n", sep="")
		
        l1 <- list()
        l1[[1]] <- x.mode
        indx <- order(x)
        ii <- 2
        flip <- TRUE
        for(j in 2:length(x)){
                if(flip && fx[indx[j]] > lev){
                        l1[[ii]] <- x[indx[j-1]]
                        flip <- FALSE
                        ii <- ii+1
                }
                else if(!flip && fx[indx[j]] < lev){
                        l1[[ii]] <- x[indx[j]]
                        flip <- TRUE
                        ii <- ii+1
                }
        }
        as.numeric(l1)
}
                        
                

tloc2plot <- function(x,y,cprob=0.5,alpha=0.5,xlim,gxlim,...)
{
        sc1 <- sqrt(var(x))
        sc2 <- sqrt(var(y))
        if(missing(xlim)) fit <- locfit(~x+y,alpha=alpha)
        else fit <- locfit(~x+y,alpha=alpha,
        xlim=xlim)
        lev <- sort(fitted(fit))[floor(cprob*length(x))]
        plot(fit,lev=lev,m=100,label=paste(as.character(100*(1-cprob)),"%",sep=""),
        xlim=gxlim,...)
}


##### from calmod.r #####
## from popABC, Aug 05, 2009 version
## Naoki fixed that $x2 prints correct categories instead of mu1, mu2,...

calmod <- function(target,x,sumstat,tol,gwt,rejmethod=T)
{
# this function uses categorical regression to estimate the posterior probability of a particular model 
#      under the ABC framework P(Y=y | S = s)
#
# target is the set of target summary stats - i.e. s, what you measured from the data.
# x is the parameter vector, Y (long vector of numbers from the simulations) and is the dependent variable for the regression
#          This is a categorical variable, taking values 1, 2, .. n where there are n categories (models)
# sumstat is an array of simulated summary stats, S (i.e. independent variables).
# tol is the required proportion of points nearest the target values
# gwt is a vector with T/F weights, weighting out any 'bad' values (determined by the simulation program - i.e. nan's etc)
# if rejmethod=T it doesn't bother with the regression, and just does rejection.


# If rejmethod=F it returns a list with the following components:-

# $x1 expected value on a logit scale, with standard errors of the estimate
# $x2 expected value on a natural scale - i.e. p(Y=y | S = s) This is what we would normally report.
# $vals - the Ys in the rejection region. The proportion of the different Ys (different models) is a Pritchard etal-style, rejection-based
#                                         estimate of the posterior probability of the different models. You might get some improvement by 
#                                          weighting the frequencies with $wt.
# $wt - the Epanechnikov weight. 
# $ss - the sumstats corresponding to the points in $vals. 

  if(missing(gwt))gwt <- rep(T,length(sumstat[,1]))
  
  nss <- length(sumstat[1,])
  # scale everything 
  scaled.sumstat <- sumstat
    
  for(j in 1:nss){
    
    scaled.sumstat[,j] <- normalise(sumstat[,j],sumstat[,j][gwt])
  }
  target.s <- target
  
  for(j in 1:nss){
    target.s[j] <- normalise(target[j],sumstat[,j][gwt])
  }
    
  # calc euclidean distance
  sum1 <- 0
  for(j in 1:nss){
    sum1 <- sum1 + (scaled.sumstat[,j]-target.s[j])^2
  }
  dst <- sqrt(sum1)
  
  # includes the effect of gwt in the tolerance
  dst[!gwt] <- floor(max(dst[gwt])+10)
  
  # wt1 defines the region we're interested in 
  abstol <- quantile(dst,tol)
  # making sure all simulated results are included when tol = 1.
  if(tol == 1) {abstol <- abstol * 1.1}
  wt1 <- dst < abstol
  nit <- sum(wt1)
  
  if(rejmethod){
    l1 <- list(x=x[wt1],wt=0)
  }
  else{
    regwt <- 1-dst[wt1]^2/abstol^2
    
    catx <- as.numeric(names(table(x[wt1])))
    ncat <- length(catx)
    yy <- matrix(nrow=nit,ncol=ncat)
    for(j in 1:ncat){
      yy[,j] <- as.numeric(x[wt1] == catx[j])
    }
    # Used to print generic column names: mu1, mu2, ... for $x2
    # now it prints the "cagegory names" from x, Naoki
    colnames(yy)  <- names(table(x[wt1]))
    
    tr <- list()
    
    for(j in 1:nss){
      tr[[j]] <- scaled.sumstat[wt1,j]
    }
    
    xvar.names <- paste("v",as.character(c(1:nss)),sep="")
    
    names(tr) <- xvar.names
    
    fmla <- as.formula(paste("yy ~ ", paste(xvar.names, collapse= "+")))
    
#        fit1 <- vglm(fmla,data=tr,multinomial) this is the version described in the 
#               manuscript,which did not use the Epanechnikov weights. 
        

    fit1 <- vglm(fmla,data=tr,multinomial,weights=regwt)
        
    target.df <- list()
    for(j in 1:nss){
      target.df[[j]] <- target.s[j]
    }
    names(target.df) <- xvar.names
    
    # predict.vglm renamed to predictvglm in VGAM version 0.8-4
    if (exists('predict.vglm')) {
        prediction1 <- predict.vglm(fit1,target.df,se.fit=T)
        prediction2 <- predict.vglm(fit1,target.df,type="response")
    } else {
        prediction1 <- predictvglm(fit1,target.df,se.fit=T)
        prediction2 <- predictvglm(fit1,target.df,type="response")
    }

    colnames(prediction2) <- colnames(yy) # added by Naoki
    
    l1 <- list(x1=prediction1,x2=prediction2,vals=x[wt1],wt=regwt,ss=sumstat[wt1,])
    
  }
  l1
}

normalise <- function(x,y){
  if(var(y) == 0)return (x - mean(y))
  (x-(mean(y)))/sqrt(var(y))
}

get_indices_of_patterns = function(target_vector, patterns,
        ignore_case=FALSE,
        sort_indices=TRUE) {
    indices = c()
    for (p in patterns) {
        indices = c(indices,
                    grep(p, target_vector, ignore.case=ignore_case))
    }
    if (sort_indices) {
        return(sort(indices))
    }
    return(indices)
}


##############################################################################
## Main CLI

## Option parsing
suppressPackageStartupMessages(library("optparse"))
option_list = list(
    make_option(c("-t", "--tolerance"), type="double", default=1.0,
            dest="tolerance",
            help="Proportion of samples to be accepted (default: 1.0)."),
    make_option("--observed-path", type="character", dest="observed_path",
            help="Path to file with observed summary statistics."),
    make_option("--posterior-path", type="character", dest="posterior_path",
            help="Path to file with unadjusted posterior samples."),
    make_option("--summary-path", type="character",
            dest="summary_path", default="posterior-summary.txt",
            help=paste(
                "Path to posterior summary output file. This file will",
                "contain summary statistics (i.e., mean, median, mode,",
                "CIs) for the regression-adjusted posterior samples of",
                "each continuous parameter specified by",
                "`--continuous-prefixes'. It will also contain the posterior",
                "probabilities for discrete parameters specified by the",
                "'--discrete-prefixes' option. It will also indicate whether",
                "the multinomial logistic regression of each discrete",
                "parameter failed. The unadjusted probabilities are always",
                "reported, and the adjusted values are reported if the",
                "regression worked. This file is formatted as a Python",
                "config file for easy parsing.",
                sep = '\n\t\t')),
    make_option(c("--adjusted-path"), type="character",
            dest="adjusted_path", default="adjusted-posterior-samples.txt",
            help=paste(
                "Path to output file for the regression-adjusted",
                "posterior samples. This will contain the adjusted values",
                "for continuous parameters specified with the",
                "'--continuous-prefixes' option. The default is",
                "'./adjusted-posterior-samples.txt'.",
                sep = '\n\t\t')),
    make_option("--continuous-prefixes", type="character",
            default="PRI.omega,PRI.E.t",
            dest="continuous_prefixes",
            help=paste(
                "The comma-separated prefixes of continuous parameters",
                "to analyze. The default is 'PRI.omega,PRI.E.t'. If you",
                "specify an empty string, no continuous parameters will",
                "will be analyzed.",
                sep = '\n\t\t')),
    make_option("--discrete-prefixes", type="character",
            default="PRI.Psi,PRI.model",
            dest="discrete_prefixes",
            help=paste(
                "The comma-separated prefixes of discrete parameters",
                "to analyze. The default is 'PRI.Psi,PRI.model'. If you",
                "specify an empty string, no discrete parameters will be",
                "analyzed.",
                sep = '\n\t\t')),
    make_option("--stat-prefixes", type="character",
            default="pi,wattTheta,pi.net,tajD.denom",
            dest="stat_prefixes",
            help=paste(
                "The comma-separated prefixes of summary statistics",
                "to use in the analysis. The default is",
                "'pi,wattTheta,pi.net,tajD.denom'.",
                sep = '\n\t\t')),
    make_option(c("-c", "--continuous-indices"), type="character",
            default=NULL,
            dest="continuous_indices",
            help=paste(
                "The comma-separated column indices of continuous parameters",
                "to analyze. The default is 'NULL' (i.e.",
                "'--continuous-prefixes' is used). If you specify this",
                "option, it will override the prefix option",
                sep = '\n\t\t')),
    make_option(c("-d", "--discrete-indices"), type="character",
            default=NULL,
            dest="discrete_indices",
            help=paste(
                "The comma-separated column indices of discrete parameters",
                "to analyze. The default is 'NULL' (i.e.",
                "'--discrete-prefixes' is used). If you specify this",
                "option, it will override the prefix option",
                sep = '\n\t\t')),
    make_option(c("-s", "--stat-indices"), type="character",
            default=NULL,
            dest="stat_indices",
            help=paste(
                "The comma-separated indices of summary statistics",
                "to use in the analysis. The default is 'NULL' (i.e.",
                "'--stat-prefixes' is used). If you specify this",
                "option, it will override the prefix option",
                sep = '\n\t\t')))

options = parse_args(OptionParser(option_list=option_list))

library(VGAM, warn.conflicts=F)
library(KernSmooth)
library(locfit) 

stat_prefixes = strsplit(options$stat_prefixes, ",")[[1]]
continuous_prefixes = strsplit(options$continuous_prefixes, ",")[[1]]
discrete_prefixes = strsplit(options$discrete_prefixes, ",")[[1]]

stat_indices = NULL
discrete_indices = NULL
continuous_indices = NULL
if (!is.null(options$stat_indices)){
    stat_indices = as.numeric(strsplit(options$stat_indices, ",")[[1]])
}
if (!is.null(options$continuous_indices)){
    continuous_indices = as.numeric(strsplit(options$continuous_indices, ",")[[1]])
}
if (!is.null(options$discrete_indices)){
    discrete_indices = as.numeric(strsplit(options$discrete_indices, ",")[[1]])
}

rejmethod = FALSE
if (options$tolerance < 1.0) {
    rejmethod = TRUE
}

res = stdAnalysis(obs.infile = options$observed_path,
		sim.infile = options$posterior_path,
		out.file = options$summary_path,
        adjusted_path = options$adjusted_path,
        tol = options$tolerance,
        rejmethod = rejmethod,
        stat_prefixes = stat_prefixes,
        continuous_prefixes = continuous_prefixes,
        discrete_prefixes = discrete_prefixes,
        stat_indices = stat_indices,
        discrete_indices = discrete_indices,
        continuous_indices = continuous_indices)

