
####################################
## Shell Function
####################################


tjbal <- function(
    formula = NULL,
    data, # data in long form
    Y, # outcome
    D, # treatment
    X = NULL, # time-invariant covariates
    X.avg.time = NULL, # take averages of covariates in a given time period
    index, # unit and time
    Y.match.periods = NULL,
    demean = TRUE, # take out pre-treatment unit mean
    kernel = FALSE, # kernel method
    sigma=NULL,
    maxnumdims = NULL,
    method="ebal",
    whiten=FALSE,
    test = FALSE, ## test different sigmas
    nsigma = 16,
    kbal.step = 1,
    print.baltable = TRUE, # print out table table
    bootstrap = FALSE, ## uncertainty via bootstrap
    conf.lvl = 0.95, ## confidence interval
    nboots = 500, ## number of bootstrap runs
    parallel = TRUE, ## parallel computing
    cores = 4,
    seed = 1234
    ) {
    UseMethod("tjbal")
}   

####################################
## Main Functions
####################################

tjbal.formula <- function(
    formula = NULL,
    data, # data in long form
    X.avg.time = NULL, # take averages of covariates in a given time period
    index, # unit and time
    Y.match.periods = NULL,
    demean = TRUE, # take out pre-treatment unit mean
    kernel = FALSE, # kernel method
    sigma=NULL,
    maxnumdims = NULL,
    method="ebal",
    whiten=FALSE,
    test = FALSE, ## test different sigmas
    nsigma = 16,
    kbal.step = 1,
    print.baltable = TRUE, # print out table table
    bootstrap = FALSE, ## uncertainty via bootstrap
    conf.lvl = 0.95, ## confidence interval
    nboots = 500, ## number of bootstrap runs
    parallel = TRUE, ## parallel computing
    cores = 4,
    seed = 1234
    ) {

    ## parsing
    varnames <- all.vars(formula)
    Yname <- varnames[1]
    Dname <- varnames[2]
    if (length(varnames) > 2) {
        Xname <- varnames[3:length(varnames)]
    } else {
        Xname <- NULL
    }

    namesData <- colnames(data)
    for (i in 1:length(varnames)) {
        if(!varnames[i] %in% namesData) {
            stop(paste0("variable \"", varnames[i],"\" is not in the dataset."))
        }
    }

    ## run the model
    out <- tjbal.default(data = data, Y = Yname,
                          D = Dname, X = Xname,
                          X.avg.time, index, Y.match.periods, demean, kernel, sigma,
                          maxnumdims, method, whiten, test, nsigma, print.baltable, 
                          bootstrap, conf.lvl, nboots, parallel, cores)
    
    out$call <- match.call()
    out$formula <- formula
    ## print(out)
    return(out)

}

tjbal.default <- function(
    data, # data in long form
    Y, # outcome
    D, # treatment
    X = NULL, # time-invariant covariates
    X.avg.time = NULL, # take averages of covariates in a given time period
    index, # unit and time
    Y.match.periods = NULL,
    demean = TRUE, # take out pre-treatment unit mean
    kernel = FALSE, # kernel method
    sigma=NULL,
    maxnumdims = NULL,
    method="ebal",
    whiten=FALSE,
    test = FALSE, ## test different sigmas
    nsigma = 16,
    kbal.step = 1,
    print.baltable = TRUE, # print out table table
    bootstrap = FALSE, ## uncertainty via bootstrap
    conf.lvl = 0.95, ## confidence interval
    nboots = 500, ## number of bootstrap runs
    parallel = TRUE, ## parallel computing
    cores = 4,
    seed = 1234
    ) {

    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------##  

    
    if (is.data.frame(data) == FALSE) {
        warning("Not a data frame.")
        data <- as.data.frame(data)
    }
    data <- droplevels(data)

    ## index
    if (length(index) != 2 | sum(index %in% colnames(data)) != 2) {
        stop("\"index\" option misspecified. Try, for example, index = c(\"unit.id\", \"time\").")
    }


    ##-------------------------------##
    ## Parsing raw data
    ##-------------------------------##  

    ##store variable names
    Yname <- Y
    Dname <- D
    Xname <- X
    
    id <- index[1];
    time <- index[2];
    TT <- length(unique(data[,time]))
    N <- length(unique(data[,id]))
    p <- length(Xname)
    
    ## check balanced panel
    if (var(table(data[,id])) + var(table(data[, time])) > 0) {
        stop("The panel is not balanced.")
    }

    ## time should be numeric
    if (is.numeric(data[,time])==FALSE) {
        stop("The time indicator must be numeric.")
    }
    
    ## check missingness
    if (sum(is.na(data[, Yname])) > 0) {
        stop(paste("Missing values in variable \"", Yname,"\".", sep = ""))
    }
    if (sum(is.na(data[, Dname])) > 0) {
        stop(paste("Missing values in variable \"", Dname,"\".", sep = ""))
    }    
    if (sum(is.na(data[, id])) > 0) {
        stop(paste("Missing values in variable \"", id,"\".", sep = ""))
    }
    if (sum(is.na(data[, time])) > 0) {
        stop(paste("Missing values in variable \"", time,"\".", sep = ""))
    } 

 
    ## sort data
    data <- data[order(data[,id], data[,time]), ]

    ## time and unit
    Ttot <- sort(unique(data[,time]))
    units <- unique(data[,id])
    
    ##treatment indicator
    D<- matrix(data[,Dname],TT,N)

    ## once treated, always treated
    D <- apply(D, 2, function(vec){cumsum(vec)})
    D <- ifelse(D > 0, 1, 0)

    ## treatment
    treat <-ifelse(D[TT,]==1, 1, 0)     # cross-sectional: treated unit
    id.tr <- which(treat == 1)
    id.co <- which(treat == 0)
    Ntr <- length(id.tr)
    Nco <- length(id.co)

    ## check the number of treated units
    if (Ntr <= 5) {
        warning("Too few treated unit(s). Uncertainty estimates not provided.")
        bootstrap <- FALSE
    }

    ## treatment timing
    T0 <- apply((D==0),2,sum) 
    T0.tr <- T0[id.tr]
    T0.min<-min(T0.tr)
    ## same timing: 
    if (Ntr==1) {
        sameT0 <- TRUE
    } else {
        if (var(T0.tr)==0) {
            sameT0 <- TRUE        
        } else {
            sameT0 <- FALSE
        }
    }    
    if (sameT0==TRUE) {
        Tpre <- Ttot[1:unique(T0.tr)]        
    }


    ## outcome variable
    outcome <- matrix(data[,Yname],N, TT, byrow = TRUE)
    colnames(outcome) <- paste0(Yname,Ttot) ## including both pre and post

    ## covariates (allow missing, but non-missing values have to be same for each unit)
    if (p > 0) {
        if (is.null(X.avg.time)==FALSE) {
            if (sameT0 == FALSE) {
                stop("\"X.avg.time\" is only allowed when the treatment starts at the same time.")
            }
            if (is.list(X.avg.time)==TRUE) {
                if (length(X.avg.time)!=p) {
                    stop("Length of \"X.avg.time\" (as a list) must equal the number of covariates.")
                }
                Xvar <- matrix(NA, N, p)
                colnames(Xvar) <- Xname
                for (i in 1:p) {
                    this.period <- X.avg.time[[i]]
                    if (sum(1 - this.period%in%Tpre)>0) {
                        stop("Elements in \"X.avg.time\" must be in the pre-treatment period.")
                    }
                    selected.row <- which(data[,time] %in% this.period)
                    X.pre <- data[selected.row, c(id,Xname[i]),drop = FALSE] 
                    covar.tmp <- ddply(X.pre, .(unit = X.pre[, id]), 
                        numcolwise(mean), na.rm = TRUE)[,-1]
                    if (length(covar.tmp)!=N) {
                        stop(paste0("Missing values in ",Xname[i]," in specified years."))
                    } else{
                        Xvar[,i] <- covar.tmp
                    }
                }
            } else { # not a list, a set of numbers only
                if (sum(1 - X.avg.time%in%Tpre)>0) {
                    stop("\"X.avg.time\" must be in the pre-treatment period.")
                }
                selected.row <- which(data[,time] %in% X.avg.time)
                X.pre <- data[selected.row, Xname, drop = FALSE] 
                Xvar <- ddply(X.pre, .(unit = data[selected.row, id]), 
                    numcolwise(mean), na.rm = TRUE)[,-1]
                ## check missingness
                if (nrow(Xvar)!= N) {
                    stop("Missing values in covariates.")
                }
            }
            ## check missingness again
            for (i in 1:p) {                
                if (sum(is.na(Xvar[, i])) > 0) {
                    stop(paste0("Missing values in variable \"", Xname[i],"\".", sep = ""))
                }
            }
        } else { # no X.avg.time is given
           Xvar <- matrix(NA, N, p);  colnames(Xvar) <- Xname
           for (i in 1:p) {
                if (sum(is.na(data[, Xname[i]])) > 0) {
                    warning(paste0("Missing values in variable \"", Xname[i],"\".", sep = ""))
                }
                ## check variation
                X.tmp <- matrix(data[,Xname[i]], N, TT, byrow = TRUE)
                X.var <- apply(X.tmp,1,var,na.rm = TRUE)
                if (sum(is.na(X.var))>0) {
                    stop(paste0("Variable \"", Xname[i], "\" is completely missing in some unit(s)."))
                } 
                if (sum(X.var)!=0) { # if not time-invariant
                    stop(paste0("\"", Xname[i],"\" is not time-invariant for some unit(s)."))
                }
                ## fill in the matrix
                Xvar[,i] <- apply(X.tmp, 1, mean, na.rm=TRUE) # the first period
            }   
        }
    }

    
    
    ## prepare "wide" form data
    if (p>0) {
        data.wide <- cbind.data.frame(id = units, treat = treat, T0 = T0, outcome, Xvar)
    } else {
        data.wide <- cbind.data.frame(id = units, treat = treat, T0 = T0, outcome)
    } 

    #######################
    ## balancing
    #######################

    if (sameT0 == TRUE) {

        bal.out <- tjbalance(data = data.wide, Y = Yname, D = "treat", X = Xname,
            Y.match.periods = Y.match.periods, Ttot = Ttot, Tpre = Tpre, unit = "id", 
            demean = demean, kernel = kernel, sigma = sigma, maxnumdims = maxnumdims, 
            method=method,whiten = whiten, test = test, nsigma = nsigma, kbal.step = kbal.step,
            bootstrap = bootstrap, nboots = nboots, conf.lvl = conf.lvl, parallel = parallel, cores = cores) 

    } else { ## will not be able to give uncertainty estimates

        bal.out <- tjbalance.mult(data = data.wide, Y = Yname, D = "treat", X = Xname,
            Y.match.periods = Y.match.periods, Ttot = Ttot, unit = "id", 
            demean = demean, kernel = kernel, sigma = sigma, kbal.step = kbal.step,
            maxnumdims = maxnumdims, method=method, whiten = whiten, test = test, nsigma = nsigma) 
    }

    out <- c(list(sameT0 = sameT0), bal.out)
    out$call <- match.call()
    class(out) <- "tjbal"
    return(out)

}



###################################################
#### trajectory balancing with one treatment timing
###################################################


linearweights  <- function(X, D, ebal.tol=1e-4){
    N <- length(D)
    Ntr <- sum(D)
    ebal.out <- ebalance(Treatment = D, X= X, constraint.tolerance=ebal.tol, max.iterations = 5000,print.level=-1)
    w <- rep(1, N)
    w[D==0] <- ebal.out$w * (-1)  ## non-convergence can occur if max.iteration is set too small
    return(w/Ntr) # treated add up to 1; controls add up to -1
}

tjbalance <- function(
    data, ## wide form
    Y,
    D,
    X = NULL,
    Y.match.periods = NULL,
    Ttot,
    Tpre,
    unit,
    demean = FALSE, # take out pre-treatment unit mean
    kernel = FALSE, # kernel method
    sigma=NULL,
    maxnumdims = NULL,
    method="ebal",
    whiten=FALSE,
    test = FALSE, ## test different sigmas
    nsigma = 16,
    kbal.step = 1,
    conf.lvl = 0.95,
    print.baltable = TRUE, # print balance table 
    bootstrap = FALSE, ## uncertainty via bootstrap
    nboots = 500, ## number of bootstrap runs
    parallel = FALSE, ## parallel computing
    cores = 4,
    seed = 1234    
    ){

    
    ## process data
    N <- dim(data)[1]
    Ntr <-  length(which(data[,D]==1))
    Nco <- N - Ntr
    TT <- length(Ttot)
    T0 <- length(Tpre)
    Tpst <- setdiff(Ttot, Tpre)

    ## default: match on all pre-treatment periods
    excludeY <- FALSE
    if (is.null(Y.match.periods)==TRUE) {
        Y.match.periods <- Tpre
    } else {
        if (length(Y.match.periods)==1) {
            if (Y.match.periods == "none") {
                Y.match.periods <- NULL
                demean <- TRUE
                excludeY <- TRUE
            } 
        }        
    }    

    ## outcome variable names (wide form)
    Y.var <- paste0(Y, Ttot)
    
    ## demean
    if (demean == TRUE) { 
        Y.dm.var <- paste0(Y,".dm",Ttot)
        Ypre.mean <- apply(data[, paste0(Y, Tpre), drop = FALSE], 1, mean) # N*1
        outcome.dm <- data[, Y.var, drop = FALSE] - matrix(Ypre.mean, N, TT) # N * TT
        colnames(outcome.dm) <- Y.dm.var
        data <- cbind.data.frame(data, outcome.dm)              
        Y.match <- paste0(Y,".dm",Y.match.periods)
        Y.target <- Y.dm.var
        Y.target.pst <- paste0(Y,".dm",Tpst)
    } else {
        Y.match <- paste0(Y, Y.match.periods)
        Y.target <- Y.var
        Y.target.pst <- paste0(Y, Tpst)
    }
    if (excludeY == TRUE) {
        Y.match <- NULL
    }
    matchvar <- c(Y.match, X)
    
    ## default weights
    w <- rep(NA, N)
    weights.tr <- rep(1/Ntr, Ntr) 
    weights.co <- rep(1/Nco, Nco)
    w[data[,D] == 1] <- weights.tr  
    w[data[,D] == 0] <- weights.co * (-1)  # controls add up to -1;

    ## have something to balance on
    if (is.null(matchvar)==FALSE) { 

        cat("\nSeek balance on:\n")
        cat(paste(matchvar, collapse = ", "),"\n\nOptimization:\n")

        ## tuning parameter
        if (is.null(sigma)){
            if (test == TRUE) {
                ## sigmas for testing
                sigma <- exp(seq(0,log(1500),length.out = nsigma)) 
            } else {
                sigma <- length(matchvar)
                nsigma <- 1
            } 
        } else {
            nsigma <- length(sigma)
        }
        if (kernel == FALSE) {
            sigma <- 1
            nsigma <- 1
        }

        ## balancing
        L1.ratios <- rep(NA, nsigma)
        L1.ratio.best <- 1
        kbal.out.best <- sigma.best <- NULL
        for (i in 1:nsigma) {
            kbal.out <- kbal(X = data[,matchvar],
                D = data[,D], method=method,
                sigma=sigma[i], maxnumdims = maxnumdims,
                linkernel = (1-kernel), incrementby = kbal.step)
            if (kbal.out$dist.record==999) {
                L1.ratio <- 1
            } else {
                L1.before <- kbal.out$L1_orig
                L1.after <- kbal.out$L1_kbal
                L1.ratio <- L1.after/L1.before            
            }
            L1.ratios[i] <- L1.ratio        
            if (L1.ratio - L1.ratio.best < 1e-5) {
                L1.ratio.best <- L1.ratio
                kbal.out.best <- kbal.out
                sigma.best <- sigma[i]                      
            } 
        }
        if (nsigma > 1) {
            test.out <- cbind(sigma, L1.ratios)
            colnames(test.out) <- c("sigma","L1.ratio")
            print(round(test.out, 2)) 
        }

        ## if success
        if (1 - L1.ratio.best >1e-5) {
            success <- 1
            ndims <- kbal.out.best$numdims
            cat(paste0("sigma* = ", round(sigma.best,1),"; L1.ratio% = ", 
                sprintf("%.3f",L1.ratio.best*100),"; num.dims = ",ndims,"\n"))

            ## weights
            weights.tr <- rep(1/Ntr, Ntr) # treated add up to 1; 
            weights.co <- kbal.out.best$w[data[,D] == 0]/Nco # controls add up to 1;
            w[data[,D] == 1] <- weights.tr  
            w[data[,D] == 0] <- weights.co * (-1)  # controls add up to -1;

        } else {
            success <- 0
            cat("\nSolution not found. Equal weights are being used.\n")            
        }

    }
    
    # ATT
    names.co <- data[data[,D]==0, unit]
    names(weights.co) <- names.co
    att <- apply(data[, Y.target] * w, 2, sum)
    att.avg <- mean(att[Y.target.pst])


    # ## treated and control data
    id.co <- which(data[,D]==0)
    id.tr <- which(data[,D]==1)
    Y.tr = data[id.tr, Y.var, drop = FALSE]
    Y.co = data[id.co, Y.var, drop = FALSE]

    Y.tr.bar <- apply(data[id.tr, Y.var], 2, mean, na.rm=TRUE)
    Y.co.bar <- apply(data[id.co, Y.var], 2, mean, na.rm=TRUE)
    Y.ct.bar <- Y.tr.bar - att
    Y.bar <- cbind(Y.tr.bar,Y.ct.bar, Y.co.bar)

    ## faster routine with perfect mean balancing
    if (is.null(matchvar)==FALSE) {
        if (kernel == FALSE & abs(L1.ratio)<1e-4) { 
            fastlinear <- TRUE
            ## deal with colinearity 
            if (demean == TRUE && length(Y.match.periods)==T0) {
                matchvar <- c(Y.match[-1],X)            
            }
            # standardize
            X.tmp <- data[,matchvar, drop = FALSE]
            for (i in 1:ncol(X.tmp)) {
                X.tmp[,i] <- X.tmp[,i]/sd(X.tmp[,i],na.rm=TRUE) 
            }
        } else {
            fastlinear <- FALSE
            K <- kbal.out.best$K        
        }
    }    

    ## bootstrap
    if (bootstrap == TRUE) {

        ## seed
        if (is.null(seed) == FALSE) {
            if (is.numeric(seed) == FALSE) {
                stop("seed should be a number.")
            }
        } else {
            seed <- 1234
        }
        set.seed(seed)            

        one.boot <- function() {
            sample.id <- c(sample(id.tr,Ntr, replace = TRUE),
                sample(id.co, Nco, replace = TRUE))
            ## weights: treated add up to 1; controls add up to -1; sum is zero
            if (is.null(matchvar) == TRUE) {
                w.boot <- rep(1/Ntr, N)
                w.boot[data[sample.id,D] == 0] <- rep(-1/Nco, Nco)
            } else {
                if (fastlinear == TRUE) {
                    w.boot <- linearweights(X = X.tmp[sample.id,], D = data[sample.id, D])
                } else {                  
                    data.tmp <- data[sample.id, ]
                    K.boot <- K[sample.id, sample.id]                                
                    kbal.boot <- kbal(X = data.tmp[, matchvar], D = data.tmp[, D], K = K.boot,
                        method=method, linkernel = (1-kernel), incrementby = kbal.step)
                    w.boot <-  rep(1/Ntr, N)
                    w.boot[data.tmp[,D] == 0] <- kbal.boot$w[data.tmp[,D] == 0]/Nco*(-1)  # controls add up to -1;   
                }
            }
            
            ## ATT
            att <- apply(data[sample.id, Y.target] * w.boot, 2, sum)
            att.avg <- mean(att[Y.target.pst])
            out <- list(att = att, att.avg = att.avg)
            return(out)            
        }

        ## Storing bootstrapped estimates
        att.boot<-matrix(0,TT,nboots)
        att.avg.boot<-matrix(0,nboots,1)  

        ## computing
        cat("\nBootstrapping...")
        if (parallel == TRUE) {
            ## prepare
            if (is.null(cores) == TRUE) {
                cores <- detectCores()
            }
            para.clusters <- makeCluster(cores)
            registerDoParallel(para.clusters)
            ## start    
            cat("Parallel computing...") 
            boot.out <- foreach(j=1:nboots, 
                .inorder = FALSE,
                .export = c("ebalance","linearweights"),
                .packages = c("KBAL")
                ) %dopar% {
                return(one.boot())
            }
            stopCluster(para.clusters)
            ## save results
            for (j in 1:nboots) { 
                att.boot[,j]<-boot.out[[j]]$att
                att.avg.boot[j,]<-boot.out[[j]]$att.avg                  
            } 
        } else { ## single core
            for (j in 1:nboots) { 
                boot <- one.boot() 
                att.boot[,j]<-boot$att
                att.avg.boot[j,]<-boot$att.avg                
                ## report progress
                if (kernel == FALSE) {
                    if (j%%50==0)  {cat(".")}
                } else {
                    cat(j,"\n")
                } 
            }  
        }
        ## end of bootstrapping
        cat("\n")

        ####################################
        ## Variance and CIs
        ####################################

        ## function to get two-sided p-values
        get.pvalue <- function(vec) {
            if (NaN%in%vec|NA%in%vec) {
                nan.pos <- is.nan(vec)
                na.pos <- is.na(vec)
                pos <- c(which(nan.pos),which(na.pos))
                vec.a <- vec[-pos]
                a <- sum(vec.a >= 0)/(length(vec)-sum(nan.pos|na.pos)) * 2
                b <- sum(vec.a <= 0)/(length(vec)-sum(nan.pos|na.pos)) * 2  
            } else {
                a <- sum(vec >= 0)/length(vec) * 2
                b <- sum(vec <= 0)/length(vec) * 2  
            }
            return(min(as.numeric(min(a, b)),1))
        }

        ## For ATT estimates
        conf.lvl.lb <- (1 - conf.lvl)/2
        conf.lvl.ub <- conf.lvl.lb + conf.lvl

        CI.att <- t(apply(att.boot, 1, function(vec) 
            quantile(vec,c(conf.lvl.lb, conf.lvl.ub), na.rm=TRUE)))
        se.att <- apply(att.boot, 1, function(vec) sd(vec, na.rm=TRUE))
        pvalue.att <- apply(att.boot, 1, get.pvalue)    
        pvalue.att[which(abs(att)<1e-5 & abs(se.att)<1e-5)] <- NA
        est.att <- cbind(att, se.att, CI.att, pvalue.att, ntreated = rep(Ntr,TT))
        est.att[abs(est.att)<1e-5] <- 0
        colnames(est.att) <- c("ATT", "S.E.", "CI.lower", "CI.upper",
         "p.value", "n.Treated")
        rownames(est.att) <- Ttot    

        ## average (over time) ATT
        CI.avg <- quantile(att.avg.boot, c(conf.lvl.lb, conf.lvl.ub), na.rm=TRUE)
        se.avg <- sd(att.avg.boot, na.rm=TRUE)
        pvalue.avg <- get.pvalue(att.avg.boot)
        est.avg <- t(as.matrix(c(att.avg, se.avg, CI.avg, pvalue.avg)))
        colnames(est.avg) <- c("ATT.avg", "S.E.", "CI.lower", "CI.upper", "p.value")

        ## storage
        out.inference <- list(
            est.att = est.att, 
            est.avg = est.avg, 
            att.boot = att.boot, 
            att.avg.boot = att.avg.boot
            )      
    }    

    ## balance table
    if (is.null(matchvar)==FALSE) {
        if (success == 1) {
            weighted.sd <- function(vec, w) {sqrt(sum(w * (vec - weighted.mean(vec,w))^2))}
            if (Ntr>1) {
                # treated
                mean.tr <- apply(data[id.tr, matchvar, drop = FALSE], 2, mean) 
                sd.tr <- apply(data[id.tr, matchvar, drop = FALSE], 2, sd)
                # control
                mean.co.pre <- apply(data[id.co, matchvar, drop = FALSE], 2, mean) 
                sd.co.pre <- apply(data[id.co, matchvar, drop = FALSE], 2, sd)
                # weighted control 
                mean.co.pst <- apply(data[id.co, matchvar, drop = FALSE], 2, weighted.mean, weights.co) 
                sd.co.pst <- apply(data[id.co, matchvar, drop = FALSE], 2, weighted.sd, weights.co) 
                ## standardized difference
                # diff.pre <- (mean.tr - mean.co.pre)/sqrt(sd.tr^2 + sd.co.pre^2)
                # diff.pst <- (mean.tr - mean.co.pst)/sqrt(sd.tr^2 + sd.co.pst^2)
                # normalize by SD of the treated 
                diff.pre <- (mean.tr - mean.co.pre)/sd.tr
                diff.pst <- (mean.tr - mean.co.pst)/sd.tr                
                bal.table <- cbind.data.frame(mean.tr, mean.co.pre, mean.co.pst, sd.tr, sd.co.pre,  sd.co.pst, diff.pre, diff.pst)
            } else {
                # treated
                mean.tr <- apply(data[id.tr, matchvar, drop = FALSE], 2, mean) 
                # control
                mean.co.pre <- apply(data[id.co, matchvar, drop = FALSE], 2, mean)
                # weighted control 
                mean.co.pst <- apply(data[id.co, matchvar, drop = FALSE], 2, weighted.mean, weights.co) 
                # difference in means
                diff.pre <- (mean.tr - mean.co.pre)/abs(mean.tr)
                diff.pst <- (mean.tr - mean.co.pst)/abs(mean.tr)
                bal.table <- cbind.data.frame(mean.tr, mean.co.pre, mean.co.pst, diff.pre, diff.pst)
            }
            if (print.baltable==TRUE) {
                cat("\nBalance Table\n")
                print(round(bal.table, 4))   
            }
        }
    }

    out <- list(data.wide = data,
            id.tr = id.tr,
            id.co = id.co,
            Y.tr = Y.tr,
            Y.co = Y.co,
            weights.co = weights.co, # Nco * 1 vector
            names.co = names.co,
            Ttot = Ttot,
            Tpre = Tpre,
            Tpst = Tpst,
            T0 = T0,
            N = N,
            Ntr = Ntr,
            Nco = Nco,                
            matchvar = matchvar,
            Y.bar = Y.bar, 
            att = att,
            att.avg = att.avg)
    
    if (is.null(matchvar)==FALSE) {

        out <- c(out,
            list(success = success,
            L1.before = kbal.out.best$L1_orig,
            L1.after = kbal.out.best$L1_kbal,
            L1.ratio = L1.ratio.best,
            min90 = kbal.out.best$min90,
            ndims = kbal.out.best$numdims,
            sigma.best = sigma.best))

        if (success==1) {
            out <- c(out, list(bal.table = bal.table))
        }

        if (nsigma>1) {
            out <- c(out, list(test.out = test.out))
        }
    } 
    
    if (bootstrap == TRUE) {
        out <- c(out, out.inference)
    }       
    return(out)
}


########################################################
#### trajectory balancing with multiple treatment timing
########################################################

tjbalance.mult <- function(
    data, ## data in wide form
    Y,
    D,
    X,
    Y.match.periods = NULL,
    Ttot,
    unit,
    demean = FALSE, # take out pre-treatment unit mean
    kernel = FALSE, # kernel method
    maxnumdims = NULL,
    method="ebal",
    whiten=FALSE,
    test = TRUE, ## test different sigmas
    sigma = NULL, ## tuning parameters
    nsigma = 16,
    kbal.step = 1  
    ) { 

    TT <- length(Ttot)
    id.tr <- which(data$tr == 1)
    id.co <- which(data$tr == 0)

    ## parse multiple treatment timing
    T0.all <- data$T0
    names(T0.all) <- data$id
    T0.tr <- T0.all[id.tr] # T0 for all treated units
    tb.T0 <- table(T0.tr)
    T0.unique <- as.numeric(names(tb.T0))
    T0.count <- as.numeric(tb.T0)
    T0.max <- max(T0.tr)
    T0.min <- min(T0.tr)
    T0.names <- paste0("T0=",T0.unique)

    ## recenter based on treatment timing
    time.adj <- c(-(T0.max -1) : (TT-T0.min))
    TT.adj <- length(time.adj)
        ## storage
    weights.co <- matrix(NA, length(id.co), length(T0.unique))
    rownames(weights.co) <- data$id[id.co]
    colnames(weights.co) <- T0.names
    sub.att <- sub.Ytr.avg <- sub.Yct.avg <- matrix(NA, TT, length(T0.unique))
    rownames(sub.att) <- rownames(sub.Ytr.avg) <- rownames(sub.Yct.avg) <- Ttot
    colnames(sub.att) <- colnames(sub.Ytr.avg) <- colnames(sub.Yct.avg) <- T0.names
    sub.ntr <- sub.att.adj <- matrix(0, TT.adj, length(T0.unique)) ## number of treated units for each subgroup
    rownames(sub.ntr) <- rownames(sub.att.adj) <- time.adj
    colnames(sub.ntr) <- colnames(sub.att.adj) <- T0.names
    sigmas.out <- L1 <- success <- rep(0, length(T0.unique)) 

        ## loop over different timings
    for (i in 1:length(T0.unique)) {

        T0 <- T0.unique[i]
        this.tr <- which(T0.all == T0)    

        Tpre <- Ttot[1:T0]
        Tpst <- Ttot[(T0+1):TT]
        T.mid <- floor(median(Tpre))

        ## remove other treated 
        data.tjbal <- data[c(id.co, this.tr),]
        
        ## tuning parameter
        if (is.null(sigma)==TRUE) {
            sigma <- NULL
        } else {
            sigma <- sigma[i]
        }

            ## trajectory balancing
        cat("Subgroup: T0 =",T0,"\n")
        tmp <- capture.output(
            tjbal.out <- tjbalance(data = data.tjbal, Y = Y, D = D, X = X,
               Y.match.periods = Tpre, Ttot = Ttot, Tpre = Tpre, method = method,
               unit = "id", demean = TRUE, sigma = sigma, nsigma = 16, kernel == kernel,
               maxnumdims = (length(id.co)-1))
            , file = NULL)        

        ## save
        L1[i] <- tjbal.out$L1.ratio
        sigmas.out[i] <- tjbal.out$sigma.best
        success[i] <- tjbal.out$success

        weights.co[,i] <- tjbal.out$weights.co
        sub.Ytr.avg[,i] <- tjbal.out$Y.bar[,1]
        sub.Yct.avg[,i] <- tjbal.out$Y.bar[,2]
        att <- tjbal.out$att
        sub.att[,i] <- att

        ## save ATT (realigned based on T0)        
        fill.start <- T0.max-T0+1
        fill.end <- fill.start + length(att) -1 
        sub.ntr[fill.start:fill.end, i] <- T0.count[i]   
        sub.att.adj[fill.start:fill.end, i] <- att
    }
    ntreated <- rowSums(sub.ntr) # how the number of units changes over adjusted time
    att <- rowSums(sub.att.adj * sub.ntr)/ntreated
    names(ntreated) <- names(att) <- time.adj
    sub.att.adj[sub.ntr==0] <- NA

    ## save results
    out <- list(
        id.tr = id.tr,
        id.co = id.co,
        Y.co = Y.co,
        Y.tr = Y.tr,
        Ttot = Ttot,
        Tpre = Tpre,
        Tpst = Tpst,
        N = N,
        Ntr = Ntr,
        Nco = Nco,            
        T0 = T0.unique, 
        T0.all = T0.all,
        T0.tr = T0.tr,
        weights.co = weights.co,
        names.co = colnames(weights.co),
        sub.Ytr.avg = sub.Ytr.avg,
        sub.Yct.avg = sub.Yct.avg,
        sub.att = sub.att,
        sub.ntr = sub.ntr,
        sub.att.adj = sub.att.adj,
        ntreated = ntreated,
        att = att,
        success = success,
        sigma.best = sigmas.out,           
        L1.ratio = L1
        )     
    
    return(out)
}

#######################################################
## METHODS
#######################################################


##########
## Print
##########
## a gsynth object
print.tjbal <- function(x,  
                         ...) {
    
    cat("Call:\n")
    print(x$call, digits = 4)
    
    if (is.null(x$est.avg) == TRUE) { # no uncertainties
        cat("\nAverage Treatment Effect on the Treated:\n")
        print(x$att.avg, digits = 4)
        cat("\n   ~ by Period (including Pre-treatment Periods):\n")
        print(x$att, digits = 4)
        cat("\nUncertainty estimates not available.\n")
    } else {
        cat("\nAverage Treatment Effect on the Treated:\n")
        print(x$est.avg, digits = 4)
        cat("\n   ~ by Period (including Pre-treatment Periods):\n")
        print(x$est.att, digits = 4)        
    }
}


##########
## Plot
##########

plot.tjbal <- function(x,  
    type = "gap", 
    xlim = NULL, 
    ylim = NULL,
    xlab = NULL, 
    ylab = NULL,
    legendOff = FALSE,
    main = NULL,
    raw = "none",
    wmin = -20, ## minimal log weights
    stat = "mean",
    trim = TRUE, ## trim control group in ct plot
    trim.wtot = 0.9, ## show controls whose weights sum up to a number
    theme.bw = TRUE, ## black/white or gray theme
    axis.adjust = FALSE,
    ...){


    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------##  

    outcome <- NULL
    ATT <- NULL
    CI.lower <- NULL
    CI.upper <- NULL
    co.lower <- NULL
    co.upper <- NULL
    tr.lower <- NULL
    tr.upper <- NULL
    group <- NULL
    L1 <- NULL
    out <- NULL

    if (class(x)!="tjbal") {
        stop("Not a \"tjbal\" object.")
    }
    if (!type %in% c("gap","counterfactual","ct","balance","weights")) {
        stop("\"type\" option misspecified.")        
    }
    if (type == "ct") {
        type <- "counterfactual"
    }
    if (is.null(xlim)==FALSE) {
        if (is.numeric(xlim)==FALSE) {
            stop("Some element in \"xlim\" is not numeric.")
        } else {
            if (length(xlim)!=2) {
                stop("xlim must be of length 2.")
            }
        }
    }
    if (is.null(ylim)==FALSE) {
        ## if (type!="missing") {
            if (is.numeric(ylim)==FALSE) {
                stop("Some element in \"ylim\" is not numeric.")
            } else {
                if (length(ylim)!=2) {
                    stop("ylim must be of length 2.")
                }
            }
    }
    if (is.null(xlab)==FALSE) {
        if (is.character(xlab) == FALSE) {
            stop("\"xlab\" is not a string.")
        } else {
            xlab <- xlab[1]
        }   
    }
    if (is.null(ylab)==FALSE) {
        if (is.character(ylab) == FALSE) {
            stop("\"ylab\" is not a string.")
        } else {
            ylab <- ylab[1]
        }   
    }
    if (is.logical(legendOff) == FALSE & is.numeric(legendOff)==FALSE) {
        stop("\"legendOff\" is not a logical flag.")
    }
    if (type == "counterfactual") {
        if (! raw %in% c("none","band","all")) {
            cat("\"raw\" option misspecified. Reset to \"none\".")
            raw <- "none" 
        }        
    }
    if (is.null(trim.wtot)==FALSE) {
        if (is.numeric(trim.wtot)==FALSE) {
            stop("\"trim.w\" is not numeric.")
        } else {
            if (trim.wtot <0 | trim.wtot > 1) {
                stop("\"trim.w\" must be between 0 and 1.")
            }
        }
    }
    if (is.null(main)==FALSE) {
        if (is.character(main) == FALSE) {
            stop("\"main\" is not a string.")
        } else {
            main <- main[1]
        }   
    }
    if (axis.adjust==TRUE) {
        angle <- 45
        x.v <- 1
        x.h <- 1
    } else {
        angle <- 0
        x.v <- 0
        x.h <- 0
    }

    ## color of axes
    if (theme.bw == TRUE) {
      line.color <- "#AAAAAA70"
    } else {
      line.color <- "white"
    }
    
    ##-------------------------------##
    ## Plotting
    ##-------------------------------##  

    Y.tr <- x$Y.tr
    Y.co <- x$Y.co
    tb <- x$est.att
    Yb <- x$Y.bar[,1:2] ## treated average and counterfactual average
    #tr <- x$tr
    #pre <- x$pre
    #post <- x$post
    TT <- length(x$Ttot)
    T0 <- x$T0 ## notice
    Ntr <- x$Ntr
    Nco <- x$Nco
    N <- x$N 
    time <- x$Ttot
    w.co <- x$weights.co
    id.co <- x$id.co
    id.tr <- x$id.tr
    
   
    ## parameters
    line.width <- c(1.2,0.5)
  
    ## type of plots
    if ((type == "gap" && x$sameT0 == TRUE) | type == "counterfactual" ) {
        if (is.numeric(time[1])==FALSE) {
            time <- 1:TT
        }
        time.bf <- time[unique(T0)]

        ## periods to show
        if (length(xlim) != 0) {
            show <- which(time>=xlim[1]& time<=xlim[2])
        } else {
            show <- 1:length(time)
        }
        nT <- length(show)
        time.label <- time[show]     
    }

    if (type == "gap" && x$sameT0 == FALSE)  { ## variable treatment timing
        time <- c(1:TT) - min(T0)
        time.bf <- 0 ## before treatment
        if (length(xlim) != 0) {
            show <- which(time>=xlim[1]& time<=xlim[2])     
        } else {
            show <- 1:length(time)    
        }
        
    }

    if (type == "gap" | type == "counterfactual") {
        nT <- length(show)
        time.label <- time[show]
        T.b <- 1:length(show)
    }


    ## legend on/off
    if (legendOff == TRUE) {
        legend.pos <- "none"
    } else {
        legend.pos <- "bottom"
    }

    ############  START  ###############
    
    if (type == "gap") { 
        
        if (x$Ntr>1) {
            maintext <- "Average Treatment Effect on the Treated"
        } else {
            maintext <- "Treatment Effect on the Treated"
        }
        ## axes labels
        if (is.null(xlab) == TRUE) {
            if (x$sameT0 == TRUE) {
                xlab <- x$index[2]
            } else {
                xlab <- paste("Time relative to Treatment")
            }
        } else if (xlab == "") {
            xlab <- NULL
        }
        if (is.null(ylab) == TRUE) {
            ylab <- "Estimate"
        } else if (ylab == "") {
            ylab <- NULL
        }

        ## construct data for plotting
        if (is.null(x$est.att)==TRUE) { 
            cat("Uncertainty estimates not available.\n")
            data <- cbind.data.frame(time, ATT = x$att)[show,]             
        } else {
            data <- cbind.data.frame(time, tb)[show,]
        }
                         
        ## plotting
        p <- ggplot(data)
        if (theme.bw == TRUE) {
          p <- p + theme_bw()
        }
        p <- p +
        geom_vline(xintercept = time.bf, colour=line.color,size = 2) +
        geom_hline(yintercept = 0, colour=line.color,size = 2) +
                ## annotate("rect", xmin= time.bf, xmax= Inf,
                ##          ymin=-Inf, ymax=Inf, alpha = .1,
                ##          fill = "yellow") +
        xlab(xlab) +  ylab(ylab) +
        theme(legend.position = legend.pos,
          plot.title = element_text(size=20,
            hjust = 0.5,
            face="bold",
            margin = margin(10, 0, 10, 0)))


        ## point estimates
        p <- p + geom_line(aes(time, ATT), size = 1.2)

        ## confidence intervals
        if (is.null(x$est.att)==FALSE) {
            p <- p + geom_ribbon(aes(x = time, ymin=CI.lower, ymax=CI.upper),alpha=0.2)
        }

        
    } else if (type=="counterfactual") { 

        if (Ntr == 1| x$sameT0==TRUE) { # same/single treatment timing

            if (trim == TRUE & raw %in% c("band","all")) {
                Nco <-sum(1 - (cumsum(sort(w.co,decreasing = TRUE))>trim.wtot)) + 1  # how many control units left
                trim.id<- order(w.co, decreasing = TRUE)[1:Nco]
                Y.co <- Y.co[trim.id, ,drop = FALSE]
                id.co <- id.co[trim.id]
                co.label <- paste0("Heavily Weighted Controls (",floor(trim.wtot*100),"% weights)")
            } else {
                co.label <- "Controls"
            }

            ## axes labels
            if (is.null(xlab)==TRUE) {
                xlab <- x$index[2]
            } else if (xlab == "") {
                xlab <- NULL
            }
            if (is.null(ylab)==TRUE) {
                ylab <- x$Yname
            } else if (ylab == "") {
                ylab <- NULL
            }

            if (Ntr == 1) {
                maintext <- "Treated and Counterfactual Trajectories"
            } else {
                maintext <- "Treated and Counterfactual Averages"
            }

            if (raw == "none") {
                data <- cbind.data.frame("time" = rep(time[show],2),
                 "outcome" = c(Yb[show,1],
                   Yb[show,2]),
                 "type" = c(rep("tr",nT),
                    rep("co",nT))) 
                        ## theme 
                p <- ggplot(data)
                if (theme.bw == TRUE) {
                  p <- p + theme_bw()
                }
                p <- p + xlab(xlab) +  ylab(ylab) +
                geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                ## shade in the post-treatment period
                #annotate("rect", xmin= time.bf, xmax= Inf,
                # ymin=-Inf, ymax=Inf, alpha = .3) +
                theme(legend.position = legend.pos,
                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                  plot.title = element_text(size=20,
                    hjust = 0.5,
                    face="bold",
                    margin = margin(10, 0, 10, 0)))
                        ## main
                p <- p + geom_line(aes(time, outcome,
                   colour = type,
                   size = type,
                   linetype = type))

                ## legend
                set.limits = c("tr","co")
                if (Ntr == 1) {
                    set.labels = c("Treated", "Estimated Y(0)")
                } else {
                    set.labels = c("Treated Average", "Estimated Y(0) Average")
                }
                set.colors = c("black","steelblue")
                set.linetypes = c("solid","longdash")
                set.linewidth = rep(line.width[1],2)
                p <- p + scale_colour_manual(limits = set.limits,
                 labels = set.labels,
                 values =set.colors) +
                scale_linetype_manual(limits = set.limits,
                  labels = set.labels,
                  values = set.linetypes) +
                scale_size_manual(limits = set.limits,
                  labels = set.labels,
                  values = set.linewidth) +
                guides(linetype = guide_legend(title=NULL, ncol=2),
                   colour = guide_legend(title=NULL, ncol=2),
                   size = guide_legend(title=NULL, ncol=2)) 

                if (!is.numeric(time.label)) {
                    p <- p + 
                    scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                }

            } else if  (raw == "band") {

                if (trim == FALSE) {
                    Y.tr.band <- t(apply(Y.tr, 2, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                    Y.co.band <- t(apply(Y.co, 2, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                    tr.band.label <- "Treated 5-95% Quantiles"
                    co.band.label <- "Controls 5-95% Quantiles"
                } else {
                    Y.tr.band <- t(apply(Y.tr, 2, quantile, prob=c(0,1),na.rm=TRUE))
                    Y.co.band <- t(apply(Y.co, 2, quantile, prob=c(0,1),na.rm=TRUE))
                    tr.band.label <- "Treated"
                    co.band.label <- paste0("Heavily Weighted Controls (",floor(trim.wtot*100),"% weights)")
                }

                
                data <- cbind.data.frame("time" = rep(time[show],2),
                   "outcome" = c(Yb[show,1], Yb[show,2]),"type" = c(rep("tr",nT),rep("co",nT)))

                if (Ntr == 1) {
                    data.band <- cbind.data.frame(time, Y.co.band)[show,]
                    colnames(data.band) <- c("time","co.lower","co.upper")
                } else {
                    data.band <- cbind.data.frame(time, Y.tr.band, Y.co.band)[show,]
                    colnames(data.band) <- c("time","tr.lower","tr.upper","co.lower","co.upper")
                }
                
                ## theme 
                p <- ggplot(data)
                if (theme.bw == TRUE) {
                  p <- p + theme_bw()
                }
                p <- p + xlab(xlab) +  ylab(ylab) +
                geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                #annotate("rect", xmin= time.bf, xmax= Inf,
                # ymin=-Inf, ymax=Inf, alpha = .3) +
                theme(legend.position = legend.pos,
                  axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                  plot.title = element_text(size=20,
                    hjust = 0.5,
                    face="bold",
                    margin = margin(10, 0, 10, 0)))
                        ## main
                p <- p + geom_line(aes(time, outcome,
                   colour = type,
                   size = type,
                   linetype = type))
                ## band
                if (Ntr == 1) {
                    p <- p + geom_ribbon(data = data.band,
                        aes(ymin = co.lower, ymax = co.upper, x=time),
                        alpha = 0.15, fill = "steelblue")
                    set.limits = c("tr","co","co.band")
                    set.labels = c("Treated", "Estimated Y(0)",co.band.label)
                    set.colors = c("black","steelblue","#4682B480")
                    set.linetypes = c("solid","longdash","solid")
                    set.linewidth = c(rep(line.width[1],2),4)                    
                } else {
                    p <- p + geom_ribbon(data = data.band,
                        aes(ymin = co.lower, ymax = co.upper, x=time),
                        alpha = 0.15, fill = "steelblue") +
                            geom_ribbon(data = data.band,
                               aes(ymin = tr.lower, ymax = tr.upper, x=time),
                               alpha = 0.15, fill = "black")
                    set.limits = c("tr","co","tr.band","co.band")
                    set.labels = c("Treated Average", "Estimated Y(0) Average",
                     tr.band.label, co.band.label)
                    set.colors = c("black","steelblue","#77777750","#4682B480")
                    set.linetypes = c("solid","longdash","solid","solid")
                    set.linewidth = c(rep(line.width[1],2),4,4)        
                }                

                p <- p + scale_colour_manual(limits = set.limits,
                 labels = set.labels,
                 values =set.colors) +
                scale_linetype_manual(limits = set.limits,
                  labels = set.labels,
                  values = set.linetypes) +
                scale_size_manual(limits = set.limits,
                  labels = set.labels,
                  values = set.linewidth) +
                guides(linetype = guide_legend(title=NULL, ncol=2),
                   colour = guide_legend(title=NULL, ncol=2),
                   size = guide_legend(title=NULL, ncol=2)) 

                if (!is.numeric(time.label)) {
                    p <- p + 
                    scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                }

            } else if (raw == "all") { ## plot all the raw data

                if (Ntr == 1) {
                    data <- cbind.data.frame("time" = rep(time[show],(2 + Nco)),
                     "outcome" = c(Yb[show,1],Yb[show,2], c(t(Y.co[,show]))),
                     "type" = c(rep("tr",nT), rep("co",nT), rep("raw.co",(Nco * nT))),
                     "id" = c(rep("tr",nT), rep("co",nT), rep(id.co, each = nT))) 
                } else {
                    data <- cbind.data.frame("time" = rep(time[show],(2 + N)),
                     "outcome" = c(Yb[show,1],Yb[show,2],c(t(Y.tr[,show])),c(t(Y.co[,show]))),
                     "type" = c(rep("tr",nT), rep("co",nT), rep("raw.tr",(Ntr * nT)), rep("raw.co",(Nco * nT))),
                     "id" = c(rep("tr",nT), rep("co",nT), rep(c(id.tr,id.co), each = nT))) 
                }

                
                ## theme
                p <- ggplot(data)
                if (theme.bw == TRUE) {
                  p <- p + theme_bw()
                }
                p <- p + xlab(xlab) +  ylab(ylab) +
                        geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                        #annotate("rect", xmin= time.bf, xmax= Inf,
                        #   ymin=-Inf, ymax=Inf, alpha = .3) +
                        theme(legend.position = legend.pos,
                          axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                          plot.title = element_text(size=20,
                            hjust = 0.5,
                            face="bold",
                            margin = margin(10, 0, 10, 0))) 
                
                ## main
                p <- p + geom_line(aes(time, outcome, colour = type,size = type,
                   linetype = type, group = id))

                ## legend
                if (Ntr == 1) {
                    set.limits = c("tr","co","raw.co")
                    set.labels = c("Treated","Estimated Treated Y(0)",co.label)                
                    set.colors = c("black","steelblue","#4682B430")
                    set.linetypes = c("solid","longdash","solid")
                    set.linewidth = line.width[c(1,1,2)]
                } else {
                    set.limits = c("tr","co","raw.tr","raw.co")
                    set.labels = c("Treated Average",
                     "Estimated Y(0) Average",
                     "Treated", co.label)
                    set.colors = c("black","steelblue","#77777750","#4682B430")
                    set.linetypes = c("solid","longdash","solid","solid")
                    set.linewidth = rep(line.width,each=2)
                }
                
                p <- p + scale_colour_manual(limits = set.limits,
                   labels = set.labels,
                   values =set.colors) +
                scale_linetype_manual(limits = set.limits,
                  labels = set.labels,
                  values = set.linetypes) +
                scale_size_manual(limits = set.limits,
                  labels = set.labels,
                  values = set.linewidth)
                if (Ntr ==1) {
                    p <- p + 
                    guides(linetype = guide_legend(title=NULL, ncol=3),
                       colour = guide_legend(title=NULL, ncol=3),
                       size = guide_legend(title=NULL, ncol=3)) 
                } else {
                    p <- p + 
                    guides(linetype = guide_legend(title=NULL, ncol=2),
                       colour = guide_legend(title=NULL, ncol=2),
                       size = guide_legend(title=NULL, ncol=2)) 
                }                

                if (!is.numeric(time.label)) {
                    p <- p + 
                    scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
                }
            }

            # key width in legend
            p <- p + theme(legend.key.width = unit(2.5,"line"))   
            
        } else { ## different treatment timing
            maintext <- "Treated and Counterfactual Averages"

            ## axes labels
            if (is.null(xlab)==TRUE) {
                xlab <- paste("Time relative to Treatment")
            } else if (xlab == "") {
                xlab <- NULL
            }
            if (is.null(ylab)==TRUE) {
                ylab <- x$Yname
            } else if (ylab == "") {
                ylab <- NULL
            }
            
            xx <- ct.adjsut(x$Y.tr, x$Y.ct, x$T0)

            time <- xx$timeline
            Yb <- xx$Yb
            Y.tr.aug <- xx$Y.tr.aug
            ## Y.ct.aug <- xx$Y.ct.aug
            time.bf <- 0 ## before treatment

            if (!is.null(xlim)) {
                show <- which(time>=xlim[1]& time<=xlim[2])
            } else {
                show <- 1:length(time)
            }
            nT <- length(show)

            if (raw == "none") {
                data <- cbind.data.frame("time" = rep(time[show],2),
                                         "outcome" = c(Yb[show,1],
                                                       Yb[show,2]),
                                         "type" = c(rep("tr",nT),
                                                    rep("co",nT))) 
                ## theme 
                p <- ggplot(data) 
                if (theme.bw == TRUE) {
                  p <- p + theme_bw()
                }
                p <- p + xlab(xlab) +  ylab(ylab) +
                    geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                    annotate("rect", xmin= time.bf, xmax= Inf,
                                ymin=-Inf, ymax=Inf, alpha = .3) +
                    theme(legend.position = legend.pos,
                          axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                          plot.title = element_text(size=20,
                                                    hjust = 0.5,
                                                    face="bold",
                                                    margin = margin(10, 0, 10, 0)))
                ## main
                p <- p + geom_line(aes(time, outcome,
                                       colour = type,
                                       size = type,
                                       linetype = type))

                ## legend
                set.limits = c("tr","co")
                set.labels = c("Treated Average",
                               "Estimated Y(0) Average")
                set.colors = c("red","steelblue")
                set.linetypes = c("solid","longdash")
                set.linewidth = rep(line.width[1],2)
                p <- p + scale_colour_manual(limits = set.limits,
                                             labels = set.labels,
                                             values =set.colors) +
                    scale_linetype_manual(limits = set.limits,
                                          labels = set.labels,
                                          values = set.linetypes) +
                    scale_size_manual(limits = set.limits,
                                      labels = set.labels,
                                      values = set.linewidth) +
                    guides(linetype = guide_legend(title=NULL, ncol=2),
                            colour = guide_legend(title=NULL, ncol=2),
                            size = guide_legend(title=NULL, ncol=2)) 
                    
            } else if  (raw == "band") {
                    
                Y.tr.band <- t(apply(Y.tr.aug, 1, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                ## Y.co.band <- t(apply(Y.co, 1, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                    
                data <- cbind.data.frame("time" = rep(time[show],2),
                                         "outcome" = c(Yb[show,1],
                                                       Yb[show,2]),
                                         "type" = c(rep("tr",nT),
                                                    rep("co",nT)))

                data.band <- cbind.data.frame(time, Y.tr.band)[show,]
                colnames(data.band) <- c("time","tr.lower","tr.upper")
                    
                ## theme 
                p <- ggplot(data)
                if (theme.bw == TRUE) {
                  p <- p + theme_bw()
                }
                p <- p + xlab(xlab) +  ylab(ylab) +
                    geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                    annotate("rect", xmin= time.bf, xmax= Inf,
                             ymin=-Inf, ymax=Inf, alpha = .3) +
                    theme(legend.position = legend.pos,
                          axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                          plot.title = element_text(size=20,
                                                    hjust = 0.5,
                                                    face="bold",
                                                    margin = margin(10, 0, 10, 0)))
                ## main
                p <- p + geom_line(aes(time, outcome,
                                       colour = type,
                                       size = type,
                                       linetype = type))
                ## band
                p <- p + geom_ribbon(data = data.band,
                                     aes(ymin = tr.lower, ymax = tr.upper, x=time),
                                         alpha = 0.15, fill = "red")

                set.limits = c("tr","co","tr.band")
                set.labels = c("Treated Average",
                               "Estimated Y(0) Average",
                                "Treated 5-95% Quantiles")
                set.colors = c("red","steelblue","#FF000030")
                set.linetypes = c("solid","longdash","solid")
                set.linewidth = c(rep(line.width[1],2),4)

                p <- p + scale_colour_manual(limits = set.limits,
                                             labels = set.labels,
                                             values =set.colors) +
                    scale_linetype_manual(limits = set.limits,
                                          labels = set.labels,
                                          values = set.linetypes) +
                    scale_size_manual(limits = set.limits,
                                      labels = set.labels,
                                      values = set.linewidth) +
                    guides(linetype = guide_legend(title=NULL, ncol=2),
                           colour = guide_legend(title=NULL, ncol=2),
                           size = guide_legend(title=NULL, ncol=2)) 
                    
            } else if (raw == "all") { ## plot all the raw data
                    
                data <- cbind.data.frame("time" = rep(time[show],(2 + Ntr)),
                                         "outcome" = c(Yb[show,1],
                                                       Yb[show,2],
                                                       c(Y.tr.aug[show,])),
                                         "type" = c(rep("tr",nT),
                                                    rep("co",nT),
                                                    rep("raw.tr",(Ntr * nT))),
                                          "id" = c(rep("tr",nT),
                                                  rep("co",nT),
                                                  rep(c(x$id.tr),
                                                      each = nT))) 
                ## theme
                p <- ggplot(data)
                if (theme.bw == TRUE) {
                  p <- p + theme_bw()
                }
                p <- p + xlab(xlab) +  ylab(ylab) +
                    geom_vline(xintercept=time.bf,colour=line.color,size = 2) +
                    annotate("rect", xmin= time.bf, xmax= Inf,
                             ymin=-Inf, ymax=Inf, alpha = .3) +
                    theme(legend.position = legend.pos,
                          axis.text.x = element_text(angle = angle, hjust=x.h, vjust=x.v),
                          plot.title = element_text(size=20,
                                                    hjust = 0.5,
                                                    face="bold",
                                                    margin = margin(10, 0, 10, 0))) 
                ## main
                p <- p + geom_line(aes(time, outcome,
                                       colour = type,
                                       size = type,
                                       linetype = type,
                                       group = id))
                ## legend
                set.limits = c("tr","co","raw.tr")
                set.labels = c("Treated Average",
                               "Estimated Y(0) Average",
                               "Treated Raw Data")
                set.colors = c("red","steelblue","#FC8D6280")
                set.linetypes = c("solid","longdash","solid")
                set.linewidth = c(rep(line.width[1],2),line.width[2])
                    
                p <- p + scale_colour_manual(limits = set.limits,
                                             labels = set.labels,
                                             values =set.colors) +
                    scale_linetype_manual(limits = set.limits,
                                          labels = set.labels,
                                          values = set.linetypes) +
                    scale_size_manual(limits = set.limits,
                                      labels = set.labels,
                                      values = set.linewidth) +
                    guides(linetype = guide_legend(title=NULL, ncol=2),
                           colour = guide_legend(title=NULL, ncol=2),
                           size = guide_legend(title=NULL, ncol=2)) 
            }

            ## title
            if (is.null(main) == TRUE) {
                p <- p + ggtitle(maintext)
            } else if (main!="") {
                p <- p + ggtitle(main)
            }

            ## ylim
            if (is.null(ylim) == FALSE) {
                p <- p + coord_cartesian(ylim = ylim)
            }
            suppressWarnings(print(p))

            
        }

    } else if (type == "weights") {

        maintext <- "Weights of the Control Units"
        logw <- log(x$weights.co * x$Ntr)
        logw[logw < wmin] <- wmin
        w.min <- min(logw)
        w.max <- max(logw)
        if (w.min!=w.max) {
            bw <- abs(w.max - w.min)/30
        } else {
            bw <- abs(w.min)/20            
        }
        p <- qplot(logw, col = I("gray70"), binwidth = bw,
          xlab = "Weight (logarithmic)", ylab = "Counts")
        if (theme.bw == TRUE) {
          p <- p + theme_bw()
        }
        p <- p + theme(plot.title = element_text(hjust = 0.5))
                

        if (is.null(xlim) == FALSE) {
            p <- p + coord_cartesian(xlim = xlim)
        } else if (w.min==w.max) {
            p <- p + coord_cartesian(xlim = c(w.min*2,0))
        }
        if (is.null(ylim) == FALSE) {
                p <- p + coord_cartesian(ylim = ylim)
        }

    }  else if (type == "balance") {

        if (is.null(x$matchvar)==TRUE) {
            stop("No covariates being matched on.")            
        }
        if (!stat %in% c("mean","sd")) {
            stop("Wrong specification for \"stat\".")
        }
        if (Ntr == 1 & stat =="sd") {
            stop("Standard deviation does not exist with one treated unit.")
        }


        maintext <- "Covariate Balance"
        bal <- x$bal.table
        nvar <- nrow(bal)
        var <- rep(rownames(bal),2)
        var <- factor(var, levels = var[nvar:1])
        if (stat == "mean") {
           diff <- c(bal$diff.pre, bal$diff.pst)
           if (Ntr > 1) {
              xlab <- "Standardized Difference in Means"
           } else if (Ntr == 1) {
              xlab <- "Difference in Means / |Treated Mean|"
           }
        } else if (stat == "sd") {
           diff <- c((bal$sd.co.pre-bal$sd.tr)/bal$sd.tr, (bal$sd.co.pst-bal$sd.tr)/bal$sd.tr)
           xlab <- "Standardized Difference in Standard Deviations"
        }
        group <- c(rep("Unweighted",nvar), rep("Weighted",nvar))
        newbal <- cbind.data.frame(var,diff,group)
        size <- 1
        stroke <- .8*size
        colors <- c("#E69F00", "#56B4E9")
        p <- ggplot(data = newbal, aes(y = var, x = diff, group = group))        
        if (is.null(xlim) == TRUE) {
            maxv <- max(abs(diff))
            xlim <- c(-maxv, maxv)
        } 
        p <- p + coord_cartesian(xlim = xlim)
        # p <- p + 
        # theme(panel.grid.major = element_line(color = "gray87"),
        #   panel.grid.minor = element_line(color = "gray90"),
        #   panel.background = element_rect(fill = "white", color = "black"),
        #   axis.text.x = element_text(color = "black"),
        #   axis.text.y = element_text(color = "black"),
        #   plot.title = element_text(hjust = 0.5),
        #   legend.title = element_blank(),
        #   legend.position = legend.pos
        #   )
        if (theme.bw == TRUE) {
          p <- p + theme_bw()
        }
        p <- p + theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.position = legend.pos
          )
        p <- p + labs(y = "", x = xlab) + scale_color_manual(values = colors) 
        p <- p + geom_vline(xintercept = 0, linetype = 1, colour = line.color, size = 2)
        p <- p + geom_point(aes(shape = group, color = group),
          size = 4*size, stroke = stroke, na.rm = TRUE)

    }

    ## title
    if (is.null(main) == TRUE) {
        p <- p + ggtitle(maintext)
    } else if (main!="") {
        p <- p + ggtitle(main)
    }

    ## ylim
    if (is.null(ylim) == FALSE) {
        p <- p + coord_cartesian(ylim = ylim)
    }
    suppressWarnings(print(p))
}
