
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
    trim.npre = 0, # drop units with <= certain periods of pre-treatment data
    Y.match.time = NULL,
    Y.match.npre = NULL, # fix the number of pre-periods for balancing when T0s are different
    demean = TRUE, # take out pre-treatment unit mean
    estimator = "meanfirst",  # mean, meanfirst, kernel
    sigma=NULL,
    print.baltable = TRUE, # print out table table
    vce = "jackknife", ## uncertainty estimates
    conf.lvl = 0.95, ## confidence interval
    nsims = 200, ## number of bootstrap runs
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
    trim.npre = 0, # drop units with <= certain periods of pre-treatment data
    Y.match.time = NULL,
    Y.match.npre = NULL, # fix the number of pre-periods for balancing when T0s are different
    demean = TRUE, # take out pre-treatment unit mean
    estimator = "meanfirst",  # mean, meanfirst, kernel
    sigma=NULL,
    print.baltable = TRUE, # print out table table
    vce = "jackknife", ## uncertainty via bootstrap
    conf.lvl = 0.95, ## confidence interval
    nsims = 200, ## number of bootstrap runs
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
                          X.avg.time = X.avg.time, index = index, trim.npre = trim.npre,
                          Y.match.time= Y.match.time, Y.match.npre = Y.match.npre,  
                          demean = demean, estimator = estimator, 
                          sigma = sigma, 
                          print.baltable = print.baltable, 
                          vce = vce, conf.lvl = conf.lvl, nsims = nsims, 
                          parallel = parallel, cores = cores, 
                          seed = seed)
    
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
    trim.npre = 0, # drop units with <= certain periods of pre-treatment data
    Y.match.time = NULL,
    Y.match.npre = NULL, # fix the number of pre-periods for balancing when T0s are different
    demean = TRUE, # take out pre-treatment unit mean
    estimator = "meanfirst",  # mean, meanfirst, kernel
    sigma=NULL,
    print.baltable = TRUE, # print out table table
    vce = "jackknife", ## uncertainty via bootstrap
    conf.lvl = 0.95, ## confidence interval
    nsims = 200, ## number of bootstrap runs
    parallel = TRUE, ## parallel computing
    cores = 4,
    seed = 1234
    ) {

    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------##  

    
    if (class(data)[1] == "tbl_df") {
        #warning("Transforming a tibble into a data frame.")
        data <- as.data.frame(data)
    }
    if (is.data.frame(data)==FALSE) {
        stop("Not a data frame")
    }
    data <- droplevels(data)

    ## index
    if (length(index) != 2 | sum(index %in% colnames(data)) != 2) {
        stop("\"index\" option misspecified. Try, for example, index = c(\"unit.id\", \"time\").")
    }   

    if (vce == "boot") {vce <- "bootstrap"}
    if (vce == "jack") {vce <- "jackknife"}
    if (vce == "fixed") {vce <- "fixed.weights"}

    

    if (is.null(Y.match.time)==FALSE) {
        if (Y.match.time[1] == "none") {
            Y.match.pre <- 0
            Y.match.time <- NULL
        }         
    }

    if (is.null(sigma)==FALSE) {
        if (is.numeric(sigma)==FALSE) {
            stop("\"sigma\" needs to be numeric; the default is 2.")
        }
    }

    if (is.null(nsims)==TRUE) {
        nsims <- 200
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

    ## check balanced panel
    if (nrow(data) != TT*N) {
        stop("Data are not balanced or \"index\" does not uniquely identity an observation.")
    }
    
    ##treatment indicator
    D.sav <- D<- matrix(data[,Dname],TT,N)

    ## once treated, always treated
    D <- apply(D, 2, function(vec){cumsum(vec)})
    T0 <- TT - D[TT,] # a vector, number of pre-treatment periods for each unit

    ## drop units with too few pre-treatment periods
    id.drop <- which(T0 <= trim.npre)
    N.drop <- length(id.drop)
    D <- ifelse(D > 0, 1, 0)
    if (sum(abs(D-D.sav))!=0) {
        cat("\nTreatment status changed to \"treated\" after a unit has even been treated; in other words, no switch on-and-off is allowed.\n")
    }
    if (N.drop>0) {
        N <- N - N.drop
        D <- D[,-id.drop, drop = FALSE]
        data <- data[rep(T0,each = TT)>trim.npre,]
        units <- units[-id.drop]
        T0 <- T0[-id.drop]
        cat(paste0("\nDrop ",length(id.drop)," units with ",trim.npre," or fewer pre-treatment periods.\n"))
    }
    

    ## treatment
    treat <-ifelse(D[TT,]==1, 1, 0)     # cross-sectional: treated unit
    id.tr <- which(treat == 1)
    id.co <- which(treat == 0)
    Ntr <- length(id.tr)
    Nco <- length(id.co) 
    if (Ntr == 0) {
        stop("No treated units remain.")
    } 
    if (Nco == 0) {
        stop("No control units remain.")
    }

    ## check the number of treated units
    if (Ntr <= 5) {
        cat("Too few treated unit(s). Uncertainty estimates not provided.\n")
        vce <- "none"
    }    

    ## treatment timing
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
    Y.var <- paste0(Yname, Ttot) ## outcome variable names (wide form)
    colnames(outcome) <- Y.var ## including both pre and post

    ## covariates (allow missing, but non-missing values have to be same for each unit)
    if (class(data[,id])!="factor") { ## to avoid an error with ddply
        data[,id] <- as.factor(data[,id])        
    }
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
        data.wide <- cbind.data.frame(id = 1:N, unit = units, treat = treat, T0 = T0, outcome, Xvar)
    } else {
        data.wide <- cbind.data.frame(id = 1:N, unit = units, treat = treat, T0 = T0, outcome)
    } 


    #######################
    ## balancing
    #######################

    if (sameT0 == TRUE) {
        bal.out <- tjbal.single(data = data.wide, Y = Yname, D = "treat", X = Xname,
            Y.match.time = Y.match.time, Y.match.npre = Y.match.npre, 
            Ttot = Ttot, unit = "id", 
            demean = demean, estimator = estimator, sigma = sigma, 
            print.baltable = print.baltable,
            vce = vce, conf.lvl = conf.lvl,
            nsims = nsims, parallel = parallel, cores = cores, seed = seed)         
    } else {        
        bal.out <- tjbal.multi(data = data.wide, Y = Yname, D = "treat", X = Xname,
            Y.match.time = Y.match.time, Y.match.npre = Y.match.npre, 
            Ttot = Ttot, unit = "id", 
            demean = demean, estimator = estimator, sigma = sigma, 
            vce = vce, conf.lvl = conf.lvl,
            nsims = nsims, parallel = parallel, cores = cores, seed = seed)  
    } 

    out <- c(list(sameT0 = sameT0, index = index, Yname = Yname), bal.out)
    out$call <- match.call()
    class(out) <- "tjbal"
    return(out)
}


###################################################
#### tjbal.multi: allow multiple treatment timing
###################################################


tjbal.multi <- function(
    data, ## data in wide form
    Y,
    D,
    X,
    Y.match.time = NULL,
    Y.match.npre = NULL, # fix the number of pre-periods for balancing when T0s are different
    Ttot,
    unit,
    demean = FALSE, # take out pre-treatment unit mean
    estimator = "mean",
    sigma = NULL, ## tuning parameters   
    vce = "jackknife", ## uncertainty via jackknife
    conf.lvl = 0.95, ## confidence interval
    nsims = 200, ## number of bootstrap runs
    parallel = TRUE, ## parallel computing
    cores = 4,
    seed = 1234  
    ) { 

    set.seed(seed)

    TT <- length(Ttot)
    id.tr <- which(data$tr == 1)
    id.co <- which(data$tr == 0)
    Ntr <- length(id.tr)
    Nco <- length(id.co)
    N <- Ntr + Nco
    Y.var  <- paste0(Y, Ttot)
    units <- as.character(unique(data$unit))

    ## parse multiple treatment timing
    T0.all <- data$T0
    names(T0.all) <- units
    T0.tr <- T0.all[id.tr] # T0 for all treated units
    tb.T0 <- table(T0.tr)
    T0.unique <- as.numeric(names(tb.T0))
    T0.count <- as.numeric(tb.T0)
    T0.max <- max(T0.tr)
    T0.min <- min(T0.tr)
    T0.names <- paste0("T0 = ",Ttot[T0.unique])
    nT0 <- length(T0.unique)

    ## recenter based on treatment timing
    time.adj <- c(-(T0.max -1) : (TT-T0.min))
    time.adj.max <- TT-T0.min
    TT.adj <- length(time.adj)

    ## storage
    sub.weights.co <- matrix(NA, length(id.co), nT0) # id.co * T0.unique
    rownames(sub.weights.co) <- units[id.co]
    colnames(sub.weights.co) <- T0.names
    sub.att <- sub.Ytr.avg <- sub.Yct.avg <- matrix(NA, TT, nT0)
    sub.att.avg <- rep(NA, nT0)
    rownames(sub.att) <- rownames(sub.Ytr.avg) <- rownames(sub.Yct.avg) <- Ttot
    names(sub.att.avg) <- colnames(sub.att) <- colnames(sub.Ytr.avg) <- colnames(sub.Yct.avg) <- T0.names
    ## number of treated units for each subgroup
    sub.Ytr.adj <- sub.ntr.pst <- sub.ntr <- sub.att.adj <- matrix(0, TT.adj, nT0) 
    se.sub.att.adj <- matrix(NA, TT.adj, nT0) 
    # naming
    rownames(sub.Ytr.adj) <-rownames(sub.ntr.pst) <- rownames(sub.ntr) <- rownames(se.sub.att.adj) <- rownames(sub.att.adj) <- time.adj
    colnames(sub.Ytr.adj) <-colnames(sub.ntr.pst) <- colnames(sub.ntr) <- colnames(se.sub.att.adj) <- colnames(sub.att.adj) <- T0.names
     

    ## save subsets of the original data, K matrix, and variables to be balanced on
    constraint.list <- bal.table.list <- matchvar.list <- K.list <- data.list <- vector("list", length = nT0)  
    ## other storage
    success <- sigmas <- ndims <- rep(0, nT0) 
    bias.ratios <- rep(1, nT0) 
    kernels <- rep(NA, nT0)
    

    cat("Balancing...\n")        
    for (i in 1:nT0) {

        T0 <- T0.unique[i]
        cat(paste0("Subgroup T0 = ",T0,": "))
        id.tr.one <- which(T0.all == T0)    

        out <- tjbal.core(data = data, Y = Y, X = X, Y.match.time = Y.match.time, 
            Y.match.npre = Y.match.npre, Ttot = Ttot, T0 = T0, id.tr = id.tr.one, id.co = id.co,
            demean = demean, estimator = estimator, sigma = sigma, 
            info = FALSE) 

        # saved results
        att <- out$att
        Y.var <- out$Y.var
        w <- out$w
        matchvar <- out$matchvar
        Y.target <- out$Y.target
        Y.target.pst <- out$Y.target.pst
        data.list[[i]] <- out$data.wide
        if (is.null(matchvar)==FALSE) {
            matchvar.list[[i]] <- matchvar
            success[i] <- out$success
            if (success[i] == TRUE) {
                bias.ratios[i] <- out$bias.ratio
                ndims[i] <- out$ndims
                if (out$bal.type == "kbal") {
                    K.list[[i]] <- out$kbal.out$K           
                    sigmas[i] <- out$b                    
                    kernels[i] <- 1 
                    if (estimator == "meanfirst") {
                        constraint.list[[i]] <- out$constriant
                    }                   
                } else {
                    kernels[i] <- 0
                }                
            }
        }

        # cat(paste0("bias.ratio = ",sprintf("%.4f",bias.ratios[i]),"; num.dims = ",ndims[i]))
        # if (ndims[i]==0) {cat(" (Fail)")}
        # cat("\n")                        
        
        sub.weights.co[,i] <- out$weights.co
        sub.att[,i] <- att
        sub.att.avg[i] <- out$att.avg
        sub.Ytr.avg[,i] <- out$Y.bar[,1]
        sub.Yct.avg[,i] <- out$Y.bar[,2]     
       

        ## save ATT (realigned based on T0)        
        fill.start <- T0.max-T0+1
        fill.end <- fill.start + length(att) -1 
        sub.ntr[fill.start:fill.end, i] <- T0.count[i]   
        sub.att.adj[fill.start:fill.end, i] <- att
        sub.Ytr.adj[fill.start:fill.end, i] <- sub.Ytr.avg[,i]

        ## balance table
        if (is.null(matchvar)==FALSE) {          
            bal.table.list[[i]] <- out$bal.table
        }       
    }
    ntreated <- rowSums(sub.ntr) # how the number of units changes over adjusted time
    att <- rowSums(sub.att.adj * sub.ntr)/ntreated
    names(ntreated) <- names(att) <- time.adj

    # average Y and average counterfactual (realigned)
    Y.bar.tr <- rowSums(sub.Ytr.adj * sub.ntr)/ntreated
    Y.bar.ct <- Y.bar.tr - att
    Y.bar <- cbind(Y.bar.tr, Y.bar.ct)

    # this matrix count the number of treated units for each T0 (post-treatment)
    sub.ntr.pst <- sub.ntr
    sub.ntr.pst[time.adj<=0,] <- 0 
    
    ## average effect (weighted by obs)
    att.avg <- sum(sub.att.adj * sub.ntr.pst, na.rm = TRUE)/sum(sub.ntr.pst)

    ## weights.co
    weights.co <- apply(sub.weights.co,1,weighted.mean,T0.count)
    names(weights.co) <- units[id.co]

    ## group statistics
    group.stats <- cbind(T0 = T0.unique, time = Ttot[T0.unique], Ntr = T0.count, 
        success = success, sigma = floor(sigmas), bias.ratio = bias.ratios)

    
    #######################################
    ## Uncertainty Estimates via Jackknife
    #######################################    

    #######################
    ## Jackknife
    #######################

    if (vce == "jackknife") {

        # number of jackknife runs
        if (nsims > Ntr) {
            njacks <- Ntr
        } else {
            njacks <- nsims
        }       

        cat("\nJackknife...\n")
        drop.id.pos <- sample(1:Ntr, njacks, replace = FALSE)
        drop.id <- id.tr[drop.id.pos]
        drop.id.T0s <- T0.tr[drop.id.pos]   
        drop.id.list <- split(drop.id, drop.id.T0s) # put drop id into different T0 groups

        # save att and att.avg for each subgroup when one of the units from the group is dropped
        sub.att.jack <- vector("list", length = nT0)
        sub.att.avg.jack <- vector("list", length = nT0)
        names(sub.att.avg.jack) <- names(sub.att.jack) <- colnames(sub.att)
      
        ## prepare for parallel computing
        if (parallel == TRUE) {                        
            if (is.null(cores) == TRUE) {
                cores <- detectCores()
            }
            para.clusters <- makeCluster(cores)
            registerDoParallel(para.clusters)            
        }

        k <- 1
        for (i in 1:nT0) {
            
            T0 <- T0.unique[i]
            cat("\nDropping units from Subgroup T0 =",T0)
            drop.id.oneT0 <- drop.id.list[[as.character(T0)]]

            # this matrix count the number of treated units for each T0
            sub.ntr.tmp <- sub.ntr
            sub.ntr.tmp[,i] <- max(0,sub.ntr[,i]-1)
            # this matrix count the number of treated units for each T0 (post-treatment)
            sub.ntr.pst.tmp <- sub.ntr.pst
            sub.ntr.pst.tmp[,i] <- max(0,sub.ntr.pst.tmp[,i]-1)

            if (length(drop.id.oneT0) == 1) {
                tmp <- matrix(NA, 1, 1)
                colnames(tmp) <- paste("Drop:",units[drop.id.oneT0])
                sub.att.jack[[i]] <- tmp
                sub.att.avg.jack[[i]] <- tmp
                # att.adj for this run
                sub.att.adj.tmp <- sub.att.adj
                sub.att.adj.tmp[,i] <- 0
                # counter
                k <- k + 1
            } else { # more than one unit in this group 

                one.jack <- function(data, K, id, ndims, matchvar, sigma, kernel, constraint) {
                        N <- nrow(data)
                        Ntr <- N - Nco
                        id <- which(data$id == id) # translate original id to new id
                        N <- Ntr + Nco
                        # weights: treated add up to 1; controls add up to -1; sum is zero
                        w.jack <-  rep(1/(Ntr-1), (N-1))
                        if (is.null(matchvar) == TRUE | ndims == 0) { # no reweighting                         
                            w.jack[(Ntr+1):N] <- rep(-1/Nco, Nco)
                        } else {
                            tmp <- capture.output(                               
                                kbal.jack <- suppressWarnings(kbal(
                                    allx = data[-id, matchvar], 
                                    K = K[-id,],
                                    constraint = constraint[-id,], 
                                    treatment = data[-id, D],
                                    linkernel = (1-kernel), b = sigma, 
                                    printprogress = FALSE,                    
                                    minnumdims = max(0,ndims-5), maxnumdims = ndims,
                                    sampledinpop=FALSE))
                            , file = NULL)                
                            w.jack[Ntr:(N-1)] <- kbal.jack$w[Ntr:(N-1)]/Nco*(-1)  # controls add up to -1;   
                        }                      
                        att <- apply(data[-id, Y.target] * w.jack, 2, sum)
                        att.avg <- mean(att[(T0+1):TT])
                        out <- list(att = att, att.avg = att.avg)
                        return(out)            
                }

                nid <- length(drop.id.oneT0) # number of treated units with one T0

                # storage att and att.avg for this T0
                sub.att.oneT0 <- matrix(NA,TT,nid)
                rownames(sub.att.oneT0) <- Ttot
                colnames(sub.att.oneT0) <- paste("Drop:",units[drop.id.oneT0])
                sub.att.avg.oneT0 <- rep(NA,nid) 
                names(sub.att.avg.oneT0) <- paste("Drop:",units[drop.id.oneT0])
                # copy from previous result (will change one column)
                sub.att.adj.tmp <- sub.att.adj
                sub.att.adj.tmp[,i] <- 0
                

                ## computing
                if (parallel == TRUE & nid >= 8) {
                    
                    ## start    
                    cat(" -- Parallel computing...") 
                    jack.out <- foreach(j=1:nid, 
                        .inorder = FALSE,
                        .packages = c("kbal")
                        ) %dopar% {
                        return(
                            one.jack(data = data.list[[i]], 
                            K = K.list[[i]], id = drop.id.oneT0[j], 
                            ndims = ndims[i], matchvar = matchvar.list[[i]], 
                            sigma = sigmas[i], kernel = kernels[i], 
                            constraint = constraint.list[[i]])
                            )
                    }                    
                    ## save results
                    for (j in 1:nid) { 
                        sub.att.oneT0[,j] <- jack.out[[j]]$att
                        sub.att.avg.oneT0[j] <- jack.out[[j]]$att.avg 
                        # counter
                        k <- k + 1                
                    }
                } else { ## single core
                    for (j in 1:nid) {  

                        jack <- one.jack(data = data.list[[i]], 
                            K = K.list[[i]], id = drop.id.oneT0[j], 
                            ndims = ndims[i], matchvar = matchvar.list[[i]], 
                            sigma = sigmas[i], kernel = kernels[i], 
                            constraint = constraint.list[[i]])

                        sub.att.oneT0[,j] <- jack$att
                        sub.att.avg.oneT0[j] <- jack$att.avg 
                        # counter
                        k <- k + 1 
                        ## report progress
                        if (j%%50==0)  {cat(".")} 
                    }  # end of single core loop
                }  # end of T0 with more than one treated unit
                sub.att.jack[[i]] <- sub.att.oneT0
                sub.att.avg.jack[[i]] <- sub.att.avg.oneT0   

            } # end of more than one unit case            
                      

        } ## end of loop around different T0s

        ## saved objects
        # sub.att.jack -- for treated units in the group
        # sub.att.avg.jack -- average within group
        # att.jack (adjusted) -- all treated units
        # att.avg.jack -- average over all treated units

        if (parallel == TRUE) {
            stopCluster(para.clusters)
        }


        ##################################
        ## SE, z-scores p values, and CI
        ##################################

        ## critical value for a two-sided test
        c.value <- qnorm(0.5 + conf.lvl/2)

        ## uncertainty for subgroups
        CI.sub.att.low <- CI.sub.att.high <- pvalue.sub.att <- z.sub.att <- se.sub.att <- matrix(NA, TT, nT0)
        se.sub.att.avg <- rep(NA, nT0)
        rownames(se.sub.att) <- Ttot
        names(se.sub.att.avg) <- colnames(se.sub.att) <- T0.names
        for (i in 1:nT0) {
            Ntr.oneT0 <- ncol(sub.att.jack[[i]])
            if (Ntr.oneT0>1) {
                se.sub.att[,i] <- apply(sub.att.jack[[i]], 1, function(vec) sd(vec, na.rm=TRUE)) * (Ntr.oneT0-1)/sqrt(Ntr.oneT0)
                se.sub.att.avg[i] <- sd(sub.att.avg.jack[[i]], na.rm=TRUE) * (Ntr.oneT0-1)/sqrt(Ntr.oneT0)
                # z, pvalue, and CIs for sub.att
                z.sub.att[,i] <- sub.att[,i]/se.sub.att[,1]
                pvalue.sub.att[,i] <- (1 - pnorm(abs(z.sub.att[,i])))*2
                CI.sub.att.low[,i] <- sub.att[,i] - c.value * se.sub.att[,1]
                CI.sub.att.high[,i] <- sub.att[,i] + c.value * se.sub.att[,1]
            }
        }

        # se.att
        group.above1 <- which(is.na(se.sub.att.avg)==FALSE) # more than 1 unit        
        for (i in group.above1) {
            T0 <- T0.unique[i]
            fill.start <- T0.max-T0+1
            fill.end <- fill.start + TT -1 
            se.sub.att.adj[fill.start:fill.end, i] <- se.sub.att[,i]            
        } # now we have: sub.att.adj and se.sub.att.adj 
        vce.sub.att.adj <- se.sub.att.adj^2
        sum.vce.bytime <- apply(vce.sub.att.adj[,group.above1] * sub.ntr[,group.above1]^2, 1, sum, na.rm = TRUE)
        sum.obs.bytime <- apply(sub.ntr[,group.above1], 1, sum)
        sum.vce.bytime[which(sum.vce.bytime == 0)] <- sum.obs.bytime[which(sum.obs.bytime == 0)] <- NA
        se.att <- sqrt(sum.vce.bytime/sum.obs.bytime^2)


        ## se.att.avg
        sub.obs <- apply(sub.ntr.pst[,group.above1], 2, sum) # total number of obs by group
        #att.avg <- sum(sub.att.avg[group.above1] * sub.obs)/sum(sub.obs) # different from using all observations
        vce.sub.att.avg <- se.sub.att.avg[group.above1]^2
        se.att.avg <- sqrt(sum(vce.sub.att.avg * sub.obs^2)/(sum(sub.obs)^2))

        ## z-score
        z.att <- att/se.att
        z.att.avg <- att.avg/se.att.avg
        z.sub.att.avg <- sub.att.avg/se.sub.att.avg

        ## two-sided p-value
        pvalue.att <- (1 - pnorm(abs(z.att)))*2
        pvalue.att.avg <- (1 - pnorm(abs(z.att.avg)))*2
        pvalue.sub.att.avg <- (1 - pnorm(abs(z.sub.att.avg)))*2

        ## confidence intervals
        CI.att <- cbind(att - c.value * se.att, att + c.value * se.att)
        CI.att.avg <- c(att.avg - c.value * se.att.avg, att.avg + c.value * se.att.avg)
        CI.sub.att.avg <- cbind(sub.att.avg - c.value * se.sub.att.avg, sub.att.avg + c.value * se.sub.att.avg)


        ## put everything together
        est.names <- c("ATT", "S.E.", "z-score", "CI.lower", "CI.upper","p.value", "n.Treated")
        est.att <- cbind(att, se.att, z.att, CI.att, pvalue.att, ntreated = ntreated)
        est.att[abs(est.att)<1e-5] <- 0
        colnames(est.att) <- est.names
        
        ## average effect
        est.att.avg <- t(as.matrix(c(att.avg, se.att.avg, z.att.avg, CI.att.avg, pvalue.att.avg)))
        colnames(est.att.avg) <- est.names[1:6]

        ## average effect for subgroups
        est.sub.att.avg <- cbind(sub.att.avg, se.sub.att.avg, z.sub.att.avg, 
            CI.sub.att.avg, pvalue.sub.att.avg, ntreated = T0.count)
        est.sub.att.avg[abs(est.sub.att.avg)<1e-5] <- 0
        colnames(est.sub.att.avg) <- est.names
        
        ## average effect over time for subgroups
        est.sub.att <- array(NA, dim = c(TT,7,nT0),
            dimnames = list(Ttot, est.names, paste("T0 =",Ttot[T0.unique]))) 
        for (i in 1:nT0) {
            est.sub.att[,,i] <- t(as.matrix(c(sub.att[,i], se.sub.att[,i], z.sub.att[,i], 
                CI.sub.att.low[,i], CI.sub.att.high[,i], pvalue.sub.att[,i], rep(T0.count[i],TT))))
        }
        
        ## storage
        out.inference <- list(
            est.att = round(est.att,4), 
            est.att.avg = round(est.att.avg,4), 
            est.sub.att = round(est.sub.att,4),
            est.sub.att.avg = round(est.sub.att.avg,4)
            ) 
    } 

    # replace 0 with NA for att.adj
    sub.att.adj[sub.ntr==0] <- NA

    
    #####################
    ## Save Results
    #####################


    ## save results
    out <- list(
        data.wide = data,
        id.tr = id.tr,
        id.co = id.co,
        Y.var = Y.var,
        matchvar.list = matchvar.list,
        Ttot = Ttot,
        N = N,
        Ntr = Ntr,
        Nco = Nco,            
        T0 = T0.unique, 
        T0.all = T0.all,
        T0.tr = T0.tr,
        weights.co = weights.co,
        Y.bar = Y.bar,
        sub.weights.co = sub.weights.co,
        sub.Ytr.avg = sub.Ytr.avg,
        sub.Yct.avg = sub.Yct.avg,
        sub.att = sub.att,
        sub.ntr = sub.ntr,
        sub.att.adj = sub.att.adj,
        ntreated = ntreated,
        att = att,
        att.avg = att.avg,
        success = success,
        bias.ratios = bias.ratios,
        group.stats = group.stats,
        bal.table.list = bal.table.list
        )

    if (vce == "jackknife") {
        out <- c(out, out.inference)
    }     
    
    return(out)
}





###################################################
#### tjbal single: same treatment timing
###################################################


tjbal.single <- function(
    data, ## data in wide form
    Y,
    D,
    X,
    Y.match.time = NULL,
    Y.match.npre = NULL, # fix the number of pre-periods for balancing when T0s are different
    Ttot,
    unit,
    demean = FALSE, # take out pre-treatment unit mean
    estimator = "mean",
    sigma = NULL, ## tuning parameters
    test = TRUE, ## test different sigmas
    nsigma = 5,
    print.baltable = FALSE,
    vce = "jackknife", ## uncertainty via jackknife
    conf.lvl = 0.95, ## confidence interval
    nsims = 500, ## number of bootstrap runs
    parallel = TRUE, ## parallel computing
    cores = 4,
    seed = 1234  
    ) { 

    set.seed(seed)

    TT <- length(Ttot)
    id.tr <- which(data$treat == 1)
    id.co <- which(data$treat == 0)
    Ntr <- length(id.tr)
    Nco <- length(id.co)
    N <- Ntr + Nco 
    T0 <- unique(data$T0[id.tr])
    Tpre <- Ttot[1:T0]
    Tpst <- Ttot[(T0+1):TT] 

    if (vce == "jackknife") {
        if (nsims > Ntr) {
            njacks <- Ntr
        } else {
            njacks <- nsims
        }
    } 

    out <- tjbal.core(data = data, Y = Y, X = X, Y.match.time = Y.match.time, 
        Y.match.npre = Y.match.npre, Ttot = Ttot, T0 = T0, id.tr = id.tr, id.co = id.co,
        demean = demean, estimator = estimator, sigma = sigma, 
        info = TRUE) 
   
    # saved results
    att <- out$att
    att.avg <- out$att.avg  
    Y.var <- out$Y.var
    Y.bar <- out$Y.bar
    w <- out$w
    weights.co <- out$weights.co
    matchvar <- out$matchvar
    Y.target <- out$Y.target
    Y.target.pst <- out$Y.target.pst
    ndims <- out$ndims
    K <- out$kbal.out$K
    constraint <- out$constraint
    if (out$bal.type == "mbal") {
        kernel <- 0
    } else if (out$bal.type == "kbal") {
        kernel <- 1
    }
    data <- out$data.wide
    id.tr <- 1:Ntr
    id.co <- (Ntr + 1): N

    if (print.baltable == TRUE && is.null(matchvar)==FALSE) {
        cat("\nBalance Table\n")
        print(round(out$bal.table, 4))   
    }     
       

    ###########################
    ## Uncertainty Estimates
    ###########################    

    
    #################
    ## Fixed weights
    #################

    if (vce == "fixed.weights") {
        ## Storing estimates
        att.sims<-matrix(0,TT,nsims)
        att.avg.sims<-matrix(0,nsims,1)          
        for (j in 1:nsims) { 
            sample.id <- c(sample(id.tr, Ntr, replace = TRUE),
                sample(id.co, Nco, replace = TRUE))
            w.boot <- w[sample.id]
            w.boot[1:Ntr] <- w.boot[1:Ntr]/sum(w.boot[1:Ntr]) # add up to 1
            w.boot[(Ntr+1):N] <- w.boot[(Ntr+1):N]/sum(w.boot[(Ntr+1):N]) * (-1) # add up to -1
            att.sims[,j]<-apply(data[sample.id, Y.target] * w.boot, 2, sum)
            att.avg.sims[j,]<-mean(att.sims[(T0+1):TT,j])
        }
        cat("\n")
        ## standard errors
        se.att <- apply(att.sims, 1, function(vec) sd(vec, na.rm=TRUE))
        se.att.avg <- sd(att.avg.sims, na.rm=TRUE)
    }

    ###############
    ## Bootstrap
    ###############

    if (vce == "bootstrap") {

        
        ## Simulation
        cat("\nBootstrapping... \n") 

        one.boot <- function() {
            ## weights: treated add up to 1; controls add up to -1; sum is zero
            sample.id <- c(sample(id.co, Nco, replace = TRUE),sample(id.tr,Ntr, replace = TRUE))
            w.boot <- rep(1/Ntr, N)
            if (is.null(matchvar) == TRUE) { # no reweighting                         
                w.boot[1:Nco] <- rep(-1/Nco, Nco)
            } else {
                data.tmp <- data[sample.id, ]
                tmp <- capture.output(                                
                    kbal.boot <- suppressWarnings(kbal(allx = data.tmp[, matchvar], 
                        treatment = data.tmp[, D],
                        K = K[sample.id, ],
                        constraint = constraint[sample.id,], 
                        linkernel = (1-kernel),  
                        printprogress = FALSE, sampledinpop = FALSE,
                        minnumdims = max(0,ndims-5), maxnumdims = ndims))
                , file = NULL)                
                w.boot[1:Nco] <- kbal.boot$w[1:Nco]/Nco*(-1)  # controls add up to -1;   
            }
            att <- apply(data[sample.id, Y.target] * w.boot, 2, sum)
            att.avg <- mean(att[(T0+1):TT])
            out <- list(att = att, att.avg = att.avg)
            return(out)            
        }

        ## Storing bootstrapped estimates
        att.sims<-matrix(0,TT,nsims)
        att.avg.sims<-matrix(0,nsims,1)  

        ## computing
        if (parallel == TRUE) {
            ## prepare
            if (is.null(cores) == TRUE) {
                cores <- detectCores()
            }
            para.clusters <- makeCluster(cores)
            registerDoParallel(para.clusters)
            ## start    
            cat("Parallel computing...") 
            boot.out <- foreach(j=1:nsims, 
                .inorder = FALSE,                
                .packages = c("kbal")
                ) %dopar% {
                return(one.boot())
            }
            stopCluster(para.clusters)
            ## save results
            for (j in 1:nsims) { 
                att.sims[,j]<-boot.out[[j]]$att
                att.avg.sims[j,]<-boot.out[[j]]$att.avg                  
            } 
        } else { ## single core
            for (j in 1:nsims) { 
                boot <- one.boot() 
                att.sims[,j]<-boot$att
                att.avg.sims[j,]<-boot$att.avg                
                ## report progress
                if (kernel == FALSE) {
                    if (j%%50==0)  {cat(".")}
                } else {
                    cat(j,"\n")
                } 
            }  
        }
        # end of bootstrapping
        cat("\n")      
      

        ## standard errors
        se.att <- apply(att.sims, 1, function(vec) sd(vec, na.rm=TRUE))
        se.att.avg <- sd(att.avg.sims, na.rm=TRUE)       
        

    }

    #######################
    ## Jackknife
    #######################

    if (vce == "jackknife") {

        cat("\nJackknife... \n")
        drop.id <- sample(id.tr, njacks, replace = FALSE)                 
        
        one.jack <- function(id) {
            ## weights: treated add up to 1; controls add up to -1; sum is zero
            if (vce == "jackknife") {
                sample.id <- c(id.co, setdiff(id.tr,id)) # drop one treated unit each time 
                w.jack <-  rep(1/(Ntr-1), (N-1))                  
            }
            if (is.null(matchvar) == TRUE) { # no reweighting                         
                w.jack[1:Nco] <- rep(-1/Nco, Nco)
            } else {
                data.tmp <- data[sample.id, ]
                tmp <- capture.output(
                    kbal.jack <- suppressWarnings(kbal(allx = data.tmp[, matchvar], 
                        K = K[sample.id, ],
                        constraint = constraint[sample.id,], 
                        treatment = data.tmp[, D],
                        linkernel = (1-kernel), 
                        printprogress = FALSE,                    
                        minnumdims = max(0,ndims-5), 
                        maxnumdims = ndims))            
                , file = NULL)    
                w.jack[1:Nco] <- kbal.jack$w[1:Nco]/Nco*(-1)  # controls add up to -1;   
            }                      
            ## ATT
            att <- apply(data[sample.id, Y.target] * w.jack, 2, sum)
            att.avg <- mean(att[(T0+1):TT])
            out <- list(att = att, att.avg = att.avg)
            return(out)            
        }

        ## Storing jackknife estimates
        att.sims<-matrix(0,TT,njacks)
        att.avg.sims<-matrix(0,njacks,1)  

        ## computing
        if (parallel == TRUE) {
            ## prepare
            if (is.null(cores) == TRUE) {
                cores <- detectCores()
            }
            para.clusters <- makeCluster(cores)
            registerDoParallel(para.clusters)
            ## start    
            cat("Parallel computing...") 
            jack.out <- foreach(j=1:njacks, 
                .inorder = FALSE,
                .packages = c("kbal")
                ) %dopar% {
                return(one.jack(drop.id[j]))
            }
            stopCluster(para.clusters)
            ## save results
            for (j in 1:njacks) { 
                att.sims[,j]<-jack.out[[j]]$att
                att.avg.sims[j,]<-jack.out[[j]]$att.avg                  
            } 
        } else { ## single core
            for (j in 1:njacks) { 
                jack <- one.jack(drop.id[j]) 
                att.sims[,j]<-jack$att
                att.avg.sims[j,]<-jack$att.avg               
            }  
        }
        ## end of bootstrapping
        cat("\n")
        
        #### SE and CIs ####
        se.att <- apply(att.sims, 1, function(vec) sd(vec, na.rm=TRUE)) * (Ntr-1)/sqrt(Ntr)
        se.att.avg <- sd(att.avg.sims, na.rm=TRUE) * (Ntr-1)/sqrt(Ntr)        

    } # end of jackknife

    #############################
    ## z-scores p values, and CI
    #############################

    if (vce %in% c("fixed.weights","bootstrap","jackknife")) {
    
        ## z-score
        z.att <- att/se.att
        z.att.avg <- att.avg/se.att.avg
        
        ## two-sided p-value
        pvalue.att <- (1 - pnorm(abs(z.att)))*2
        pvalue.att.avg <- (1 - pnorm(abs(z.att.avg)))*2

        ## critical value for a two-sided test
        c.value <- qnorm(0.5 + conf.lvl/2)

        ## confidence intervals
        CI.att <- cbind(att - c.value * se.att, att + c.value * se.att)
        CI.att.avg <- c(att.avg - c.value * se.att.avg, att.avg + c.value * se.att.avg)


        ## put everything together
        est.att <- cbind(att, se.att, z.att, CI.att, pvalue.att, ntreated = rep(Ntr,TT))
        est.att[abs(est.att)<1e-5] <- 0
        colnames(est.att) <- c("ATT", "S.E.", "z-score", "CI.lower", "CI.upper","p.value", "n.Treated")
        rownames(est.att) <- Ttot    

        ## average effect
        est.att.avg <- t(as.matrix(c(att.avg, se.att.avg, z.att.avg, CI.att.avg, pvalue.att.avg)))
        colnames(est.att.avg) <- c("ATT", "S.E.", "z-score", "CI.lower", "CI.upper", "p.value")

        ## storage
        out.inference <- list(
            est.att = round(est.att,4), 
            est.att.avg = round(est.att.avg,4), 
            att.sims = att.sims, 
            att.avg.sims = att.avg.sims
            ) 
    }  

    
    #####################
    ## Save Results
    #####################

    out <- c(out, list(
            id.tr = id.tr,
            id.co = id.co,
            Y.var = Y.var,
            Ttot = Ttot,
            Tpre = Tpre,
            Tpst = Tpst,
            T0 = T0,
            N = N,
            Ntr = Ntr,
            Nco = Nco,                
            ntreated = rep(Ntr,TT)
            ))    
    
    if (vce %in% c("fixed.weights","bootstrap","jackknife")) {
        out <- c(out, out.inference)
    }       
    return(out)

}

###################################################
#### tjbal core function
###################################################


tjbal.core <- function(
    data, ## data in wide form
    Y,
    X,
    Y.match.time = NULL,
    Y.match.npre = NULL, # fix the number of pre-periods for balancing when T0s are different
    Ttot,
    T0,
    id.tr,
    id.co,
    demean = FALSE, # take out pre-treatment unit mean
    estimator,
    sigma = NULL, ## tuning parameters
    info ## show information
    ) { 


    TT <- length(Ttot)
    Ntr <- length(id.tr)
    Nco <- length(id.co)
    N <- Ntr + Nco 
    Tpre <- Ttot[1:T0]
    Tpst <- Ttot[(T0+1):TT] 

    if (is.null(Y.match.npre)==FALSE) {
        if (Y.match.npre == 0) {
            Y.match.time <- NULL
        } else {
            Y.match.time <- Ttot[max(1,(T0-Y.match.npre+1)):T0]
        }            
    } else {
        if (is.null(Y.match.time)==FALSE) {
            Y.match.time <- intersect(Tpre, Y.match.time)                
        } else {
            Y.match.time <- Tpre                
        }
    } 



    ## remove other treated (in case of multiple timing) 
    data <- data[c(id.tr, id.co),]
    id.tr <- 1:Ntr
    id.co <- (Ntr + 1): N           


    if (demean == TRUE) { 
        Y.dm.var <- paste0(Y,".dm",Ttot)
        Ypre.mean <- apply(data[, paste0(Y, Tpre), drop = FALSE], 1, mean) # N*1
        outcome.dm <- data[, paste0(Y, Ttot), drop = FALSE] - matrix(Ypre.mean, N, TT) # N * TT
        colnames(outcome.dm) <- Y.dm.var
        data <- cbind.data.frame(data, outcome.dm) 
        Y.match <- paste0(Y,".dm",Y.match.time)
        Y.target <- paste0(Y,".dm",Ttot)
        Y.target.pst <- paste0(Y,".dm",Tpst)
    } else {
        Y.match <- paste0(Y, Y.match.time)
        Y.target <- paste0(Y, Ttot)
        Y.target.pst <- paste0(Y, Tpst)
    }
    if (is.null(Y.match.time) == TRUE) { # do not match on pre-treatment Y
        Y.match <- NULL
    } 
    matchvar <- c(Y.match, X)

    ## tuning parameter
    if (is.null(sigma)==TRUE){
        sigma <- 2 ## sigmas for testing
    } 
    b <- length(matchvar) * sigma


    ## default weights (control first, then treated)
    weights.tr <- rep(1/Ntr, Ntr) 
    weights.co <- rep(1/Nco, Nco)
    w <- c(weights.tr, weights.co*(-1))
    att <- apply(data[, Y.target] * w, 2, sum); names(att) <- Ttot
    att.avg <- mean(att[(T0+1):TT])
    pre.mae.org <- pre.mae <- mean(abs(att[1:T0]))
    
    ## trajectory balancing
    if (is.null(matchvar)==FALSE) { ## have something to balance on

        if (info == TRUE) {
            cat("Seek balance on:\n")
            cat(paste(matchvar, collapse = ", "),"\n")
            cat("\nOptimization:\n")
        }

        # mean balancing
        if (estimator == "mean") {            
            bal.type <- "mbal"
            tmp <- capture.output(
                kbal.out <- suppressWarnings(kbal(allx = data[,matchvar],
                    treatment = data$treat, b=NULL, 
                    linkernel = TRUE, 
                    printprogress = FALSE, sampledinpop = FALSE))
                , file = NULL) 

        } # end of mean balancing

        # kernel balancing
        if (estimator == "kernel") {
            bal.type <- "kbal"
            tmp <- capture.output(
                kbal.out <- suppressWarnings(kbal(allx = data[,matchvar],
                    treatment = data$treat, b=b, 
                    linkernel = FALSE, 
                    printprogress = FALSE, sampledinpop = FALSE))
                , file = NULL)
        } # end of kernel balancing

        # mean first
        if (estimator == "meanfirst") {
            tmp <- capture.output(
                mbal.out <- suppressWarnings(kbal(allx = data[,matchvar],
                    treatment = data$treat, b=NULL, 
                    linkernel = TRUE, 
                    printprogress = FALSE, sampledinpop = FALSE))
                , file = NULL) 
            if (is.null(mbal.out$numdims) == FALSE) {
                mbal.ndims <- mbal.out$numdims # mbal dimensions   
                mbal.svd.keep <- mbal.out$svdK$u[, 1:mbal.ndims] # mbal constraints
                mbal.bias.ratio <- mbal.out$biasbound.opt/mbal.out$biasbound.orig
                success <- TRUE                  
            } else {
                success <- FALSE                  
            }
            # output mbal result        
            if (success == FALSE) {  # mbal failed, no need to do kbal
                kbal.out <- mbal.out                                     
            } else {
                tmp <- capture.output(
                    kbal.out <- suppressWarnings(kbal(allx = data[,matchvar],
                        treatment = data$treat, b=b, 
                        constraint = mbal.svd.keep,
                        printprogress = FALSE, 
                        sampledinpop = FALSE))
                    , file = NULL)
                if (is.null(kbal.out$numdims) == FALSE) {
                    bal.type <- "kbal"
                } else {
                    bal.type <- "mbal" 
                    kbal.out <- mbal.out
                } 
            }
        } # end of mean first

        if (is.null(kbal.out$numdims) == FALSE) {
            success <- TRUE
        }

        ## if success
        if (success == TRUE) {
            weights.tr <- rep(1/Ntr, Ntr) # treated add up to 1; 
            weights.co <- kbal.out$w[id.co]/Nco # controls add up to 1;
            w <- c(weights.tr, weights.co*(-1))
            # reporting
            ndims <- kbal.out$numdims
            bias.ratio <- kbal.out$biasbound.opt/kbal.out$biasbound.orig
            if (bal.type == "mbal") {
                cat(paste0("bias.ratio = ", sprintf("%.4f",bias.ratio),
                    "; num.dims = ",ndims," (mbal)\n"))
            }            
            if (bal.type == "kbal" && estimator == "kernel") {
                cat(paste0("bias.ratio = ", sprintf("%.4f",bias.ratio),
                    "; num.dims = ",ndims," (kbal)\n"))                
            }
            if (bal.type == "kbal" && estimator == "meanfirst") {
                cat(paste0("bias.ratio = ", sprintf("%.4f",bias.ratio),
                    "; num.dims = ",mbal.ndims," + ",ndims," (mbal + kbal)\n"))
            }            
        } else {
            cat("\nSolution not found. Equal weights being used.\n")
        }
    } else {
        bal.type <- "none"    
    }

    # ATT
    att <- apply(data[, Y.target] * w, 2, sum)
    names(att) <- Ttot
    att.avg <- mean(att[(T0+1):TT])   

    
    ## treated and control data
    Y.var  <- paste0(Y, Ttot)
    Y.tr.bar <- apply(data[id.tr, Y.var], 2, mean, na.rm=TRUE)
    Y.co.bar <- apply(data[id.co, Y.var], 2, mean, na.rm=TRUE)
    Y.ct.bar <- Y.tr.bar - att
    Y.bar <- cbind(Y.tr.bar,Y.ct.bar, Y.co.bar)  

    
    #####################
    ## balance table
    #####################
    
    if (is.null(matchvar)==FALSE) {
        
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
    }

    
    #####################
    ## Save Results
    #####################

    out <- list(data.wide = data, # treated units first, then controls
            w = w, # N * 1 vector
            weights.co = weights.co, # Nco * 1 vector
            matchvar = matchvar,
            Y.var = Y.var,
            Y.target = Y.target,
            Y.target.pst = Y.target.pst,
            Y.bar = Y.bar, 
            att = att,
            att.avg = att.avg,
            ntreated = rep(Ntr,TT),
            bal.type = bal.type
            )
    
    if (is.null(matchvar)==FALSE) {
        out <- c(out, list(success = success))
        if (success == TRUE) {
            out <- c(out,
            list(bias.ratio = bias.ratio,
                ndims = ndims,            
                kbal.out = kbal.out,
                bal.table = bal.table,                
                b = b))
            if (estimator == "meanfirst") {
                out <- c(out, list(constraint = mbal.svd.keep))
            }
            if (bal.type == "kbal") {
                out <- c(out, list(b = b))
            }  
        } 

    }             
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
    
    if (is.null(x$est.att.avg) == TRUE) { # no uncertainties
        cat("\n   ~ by Period (including Pre-treatment Periods):\n")
        print(x$att, digits = 4)
        cat("\nAverage Treatment Effect on the Treated:\n")
        print(x$att.avg, digits = 4)
        cat("\nUncertainty estimates not available.\n")
    } else {
        cat("\n   ~ by Period (including Pre-treatment Periods):\n")
        print(x$est.att, digits = 4)        
        cat("\nAverage Treatment Effect on the Treated:\n")
        print(x$est.att.avg, digits = 4)       
    }
}


##########
## Plot
##########

plot.tjbal <- function(x,  
    type = "gap", # c("gap","counterfactual","weights","balance")
    subgroup = NULL, # subgroup plot when sameT0 == FALSE 
    xlim = NULL, 
    ylim = NULL,
    xlab = NULL, 
    ylab = NULL,
    count = TRUE,
    legendOff = FALSE,
    legend.pos = "bottom",
    legend.ncol = NULL,
    legend.labs = NULL,
    main = NULL,
    raw = "none",
    stat = "mean",
    trim = TRUE, ## trim control group in ct plot
    trim.wtot = 0.9, ## show controls whose weights sum up to a number
    theme.bw = TRUE, ## black/white or gray theme
    cex.main = NULL,
    cex.axis = NULL,
    cex.lab = NULL, 
    cex.legend = NULL,
    cex.text = NULL,
    log.weights = FALSE,
    wmin = NULL, ## minimal log weights
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
        if (x$sameT0==FALSE && is.null(subgroup)==TRUE) {
            if (raw %in% c("band","all")) {
                cat("Control units not shown due to different treatment timing among treated.")
            }
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
    if (is.null(subgroup)==FALSE) {
        if (is.numeric(subgroup)==FALSE) {
            stop("\"subgroup\" is not numeric.")
        }
        if (x$sameT0==TRUE) {
            warning("Treatment starts at the same time")
        } else {
            if (!subgroup %in% 1:nrow(x$group.stats)) {
                stop(paste0("There exist ",nrow(x$group.stats)," subgroups."))
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
        x.h <- 0.5
    }

    ## color of axes
    if (theme.bw == TRUE) {
      line.color <- "#AAAAAA70"
    } else {
      line.color <- "white"
    }   

    #### font size
    ## title
    if (is.null(cex.main)==FALSE) {
        if (is.numeric(cex.main)==FALSE) {
            stop("\"cex.main\" is not numeric.")
        }
        cex.main <- 18 * cex.main
    } else {
        cex.main <- 18
    }
    ## axis label
    if (is.null(cex.lab)==FALSE) {
        if (is.numeric(cex.lab)==FALSE) {
            stop("\"cex.lab\" is not numeric.")
        }
        cex.lab <- 13 * cex.lab
    } else {
        cex.lab <- 13
    }
    ## axis number
    if (is.null(cex.axis)==FALSE) {
        if (is.numeric(cex.axis)==FALSE) {
            stop("\"cex.axis\" is not numeric.")
        }
        cex.axis <- 12 * cex.axis
    }  else {
        cex.axis <- 12
    }
    ## legend
    if (is.null(cex.legend)==FALSE) {
        if (is.numeric(cex.legend)==FALSE) {
            stop("\"cex.legend\" is not numeric.")
        }
        cex.legend <- 12 * cex.legend
    }  else {
        cex.legend <- 12
    }
    ## text
    if (is.null(cex.text)==FALSE) {
        if (is.numeric(cex.text)==FALSE) {
            stop("\"cex.text\" is not numeric.")
        }
        cex.text <- 6 * cex.text
    }  else {
        cex.text <- 6
    }
    
    ##-------------------------------##
    ## Take in data
    ##-------------------------------##  

    sameT0 <- x$sameT0
    data <- x$data.wide
    Y.var <- x$Y.var
    id.co <- x$id.co
    time <- x$Ttot
    TT <- length(x$Ttot)
    Nco <- x$Nco
    

    if (is.null(subgroup)==TRUE) {
        att <- x$att
        tb <- x$est.att
        TT <- length(x$Ttot)
        T0 <- x$T0 
        Ntr <- x$Ntr
        N <- x$N 
        w.co <- x$weights.co        
        id.tr <- x$id.tr
        Y.tr <- data[id.tr, Y.var, drop = FALSE]
        Y.co <- data[id.co, Y.var, drop = FALSE]
        sameT0 <- x$sameT0
        Yb <- x$Y.bar[,1:2] ## treated average and counterfactual average  
        if (sameT0 == TRUE) {
            matchvar <- x$matchvar
            bal <- x$bal.table
            ntreated <- rep(Ntr, length(time))
        } else {
            ntreated <- x$ntreated
        } 

    } else {  # make plots for subgroup
        att <- x$sub.att[,subgroup]      
        tb <- x$est.sub.att[,,subgroup]
        T0 <- x$group.stats[subgroup,"T0"] 
        Ntr <- x$group.stats[subgroup,"Ntr"] 
        N <- Ntr + Nco        
        w.co <- x$sub.weights.co[,subgroup]
        id.tr <- x$id.tr[which(x$T0.tr == T0)] 
        Y.tr <- data[id.tr, Y.var, drop = FALSE]
        Y.co <- data[id.co, Y.var, drop = FALSE]
        Y.tr.bar <- apply(data[id.tr, Y.var], 2, mean, na.rm=TRUE)
        Y.ct.bar <- Y.tr.bar - att
        Yb <- cbind(Y.tr.bar,Y.ct.bar) ## treated average and counterfactual average  
        sameT0 <- TRUE
        matchvar <- x$matchvar.list[[subgroup]]       
        bal <- x$bal.table.list[[subgroup]]  
        ntreated <- rep(Ntr, length(time))
    }


    ##-------------------------------##
    ## Plotting
    ##-------------------------------##  
    
   
    ## parameters
    line.width <- c(1.2,0.5)
  
    ## type of plots
    if(type %in% c("gap","counterfactual") && sameT0 == TRUE) {
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

    if (type %in% c("gap","counterfactual") && sameT0 == FALSE)  { ## variable treatment timing
         ## recenter based on treatment timing
        time <- c(-(max(T0) -1) : (TT-min(T0)))
        TT <- length(time)
        time.bf <- 0 ## before treatment
        if (length(xlim) != 0) {
            show <- which(time>=xlim[1]& time<=xlim[2])     
        } else {
            show <- 1:length(time)    
            Ntr.max <- max(ntreated)
            show <- show[which(ntreated >= Ntr.max/3)]
        }        
    }

    if (type %in% c("gap","counterfactual")) {
        nT <- length(show)
        time.label <- time[show]
        T.b <- 1:length(show)
        hist.adjust <- (max(time[show])-min(time[show]))/(nT-1) # adjust for the width of the histograms
    }

    ## legend on/off
    if (legendOff == TRUE) {
        legend.pos <- "none"
    }

    if (is.null(legend.ncol)==FALSE) {
        if (is.numeric(legend.ncol)==FALSE) {
            stop("\n\"legend.ncol\" needs to be numeric.")
        }
    } 

    if (is.null(legend.labs)==FALSE) {
        legend.labs <- as.character(legend.labs)
    }

    ## confidence intervals
    if (is.null(tb)==TRUE) {
        CI <- FALSE
    } else {
        CI <- TRUE
        if (Ntr == 1) { # for subgroup
            CI <- FALSE
        }
    }

    ############  START  ###############
    
    if (type == "gap") { 

        max.count.pos <- time[min(intersect(show,which(time > time.bf)))]
        
        if (Ntr>1) {
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
        if (CI==FALSE) { 
            cat("Uncertainty estimates not available.\n")
            data <- cbind.data.frame(time, ATT = att, n.Treated = ntreated)[show,]             
        } else {
            data <- cbind.data.frame(time, tb)[show,]
        }        

        # height of the histogram
        if (CI == FALSE) {
            if (length(ylim) != 0) {
                rect.length <- (ylim[2] - ylim[1]) / 5
                rect.min <- ylim[1]
            } else {
                rect.length <- (max(data[,"ATT"]) - min(data[,"ATT"]))/4
                rect.min <- min(data[,"ATT"]) - rect.length
            } 
        } else {
            if (length(ylim) != 0) {
                rect.length <- (ylim[2] - ylim[1]) / 5
                rect.min <- ylim[1]
            } else {
                rect.length <- (max(data[,"CI.upper"]) - min(data[,"CI.lower"]))/4
                rect.min <- min(data[,"CI.lower"]) - rect.length
            }  
        }

        ## plotting
        p <- ggplot(data)
        if (theme.bw == TRUE) {
          p <- p + theme_bw()
        }
        p <- p +
        geom_vline(xintercept = time.bf, colour=line.color,size = 2) +
        geom_hline(yintercept = 0, colour=line.color,size = 2) +
        xlab(xlab) +  ylab(ylab) 

        ## histogram
        if (count == TRUE) {
            data[,"xmin"] <- data[,"time"] - 0.2 * hist.adjust
            data[,"xmax"] <- data[,"time"] + 0.2 * hist.adjust
            data[,"ymin"] <- rep(rect.min, nrow(data))
            data[,"ymax"] <- rect.min + (data$n.Treated/Ntr) * 0.8 * rect.length
            xx <- range(data$time)
            p <- p + geom_rect(data = data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                fill = "grey70", colour = "grey69", alpha = 0.4, size = 0.2)
            p <- p + annotate("text", x = max.count.pos + 0.00 * (xx[2]-xx[1]), 
                y = max(data$ymax) + 0.2 * rect.length, 
                label = paste("Ntr =",Ntr), size = cex.text * 0.8, hjust = 0.5)
        }

        ## point estimates
        p <- p + geom_line(aes(time, ATT), size = 1.2)

        ## confidence intervals
        if (CI == TRUE) {
            p <- p + geom_ribbon(aes(x = time, ymin=CI.lower, ymax=CI.upper),alpha=0.2)
        }

        
    } else if (type=="counterfactual") { 

        if (sameT0 == TRUE) {            
            if (length(matchvar)==0) { ## no balancing
                trim <- FALSE
            }
            if (trim == TRUE & raw %in% c("band","all")) {
                Nco <-sum(1 - (cumsum(sort(w.co,decreasing = TRUE))>trim.wtot)) + 1  # how many control units left
                trim.id<- order(w.co, decreasing = TRUE)[1:Nco]
                Y.co <- Y.co[trim.id, ,drop = FALSE]
                id.co <- id.co[trim.id]
                co.label <- paste0("Heavily Weighted Controls (",floor(trim.wtot*100),"% weights)")
                N <- Ntr + Nco
            } else {
                co.label <- "Controls"
            }
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

            data <- cbind.data.frame(
                "time" = rep(time[show],2),
                "outcome" = c(Yb[show,1],Yb[show,2]),
                "type" = c(rep("tr",nT),rep("co",nT)),
                "n.Treated" = rep(ntreated[show],2)) 

            # height of the histogram
            max.count.pos <- time[min(intersect(show,which(time > time.bf)))]
            if (length(ylim) != 0) {
                rect.length <- (ylim[2] - ylim[1]) / 5
                rect.min <- ylim[1]
            } else {
                rect.length <- (max(data[,"outcome"]) - min(data[,"outcome"]))/4
                rect.min <- min(data[,"outcome"]) - rect.length
            } 

            ## theme 
            p <- ggplot(data)
            if (theme.bw == TRUE) {
              p <- p + theme_bw()
            }
            p <- p + xlab(xlab) +  ylab(ylab) +
            geom_vline(xintercept=time.bf,colour=line.color,size = 2) 

            ## histogram
            if (count == TRUE) {
                data[,"xmin"] <- data[,"time"] - 0.2 * hist.adjust
                data[,"xmax"] <- data[,"time"] + 0.2 * hist.adjust
                data[,"ymin"] <- rep(rect.min, nrow(data))
                data[,"ymax"] <- rect.min + (data$n.Treated/Ntr) * 0.8 * rect.length
                xx <- range(data$time)
                p <- p + geom_rect(data = data[which(data$type == "tr"),], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                    fill = "grey70", colour = "grey69", alpha = 0.4, size = 0.2)
                p <- p + annotate("text", x = max.count.pos + 0.00 * (xx[2]-xx[1]), 
                    y = max(data$ymax) + 0.2 * rect.length, 
                    label = paste("Ntr =",Ntr), size = cex.text * 0.8, hjust = 0.5)
            }
                
            ## main
            p <- p + geom_line(data = data, aes(time, outcome,colour = type,size = type,linetype = type))

             
            ## legend
            if (is.null(legend.labs)==FALSE) {
                if (length(legend.labs)!=2) {
                    stop("\n\"legend.labs\" needs to be of length 2.")
                }
                set.labels <- legend.labs
            } else {
                if (Ntr == 1) {
                    set.labels = c("Treated", "Estimated Y(0)")
                } else {
                    set.labels = c("Treated Average", "Estimated Y(0) Average")
                }
            }
            set.limits = c("tr","co")            
            set.colors = c("black","steelblue")
            set.linetypes = c("solid","longdash")
            set.linewidth = rep(line.width[1],2)
            if (is.null(legend.ncol)==TRUE) {
                legend.ncol <- 2
            }
                

        } else if  (raw == "band") {

            if (trim == FALSE) {
                Y.tr.band <- t(apply(Y.tr, 2, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                Y.co.band <- t(apply(Y.co, 2, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                tr.band.label <- "Treated 5-95% Quantiles"
                co.band.label <- "Controls 5-95% Quantiles"
            } else {
                Y.tr.band <- t(apply(Y.tr, 2, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                Y.co.band <- t(apply(Y.co, 2, quantile, prob=c(0.05,0.95),na.rm=TRUE))
                tr.band.label <- "Treated 5-95% Quantiles"
                co.band.label <- paste0("Heavily Weighted Controls 5-95% Quantiles (",floor(trim.wtot*100),"% weights)")
            }

                
            data <- cbind.data.frame(
                "time" = rep(time[show],2),
                "outcome" = c(Yb[show,1], Yb[show,2]),
                "type" = c(rep("tr",nT),rep("co",nT)),
                "n.Treated" = rep(ntreated[show],2))

            
            if (Ntr <= 10) {
                data.band <- cbind.data.frame(time, Y.co.band)[show,]
                colnames(data.band) <- c("time","co.lower","co.upper")
                range.max <- max(data.band[,c("co.upper")])
                range.min <- min(data.band[,c("co.lower")])
            } else {
                data.band <- cbind.data.frame(time, Y.tr.band, Y.co.band)[show,]
                colnames(data.band) <- c("time","tr.lower","tr.upper","co.lower","co.upper")
                range.max <- max(data.band[,c("co.upper","tr.upper")])
                range.min <- min(data.band[,c("co.lower","tr.lower")])
            }


            # height of the histogram
            max.count.pos <- time[min(intersect(show,which(time > time.bf)))]
            if (length(ylim) != 0) {
                rect.length <- (ylim[2] - ylim[1]) / 5
                rect.min <- ylim[1]
            } else {
                rect.length <- (range.max - range.min)/4
                rect.min <- range.min - rect.length
            } 

            ## theme 
            p <- ggplot(data)
            if (theme.bw == TRUE) {
              p <- p + theme_bw()
            }
            p <- p + xlab(xlab) +  ylab(ylab) +
            geom_vline(xintercept=time.bf,colour=line.color,size = 2)

            ## histogram
            if (count == TRUE) {
                data[,"xmin"] <- data[,"time"] - 0.2 * hist.adjust
                data[,"xmax"] <- data[,"time"] + 0.2 * hist.adjust
                data[,"ymin"] <- rep(rect.min, nrow(data))
                data[,"ymax"] <- rect.min + (data$n.Treated/Ntr) * 0.8 * rect.length
                xx <- range(data$time)
                p <- p + geom_rect(data = data[which(data$type == "tr"),], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                    fill = "grey70", colour = "grey69", alpha = 0.4, size = 0.2)
                p <- p + annotate("text", x = max.count.pos + 0.00 * (xx[2]-xx[1]), 
                    y = max(data$ymax) + 0.2 * rect.length, 
                    label = paste("Ntr =",Ntr), size = cex.text * 0.8, hjust = 0.5)
            }

            ## main
            p <- p + geom_line(aes(time, outcome,colour = type,size = type,linetype = type))
            ## band
            if (Ntr <= 10) {
                p <- p + geom_ribbon(data = data.band,
                    aes(ymin = co.lower, ymax = co.upper, x=time),
                    alpha = 0.15, fill = "steelblue")
                set.limits = c("tr","co","co.band")
                # legend
                if (is.null(legend.labs)==FALSE) {
                    if (length(legend.labs)!=3) {
                        stop("\n\"legend.labs\" needs to be of length 3.")
                    }
                    set.labels <- legend.labs
                } else {
                    if (Ntr == 1) {
                        set.labels = c("Treated", "Estimated Y(0)",co.band.label)
                    } else {
                        set.labels = c("Treated Average", "Estimated Y(0)",co.band.label)
                    }
                }                 
                set.colors = c("black","steelblue","#4682B470")
                set.linetypes = c("solid","longdash","solid")
                set.linewidth = c(rep(line.width[1],2),4)  
                # legend columns
                if (is.null(legend.ncol)==TRUE) {
                    legend.ncol <- 1
                }


            } else {
                p <- p + geom_ribbon(data = data.band,
                    aes(ymin = co.lower, ymax = co.upper, x=time),
                    alpha = 0.15, fill = "steelblue") +
                geom_ribbon(data = data.band,
                   aes(ymin = tr.lower, ymax = tr.upper, x=time),
                   alpha = 0.15, fill = "black")
                # legend
                if (is.null(legend.labs)==FALSE) {
                    if (length(legend.labs)!=4) {
                        stop("\n\"legend.labs\" needs to be of length 4.")
                    }
                    set.labels <- legend.labs
                } else {
                    set.labels <- c("Treated Average", "Estimated Y(0) Average", tr.band.label, co.band.label)
                } 
                set.limits = c("tr","co","tr.band","co.band")
                set.colors = c("black","steelblue","#77777750","#4682B470")
                set.linetypes = c("solid","longdash","solid","solid")
                set.linewidth = c(rep(line.width[1],2),4,4)
                # legend columns
                if (is.null(legend.ncol)==TRUE) {
                    legend.ncol <- 2
                }        
            } 

                           
          

        } else if (raw == "all") { ## plot all the raw data

            if (Ntr == 1) {
                data <- cbind.data.frame("time" = rep(time[show],(2 + Nco)),
                   "outcome" = c(Yb[show,1],Yb[show,2], c(t(Y.co[,show]))),
                   "type" = c(rep("tr",nT), rep("co",nT), rep("raw.co",(Nco * nT))),
                   "id" = c(rep("tr",nT), rep("co",nT), rep(id.co, each = nT)),
                   "n.Treated" = rep(ntreated[show],(2 + Nco)))                
            } else {
                data <- cbind.data.frame("time" = rep(time[show],(2 + N)),
                   "outcome" = c(Yb[show,1],Yb[show,2],c(t(Y.tr[,show])),c(t(Y.co[,show]))),
                   "type" = c(rep("tr",nT), rep("co",nT), rep("raw.tr",(Ntr * nT)), rep("raw.co",(Nco * nT))),
                   "id" = c(rep("tr",nT), rep("co",nT), rep(c(id.tr,id.co), each = nT)),
                   "n.Treated" = rep(ntreated[show],(2+N))) 
            }

            # height of the histogram
            range.max <- max(data[,c("outcome")])
            range.min <- min(data[,c("outcome")]) 
            max.count.pos <- time[min(intersect(show,which(time > time.bf)))]
            if (length(ylim) != 0) {
                rect.length <- (ylim[2] - ylim[1]) / 5
                rect.min <- ylim[1]
            } else {
                rect.length <- (range.max - range.min)/4
                rect.min <- range.min - rect.length
            } 

                
            ## theme
            p <- ggplot(data)
            if (theme.bw == TRUE) {
              p <- p + theme_bw()
            }
            p <- p + xlab(xlab) +  ylab(ylab) + geom_vline(xintercept=time.bf,colour=line.color,size = 2) 

            ## histogram
            if (count == TRUE) {
                data[,"xmin"] <- data[,"time"] - 0.2 * hist.adjust
                data[,"xmax"] <- data[,"time"] + 0.2 * hist.adjust
                data[,"ymin"] <- rep(rect.min, nrow(data))
                data[,"ymax"] <- rect.min + (data$n.Treated/Ntr) * 0.8 * rect.length
                xx <- range(data$time)
                p <- p + geom_rect(data = data[which(data$type == "tr"),], aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                    fill = "grey70", colour = "grey69", alpha = 0.4, size = 0.2)
                p <- p + annotate("text", x = max.count.pos + 0.00 * (xx[2]-xx[1]), 
                    y = max(data$ymax) + 0.2 * rect.length, 
                    label = paste("Ntr =",Ntr), size = cex.text * 0.8, hjust = 0.5)
            }

            ## main
            p <- p + geom_line(aes(time, outcome, colour = type,size = type,
               linetype = type, group = id))

                    ## legend
            if (Nco >= 50) {
                co.color = "#4682B430"
            } else if (Nco >= 20) {
                co.color = "#4682B450"
            } else {
                co.color = "#4682B470"
            }     
            if (Ntr == 1) {
                # legend
                if (is.null(legend.labs)==FALSE) {
                    if (length(legend.labs)!=3) {
                        stop("\n\"legend.labs\" needs to be of length 3.")
                    }
                    set.labels <- legend.labs
                } else {
                    set.labels <- c("Treated","Estimated Treated Y(0)",co.label) 
                }       
                set.limits = c("tr","co","raw.co")         
                set.colors = c("black","steelblue",co.color)
                set.linetypes = c("solid","longdash","solid")
                set.linewidth = line.width[c(1,1,2)]
                # legend columns
                if (is.null(legend.ncol)==TRUE) {
                    legend.ncol <- 1
                }  
            } else {
                # legend
                if (is.null(legend.labs)==FALSE) {
                    if (length(legend.labs)!=4) {
                        stop("\n\"legend.labs\" needs to be of length 4.")
                    }
                    set.labels <- legend.labs
                } else {
                    set.labels <- c("Treated Average","Estimated Y(0) Average","Treated", co.label)
                }
                set.limits = c("tr","co","raw.tr","raw.co")
                set.colors = c("black","steelblue","#77777750",co.color)
                set.linetypes = c("solid","longdash","solid","solid")
                set.linewidth = rep(line.width,each=2)
                # legend columns
                if (is.null(legend.ncol)==TRUE) {
                    legend.ncol <- 2
                }  
            }
            
        }

        # legend 
        p <- p + scale_colour_manual(limits = set.limits,
               labels = set.labels,
               values =set.colors) +
            scale_linetype_manual(limits = set.limits,
              labels = set.labels,
              values = set.linetypes) +
            scale_size_manual(limits = set.limits,
              labels = set.labels,
              values = set.linewidth) 
        p <- p + guides(linetype = guide_legend(title=NULL, ncol=legend.ncol, bycol=TRUE),
             colour = guide_legend(title=NULL, ncol=legend.ncol, bycol=TRUE),
             size = guide_legend(title=NULL, ncol=legend.ncol), bycol=TRUE) 

        # x axis labels
        if (!is.numeric(time.label)) {
            p <- p + scale_x_continuous(expand = c(0, 0), breaks = show[T.b], labels = time.label[T.b])
        }

        # key width in legend
        p <- p + theme(legend.key.width = unit(2.5,"line"))   


        

    } else if (type == "weights") {

        if (sameT0==FALSE) {
            cat("Weighted by the number of treated unit in each subgroup.")
        }

        maintext <- "Weights of the Control Units"
        w <- w.co
        if (log.weights == TRUE) {
            w <- log(w)
        }
        if (is.null(wmin)==TRUE) {
            if (log.weights==TRUE) {
                wmin <- -10
            } else {
                wmin <- 0
            }
        }
        w[w < wmin] <- wmin
        if (is.null(xlab)==TRUE) {
            if (log.weights==TRUE) {
                xlab <- "Weight (logarithmic)"                
            } else {
                xlab <- "Weight"
            }            
        }
        if (is.null(ylab)==TRUE) {
            ylab <- "Counts"
        }        
        w.min <- min(w)
        w.max <- max(w)
        if (w.min!=w.max) {
            bw <- abs(w.max - w.min)/30
        } else {
            bw <- abs(w.min)/20            
        }
        p <- qplot(w, col = I("gray70"), binwidth = bw, xlab = xlab, ylab = ylab)
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

        if (sameT0==FALSE) {
            stop("Requires treatment starts at the same time or specifying a subgroup.")
        }

        
        if (!stat %in% c("mean","sd")) {
            stop("Wrong specification for \"stat\".")
        }
        if (Ntr == 1 & stat =="sd") {
            stop("Standard deviation does not exist with one treated unit.")
        }

        if (is.null(matchvar)==TRUE) {
            stop("No covariates being balanced on")            
        }

        if (is.null(subgroup)==TRUE) {
            maintext <- "Covariate Balance"
        } else {
            maintext <- paste0("Covariate Balance (T0 = ",time[T0],")")
        }
        
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

    } # end of weights

     p <- p + theme(legend.text = element_text(margin = margin(r = 10, unit = "pt"), size = cex.legend),
       legend.position = legend.pos,
       legend.background = element_rect(fill="transparent",colour=NA),
       axis.title=element_text(size=cex.lab),
       axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
       axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0)),
       axis.text = element_text(color="black", size=cex.axis),
       axis.text.x = element_text(size = cex.axis, angle = angle, hjust=x.h, vjust=x.v),
       axis.text.y = element_text(size = cex.axis),
       plot.title = element_text(size = cex.main, hjust = 0.5, face="bold", margin = margin(10, 0, 10, 0)))
       



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
