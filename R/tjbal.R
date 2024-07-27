
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
    estimator = "mean",  # mean, meanfirst, kernel
    sigma=NULL,
    print.baltable = TRUE, # print out table table
    vce = "jackknife", ## uncertainty estimates
    conf.lvl = 0.95, ## confidence interval
    nsims = 200, ## number of bootstrap runs
    parallel = TRUE, ## parallel computing
    cores = NULL,
    seed = NULL
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
    seed = NULL
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
    seed = NULL
    ) {

    ##-------------------------------##
    ## Checking Parameters
    ##-------------------------------##  

    if (is.null(seed)==FALSE) {set.seed(seed)}

    
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
            nsims = nsims, parallel = parallel, cores = cores)         
    } else { 
        bal.out <- tjbal.multi(data = data.wide, Y = Yname, D = "treat", X = Xname,
            Y.match.time = Y.match.time, Y.match.npre = Y.match.npre, 
            Ttot = Ttot, unit = "id", 
            demean = demean, estimator = estimator, sigma = sigma, 
            vce = vce, conf.lvl = conf.lvl,
            nsims = nsims, parallel = parallel, cores = cores)  
    } 

    out <- c(list(sameT0 = sameT0, index = index, Yname = Yname), bal.out)
    out$call <- match.call()
    class(out) <- "tjbal"
    return(out)
}







