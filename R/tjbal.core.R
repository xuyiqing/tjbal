
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
                kbal.out <- suppressWarnings(kbal(allx = data[,matchvar,drop = FALSE],
                    treatment = data$treat, b=NULL, 
                    linkernel = TRUE, 
                    printprogress = FALSE, sampledinpop = FALSE))
                , file = NULL) 

        } # end of mean balancing

        # kernel balancing
        if (estimator == "kernel") {
            bal.type <- "kbal"
            tmp <- capture.output(
                kbal.out <- suppressWarnings(kbal(allx = data[,matchvar,drop = FALSE],
                    treatment = data$treat, b=b, 
                    linkernel = FALSE, 
                    printprogress = FALSE, sampledinpop = FALSE))
                , file = NULL)
        } # end of kernel balancing

        # mean first
        if (estimator == "meanfirst") {
            tmp <- capture.output(
                kbal.out <- suppressWarnings(kbal(allx = data[,matchvar,drop = FALSE],
                    treatment = data$treat, b=b, 
                    linkernel = FALSE, meanfirst = TRUE,
                    printprogress = FALSE, sampledinpop = FALSE))
                , file = NULL) 
            if (is.null(kbal.out$meanfirst_dims) == FALSE) {
                mbal.ndims <- kbal.out$meanfirst_dims # mbal dimensions   
                mbal.svd.keep <- kbal.out$meanfirst_cols # mbal constraints
                kbal.bias.ratio <- kbal.out$biasbound_opt/kbal.out$biasbound_orig
                if (is.null(kbal.out$numdims) == FALSE) {
                    bal.type <- "kbal"
                } else {
                    bal.type <- "mbal" 
                } 
                success <- TRUE                  
            } else {
                success <- FALSE                  
            }
        } # end of mean first

        if (is.null(kbal.out$numdims) == FALSE) {
            success <- TRUE
        }

        ## if success
        if (success == TRUE) {
            weights.tr <- rep(1/Ntr, Ntr) # treated add up to 1; 
            weights.co <- kbal.out$w[id.co]/Nco # controls add up to 1;
            w <- c(weights.tr, weights.co * (-1))
            # reporting
            ndims <- kbal.out$numdims
            bias.ratio <- kbal.out$biasbound_opt/kbal.out$biasbound_orig
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

    

