###################################################
#### tjbal.multi: allow multiple treatment timing
###################################################


tjbal.multi <- function(
    data, ## data in wide form
    Y,
    D,
    X,
    cat.columns = NULL, # list of columns which contain categorical one-hot/dummy encoded variables
    mixed.data=FALSE, # TRUE if X contains categorical/one-hot/dummy encoded variables
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
    cores = 4
    ) { 

    
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
    success <- sigmas <- ndims <- ndims.mbal <- rep(0, nT0) 
    bias.ratios <- rep(1, nT0) 
    kernels <- rep(NA, nT0)
    

    cat("Balancing...\n")  
    for (i in 1:nT0) {

        T0 <- T0.unique[i]
        cat(paste0("Subgroup T0 = ",T0,": "))
        id.tr.one <- which(T0.all == T0)    

        out <- tjbal.core(data = data, Y = Y, X = X,
            cat.columns = cat.columns, mixed.data=mixed.data=FALSE,
            Y.match.time = Y.match.time, 
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
                cores <- max(detectCores() - 1, 1)
            }
            para.clusters <- makeClusterPSOCK(cores,verbose = FALSE)
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
                                    K = K[-id,], mixed_data = mixed_data, cat_columns=cat.columns,
                                    constraint = constraint[-id,], 
                                    treatment = data[-id, D],
                                    linkernel = (1-kernel), b = sigma, 
                                    printprogress = FALSE, fullSVD = TRUE,
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
