
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
    cores = 4 
    ) { 

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
        info = TRUE, balance.table = TRUE) 
   
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
    ## kernel type follows the estimator, not bal.type: meanfirst with no
    ## kernel dimensions reports "mbal" but its K is still a Gaussian kernel
    kernel <- as.numeric(estimator != "mean")
    reweighted <- isTRUE(out$success) # FALSE: main fit used equal weights
    data <- out$data.wide
    id.tr <- 1:Ntr
    id.co <- (Ntr + 1): N

    ## control weights for a resample (controls first in sample.id); controls
    ## add up to -1. Re-balances with the dimensions the main fit chose; if the
    ## main fit fell back to equal weights, so does the resample.
    resample.weights <- function(sample.id) {
        if (is.null(matchvar) == TRUE || reweighted == FALSE) {
            return(rep(-1/Nco, Nco))
        }
        data.tmp <- data[sample.id, ]
        args <- list(allx = data.tmp[, matchvar], treatment = data.tmp[, D],
            K = K[sample.id, , drop = FALSE],
            linkernel = (1-kernel), fullSVD = TRUE,
            printprogress = FALSE, sampledinpop = FALSE)
        if (is.null(constraint) == FALSE) {
            args$constraint <- constraint[sample.id, , drop = FALSE]
        }
        if (ndims == 0) {
            if (is.null(constraint) == TRUE) {return(rep(-1/Nco, Nco))}
            ## meanfirst kept no kernel dimensions: balance on the mean-first
            ## columns only (kbal refuses maxnumdims = 0 and numdims = 0)
            treat <- data.tmp[, D]
            getw.out <- kbal::getw(target = treat, observed = 1 - treat,
                svd.U = constraint[sample.id, , drop = FALSE])
            return(getw.out$w[1:Nco]/Nco*(-1))
        }
        args$minnumdims <- max(0,ndims-5)
        args$maxnumdims <- ndims
        tmp <- capture.output(
            kbal.re <- suppressWarnings(do.call(kbal, args))
        , file = NULL)
        return(kbal.re$w[1:Nco]/Nco*(-1))
    }

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
            w.boot[1:Nco] <- resample.weights(sample.id)
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
            w.jack[1:Nco] <- resample.weights(sample.id)
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
