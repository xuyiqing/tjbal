#######################################################
## S3 METHODS
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
