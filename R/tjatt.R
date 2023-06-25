
###################################################
#### att: obtaining ATT estimates
###################################################

att <- function(
	data, # wide format
    Y.target, # outcome variable names (including pres)
    X, # covariates
    w, # weights
    cl, # clustering variable
	method = "lm_robust",
	dr = TRUE,
	...
	){
	elpss <- list(...)
	if (dr == FALSE & method == "lin") {
		method <- "lm_robust"
	}
	dat <- data[, unique(c("unit", "treat", Y.target, cl))]
    nY <- length(Y.target)
    att <- matrix(NA, nrow = nY, ncol = 2)
    colnames(att) <- c("att", "se")
    ## Estimate ATT for each period
    for (i in 1:nY) {
        outcome <- Y.target[i]
        s <- dat[, c(outcome, cl)]
        if (dr == TRUE) {
            fml <- as.formula(paste(outcome, "~ treat + ", paste(X, collapse = " + ")))
        } else {
            fml <- as.formula(paste(outcome, "~ treat"))
        }
        fml <- as.formula(paste(outcome, "~ treat"))
        if (method == "lm_robust") {
            out <- lm_robust(fml, data = dat, weights = w, 
                se_type = "stata", cluster = dat[, cl])
        }
        if (method == "lin") {
            out <- lin(fml, data = dat, weights = w, 
                se_type = "stata", cluster = dat[, cl])
        }
        att[i, ] <- as.numeric(c(out$coefficients["treat"],out$std.error["treat"]))
    }
    ## Estimate ATT for all periods
    # step 1: tranforms form data to long format
    # step 2: estimate using estimatr

    # unfinished
   
}




