## Function to do contrasts fit without doing eBayes. Function very similar with
## contrasts.fit from limma
##
## Parameters: The same of contrasts.fit from limma
##
## Gustavo H. Esteves (addapted from contrasts.fit function from limma package)
## 27/05/07
##
##


contrastsFitM <- function (fit, contrasts) {
  
  ncoef <- NCOL(fit$coefficients)
  if (NROW(contrasts) != ncoef) 
    stop("Number of rows of contrast matrix must match number of coefficients")
  fit$contrasts <- contrasts
  cormatrix <- cov2cor(fit$cov.coefficients)
  if (is.null(cormatrix)) {
    warning("no coef correlation matrix found in fit - assuming orthogonal")
    cormatrix <- diag(ncoef)
  }
  r <- nrow(cormatrix)
  if (r < ncoef) {
    if (is.null(fit$pivot)) 
      stop("cor.coef not full rank but pivot column not found in fit")
    est <- fit$pivot[1:r]
    if (any(contrasts[-est, ])) 
      stop("trying to take contrast of non-estimable coefficient")
    contrasts <- contrasts[est, , drop=FALSE]
    fit$coefficients <- fit$coefficients[, est, drop=FALSE]
    fit$stdev.unscaled <- fit$stdev.unscaled[, est, drop=FALSE]
    ncoef <- r
  }
  fit$coefficients <- fit$coefficients %*% contrasts
  if (length(cormatrix) < 2) {
    orthog <- TRUE
  }
  else {
    orthog <- all(abs(cormatrix[lower.tri(cormatrix)]) < 
                  1e-14)
  }
  R <- chol(fit$cov.coefficients)
  fit$cov.coefficients <- crossprod(R %*% contrasts)
  fit$pivot <- NULL
  if (orthog) 
    fit$stdev.unscaled <- sqrt(fit$stdev.unscaled^2 %*% contrasts^2)
  else {
    R <- chol(cormatrix)
    ngenes <- NROW(fit$stdev.unscaled)
    ncont <- NCOL(contrasts)
    U <- matrix(1, ngenes, ncont, dimnames=list(rownames(fit$stdev.unscaled), 
                                    colnames(contrasts)))
    o <- array(1, c(1, ncoef))
    for (i in 1:ngenes) {
      RUC <- R %*% (diag(fit$stdev.unscaled[i, ]) %*% contrasts)
      U[i, ] <- sqrt(o %*% RUC^2)
    }
    fit$stdev.unscaled <- U
  }
  
  ## Adding (conventional) t and F stats
  fit$t <- fit$coefficients/fit$stdev.unscaled/fit$sigma
  fit$p.value <- 2 * pt(-abs(fit$t), df=fit$df.residual)
  
  ## Calculating F
  F.stat <- classifyTestsF(fit$t, df=fit$df.residual, fstat.only=TRUE)
  fit$F <- as.vector(F.stat)
  df1 <- attr(F.stat, "df1")
  df2 <- attr(F.stat, "df2")
  if (df2[1] > 1e+06) 
    fit$F.p.value <- pchisq(df1 * fit$F, df1, lower.tail=FALSE)
  else fit$F.p.value <- pf(fit$F, df1, df2, lower.tail=FALSE)
  
  return(fit)
  
}
