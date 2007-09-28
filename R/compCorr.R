## Function to calculate a p-value (from compcorr function) between to
## correlations values and the sample sizes.
##
## Parameters: n1 -> sample length of 1st correlation
##             r1 -> first correlation value
##             n2 -> sample length of 2nd correlation
##             r2 -> second correlation value
##
## function from web (by Christian Stratowa,
## email: christian.stratowa@vie.boehringer-ingelheim.com)
## see:
##      http://ftp.sas.com/techsup/download/stat/compcorr.html
##      http://www.fon.hum.uva.nl/Service/Statistics/Two_Correlations.html
##
## Gustavo H. Esteves
## 14/05/07
##
## Version: 1.0
##


compCorr <- function(n1, r1, n2, r2) {
    
    
    ## Changing values exactly equal to one for values close to one but not one
    ## this prevents problems at Fisher's tranformation
    r1[r1 == 1] <- 0.9999999999
    r2[r2 == 1] <- 0.9999999999
    
    
    ## Fisher Z-transform
    zf1 <- 0.5*log((1+r1)/(1-r1))
    zf2 <- 0.5*log((1+r2)/(1-r2))
    
    
    ## difference
    dz <- (zf1-zf2)/sqrt(1/(n1-3)+(1/(n2-3)))
    
    
    ## p-value
    pv <- 2*(1-pnorm(abs(dz)))
    
    ## Returning the object
    return(list(diff=dz, pval=pv))
    
}
