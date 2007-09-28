## Definition of pallete colors to use
##
## Gustavo H. Esteves
## 18/05/07
##
## Version: 1.0
##


## Defining pallete green -> black -> red
greenRed <- function() {
    
    low <- col2rgb("green")
    mid <- col2rgb("black")
    high <- col2rgb("red")
    lower <- floor(100/2)
    upper <- 100-lower
    
    red <- c(seq(low[1, 1], mid[1, 1], length=lower), seq(mid[1, 1], high[1, 1],
    length=upper))/255
    
    green <- c(seq(low[3, 1], mid[3, 1], length=lower), seq(mid[3, 1], high[3, 1],
    length=upper))/255
    
    blue <- c(seq(low[2, 1], mid[2, 1], length=lower), seq(mid[2, 1], high[2, 1],
    length=upper))/255
    
    
    return(rgb(red, blue, green))

}
    

## Defining pallete  black -> blue
blackBlue <- function() {
    
    low <- col2rgb("black")
    high <- col2rgb("blue")
    lower <- floor(100/2)
    upper <- 100-lower
    
    red <- seq(low[1, 1], high[1, 1], length=100)/255
    green <- seq(low[3, 1], high[3, 1], length=100)/255
    blue <- seq(low[2, 1], high[2, 1], length=100)/255
    
    return(rgb(red, blue, green))

}
