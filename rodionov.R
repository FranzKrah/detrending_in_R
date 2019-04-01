## Implementation of Rodionov 2004 regime-shift detection
#' @param timeseries numeric vector
#' @param L cut-off length of the regimes
#' @param p probability level
#' @author Franz-Sebastian Krah
#' @references S. N. Rodionov (2004). A sequential algorithm for testing climate regime shifts. Geophysical Research Letters 31(9):L09204 
#' @import foreach
#' @details Still a test version. Not all L values work. 


RSI <- function(timeseries, regime.start, L, diff){
  
  rs <- regime.start
  x_star <- vector()
  R_new <- timeseries[rs]
  
  i_end <- ifelse((rs+L-1) > length(timeseries), length(timeseries), (rs+L-1))
  
  for(i in (rs+1):i_end){
    
    if(timeseries[i] > R_new){
      x_star_new <- timeseries[i] - R_new
    }
    if(timeseries[i] < R_new){
      x_star_new <- R_new - timeseries[i]
    }  
    x_star_new <- x_star_new/(L*var(timeseries[1:L]))
    x_star <- append(x_star, x_star_new)
  }
  return(cumsum(x_star))
}

search.regime.start <- function(timeseries, from, L, up, lo, diff){
  
  i_start <- ifelse((from + L) > length(timeseries), length(timeseries), (from + L))
  
  for(i in i_start:length(timeseries)){
    
    exceed <- (timeseries[i] > up) | (timeseries[i] < lo)
    exceed <- any(exceed)
    
    if(!exceed){ # no change
      
      ## recalc xR1
      R_update <- mean(timeseries[(i-L+1):i])
      up <- R_update + diff
      lo <- R_update - diff
      
    }else{
      R2_j <- i
      print(paste("regime shift found:", R2_j))
      return(R2_j)
      # break
    }
  }
}


rodionov <- function(timeseries = x, L = l, p = p){
  
  require(foreach)
  
  X <- timeseries
  ## step 2
  t <- qt((0+(p/2)), df=(2*l-2), lower.tail = FALSE)
  sigma <- rep(1:ceiling(length(x)/L), each = 10)
  sigma <- sigma[1:length(x)]
  sigma <- sqrt(mean(unlist(lapply(split(x, sigma), var))))
  diff <- t*(sqrt( (2*sigma)/L ))
  
  ## step 3
  mu1 <- mean(X[1:l])
  mu1up <- mu1 + diff
  mu1lo <- mu1 - diff
  
  mus <- vector()
  mus <- append(mus, mu1)
  regimes <- vector()
  
  ## step 4 to 7
  from <- 1
  while(from < length(X)){
    ## search regime shift
    R_new <-search.regime.start(
        timeseries = X,
        from = from,
        L = L,
        up = mu1up,
        lo = mu1lo,
        diff = diff
      )
    regimes <- append(regimes, R_new)
    
    ## calculate RSI
    rsi <- RSI(timeseries = X, regime.start = R_new, L = L, diff = diff)
    
    if(!any(rsi<0)){
      
      i_end <- ifelse((R_new+L) > length(timeseries), length(timeseries), (R_new+L))
      mus <- append(mus, mean(X[R_new:i_end]))
      ## update
      R_new_update <- mean(X[R_new:(R_new+L)])
      
      mu1up <- R_new_update + diff
      mu1lo <- R_new_update - diff
      
      from <- R_new + 1
    }else{
      from <- from + 1
    }
    if((from  + L) >= length(timeseries)) 
      break
  }
  regimes <- list(regimes, mus)
  
  regimes <- rbind(c(0, regimes[[1]], length(X)), c(1, regimes[[2]]))
  
  regimes <- foreach(j = 1:(ncol(regimes)-1), .combine = "c") %do% {
    st <- (regimes[1,j]+1):regimes[1,j+1]
    rep(regimes[2,j+1], length(st))
  }
  
  return(regimes)
}

## Example data
# # rand data
x <- c(0.1, 0.9, 0.91, 0.92, 0.7, 0.8, -0.1, 1.3,
       0.2, -0.2, -1, -1.8, -0.1, 0.3, -0.3, -0.4,
       -0.5, -1, -1.2, -1.2, -1.3, -0.9, 1, 0.8, 1.2,
       -0.1, 0.4, 1, 1.2)
plot(x, type ="b")
# 
rs <- rodionov(timeseries = x, L = 10, p = 0.05)
lines(rs)
