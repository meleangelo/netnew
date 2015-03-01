#' Computes network statistics
#' 
#' This function computes all the network statistics for network \code{g}
#' @param ppp vector of indicators for the utility functions
#' @param ggg the matrix containing the network
#' @param xxx the matrix of exogenous variables
#' @param dtt indicators of which statistics to compute
#' 
netstats <- function(ppp, ggg, xxx, dtt ) {
  ans <- .Fortran("sstat", 
                ps = as.integer(ppp), 
                ns = as.integer(dim(ggg)[1]), 
                qs = as.integer(dim(xxx)[2]), 
                g = matrix(as.integer(ggg), ncol=ncol(ggg)), 
                x = matrix(as.double(xxx), ncol=ncol(xxx)),  
                dts = matrix(as.integer(dtt),ncol = ncol(dtt)), 
                ts = as.double(rep(0,dim(dtt)[1]))
                )$ts
  ans <- as.double(ans)
}
