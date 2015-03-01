#' simulate a sample of one network from stationary distribution
#' 
#' This function simulates networks from the stationary distribution
#' starting at network \code{g}, and saves the network sampled statistics
#' in data.frame \code{netsample}: each column is a network statistics.
#' @param ppp vector of indicators for the utility functions
#' @param g the matrix containing the network
#' @param x the matrix of exogenous variables
#' @param dt indicators of which statistics to compute
#' @param samplesize length of the sample
#' @param pr probability of changing a row
#' @param pc probability of changing a column
#' @param mtemp number of temperatures for simulated tempering (not implemented)
#' @param seed sets the seed for the random number generator (an integer)
#' 
#' @author Angelo Mele 
#' @examples 
#' dt  <- matrix(data = c(1,1,1,1,1,1), nrow = 3, ncol = 2)
#' p   <- c(1, 2, 3)
#' x   <- matrix(data = rnorm(5,5,1), nrow = 5, ncol = 1)
#' g   <- matrix(data = rbinom(25,1,.5), nrow = 5, ncol = 5)
#' a   <- c(-4,-1,1)
#' ans <- samplemany(g, x, a, p, dt, skip=100, samplesize = 1000)
#' ans

sampleone <- function(g, x, a, ppp, dt, samplesize, 
                       pr = 0, pc = 0, mtemp = 2, seed = 1977 ) {
        ts <- matrix(0, nrow=1, ncol=ppp[3])
        netsample <- .Fortran("sample_ast_one", 
                nn = as.integer(dim(g)[1]),
                qq = as.integer(dim(x)[2]),
                pp = as.integer(ppp),
                g = matrix(as.integer(g), ncol=ncol(g)), 
                x = matrix(as.double(x), ncol=ncol(x)),  
                aa = as.double(a),
                dtt = matrix(as.integer(dt),ncol = ncol(dt)), 
                sample = as.integer(samplesize),
                mtemp = as.integer(mtemp),
                qup = as.double(.5),
                #tout = matrix(ts,ncol=samplesize),
                tout = matrix(as.double(ts), nrow=1),
                #kout = as.integer(1),
                seed1 = as.integer(seed),
                pr = as.double(pr),
                pc = as.double(pc)
                )
  netsample <- data.frame(matrix(netsample$tout,ncol=ppp[3],nrow=1))
  rm(ts)
  return(netsample)
  #return(ans)
}
