#' estimate a network formation model using the approximate exchange algorithm
#' 
#' 
#' This function estimates the posterior of structural parameters of 
#' a network formation model, using the approximate exchange algorithm
#' used in Mele (2015). For each proposed parameter, it simulates 
#' networks from the stationary distribution
#' starting at network \code{g}, and saves the sample in data.frame
#' \code{parsample}: each column is a parameter associated with a 
#' network statistics. The 
#' algorithm can perform small steps (1 links update per iteration),
#' intermediate steps (1 row or 1 column update per iteration),
#' or large step (update a number of links determined by the control
#' variable \code{sizelargestep}). If the model is an homogenous
#' player model, we can also use the Bhamidi et al (2011) approach
#' to establish the modes of the gibbs distribution and use that
#' in the simulation, i.e. we propose to jump to one of the 
#' modes with probability \code{pb}. 
#' @param ppp vector of indicators for the utility functions
#' @param g the matrix containing the network
#' @param x the matrix of exogenous variables
#' @param a0 parameter vector at which the simulation is started
#' @param dt indicators of which statistics to compute
#' @param parsample number of parameters to sample
#' @param netsample number of networks to sample
#' @param skip how many networks to skip between consecutive sample
#' @param pr probability of changing a row
#' @param pc probability of changing a column
#' @param pf probability of large step 
#' @param pb probability of updating using Bhamidi-like step 
#' (default \code{pb = 0})
#' @param pinv probability of inversion of adjacency matrix
#' @param nmodes number of modes of gibbs distribution
#' (default \code{nmodes = 1})
#' @param modes modes of the gibbs distribution 
#' (default \code{modes = 0.5})
#' @param sizelargestep how many links to update in a large step
#' @param mtemp number of temperatures for simulated tempering (not implemented)
#' @param seed sets the seed for the random number generator (an integer)
#' @param mu0 vector of means for the prior
#' @param sigma0 vector of variances for the prior
#' @param covpars matrix of covariances among parameters (for sampler)
#' 
#' @author Angelo Mele 
#' @examples 
#' dt  <- matrix(data = c(1,1,1,1,1,1), nrow = 3, ncol = 2)
#' p   <- c(1, 2, 3)
#' x   <- matrix(data = rnorm(5,5,1), nrow = 5, ncol = 1)
#' g   <- matrix(data = rbinom(25,1,.5), nrow = 5, ncol = 5)
#' a0   <- c(-4,-1,1)
#' mu0  <- c(0,0,0)
#' sigma0 <- c(1,1,1)
#' covpars <- matrix(c(1,0,0,0,1,0,0,0,1), ncol = 3, nrow = 3) 
#' ans <- bayesest_ex(g, x, a0, p, dt, parsample = 100, netsample = 100
#'                      mu0 = mu0 , sigma0 = sigma0, covpars = covpars)
#' ans

bayesest_ex <- function(g, x, a0, ppp, dt, paramsample, netsample, covpars,
                           mu0, sigma0,
                           skip = 1, pr = 0, pc = 0, pf = 0, pinv = 0,
                           pb = 0, nmodes = 1, modes = 0.5, sizelargestep = 10, 
                           mtemp = 2, seed = 1977 ) {
    
    #step <- matrix(round(t(chol(covpars)), 10), ncol = ncol(covpars), nrow = nrow(covpars))
    #ts <- matrix(0, nrow=paramsample, ncol=ppp[3])
    parsample <- .Fortran("bayesest_ex", 
                          nn = as.integer(dim(g)[1]),
                          qq = as.integer(dim(x)[2]),
                          pp = as.integer(ppp),
                          g = matrix(as.integer(g), ncol=ncol(g)), 
                          x = matrix(as.double(x), ncol=ncol(x)),  
                          aa0 = as.double(a0),
                          dtt = matrix(as.integer(dt),ncol = ncol(dt)), 
                          skip = as.integer(skip),
                          parsamp = as.integer(paramsample),
                          samp = as.integer(netsample),
                          #step = matrix(as.integer(step), ncol = ncol(step)),
                          #step = matrix(as.double(step), ncol = ncol(step)),
                          step = matrix(as.double(covpars), ncol = ncol(covpars)),
                          mtempnet = as.integer(mtemp),
                          qup = as.double(.5),
                          #tout = matrix(ts,ncol=samplesize),
                          #tout = matrix(as.double(ts), ncol=samplesize),
                          #kout = as.integer(1),
                          seed1 = as.integer(seed),
                          pr = as.double(pr),
                          pc = as.double(pc),
                          pf = as.double(pf),
                          pb = as.double(pb),
                          pinv = as.double(pinv),
                          nmus = as.integer(nmodes),
                          mus = as.double(modes),
                          size = as.integer(sizelargestep),
                          mu0 = as.double(mu0),
                          sigma0 = as.double(sigma0),
                          parout = matrix(0,nrow=paramsample,ncol = ppp[3])
    )
    
 
    
    parsample <- data.frame(matrix(parsample$parout,ncol=ppp[3],nrow=paramsample))
    colnames(parsample)<- paste("theta", seq(1,length(a0),1), sep ="")
    #rm(ts)
    return(parsample)
    #return(ans)
}
