#' simulate selected models with few parameters
#' 
#' @description runs a simulation of simple models (erdos-renyi, mutual links, etc)
#' default is erdos-renyi model
#' 
#' @param model a string indicating which model you want to 
#' simulate. The default is \code{"erdos-renyi"}. Other choices implemented are
#' a model with mutual links (\code{"recip"}) and one with 
#' 3 parameters (\code{"3params"})
#' @param initg the initial network for the simulations
#' @param a the vector of parameters
#' @param samplesize the number of simulated networks to collect (default is
#' \code{samplesize =  100})
#' @param skip how many networks to skip between sampled networks, with
#' default \code{skip = 10}.
#' @param pr probability of swapping a row of network matrix 
#' (default \code{pr = 0})
#' @param pc probability of swapping a column of network matrix 
#' (default \code{pc = 0})
#' @param pf probability of updating a large number of links
#' (default \code{pf = 0})
#' @param pb probability of updating using Bhamidi-like step 
#' (default \code{ob = 0})
#' @param nmodes number of modes of gibbs distribution
#' (default \code{nmodes = 1})
#' @param modes modes of the gibss distribution 
#' (default \code{modes = 0.5})
#' @param sizelargestep number of links to update in large step 
#' (default \code{sizelargestep = 10})
#' @param mtemp number of temperatures for simulated tempering (not implemented)
#' @param seed seed for random number generator
#' 
#' @seealso \link[netnew]{samplemanybham, samplemodellarge}
#' 
#' @examples 
#' n  <- 100
#' dt <- matrix(data = c(1,1,1,1,1,1), nrow = 3, ncol = 2) 
#' p  <- c(1, 2, 3) 
#' x  <- matrix(data = rnorm(n,5,1), nrow = n, ncol = 1) 
#' g  <- matrix(data = rbinom(n*n,1,.5), nrow = n, ncol = n) 
#' a  <- c(-2,.5,.01) 
#' # show data 
#' p; g; x; dt
#' tt <- simulatemodelbham("3params", g, a, skip = 1000, samplesize = 10000)
#' tt


simulatemodelbham <- function(model = "erdos-renyi", 
                          initg, a, samplesize = 100, skip = 10, 
                          pr = 0, pc = 0, pf = 0, pb = 0, pinv = 0,
                          nmodes = 1, modes = 0.5, sizelargestep = 10,
                          mtemp = 2, seed = 1977) {
    if (model == "erdos-renyi") {
        p  <- c(1,1,1)
        dt <- matrix(data = c(1,1), nrow = 1, ncol = 2)
        x  <- matrix(data = 0, nrow = dim(initg)[1], ncol = 1)
        tsamp <- samplemanybham(g = initg, x = x, a = a, ppp = p, dt = dt, 
                            skip = skip , samplesize = samplesize,
                            pr = pr, pc = pc, pf = pf, pb = pb, pinv = pinv,
                            nmodes = nmodes, modes = modes, 
                            sizelargestep = sizelargestep,
                            mtemp = mtemp, seed = seed)
    
        colnames(tsamp) <- "links"
    } else if (model == "recip") {
        p  <- c(1,2,2)
        dt <- matrix(data = c(1,1,1,1), nrow = 2, ncol = 2)
        x  <- matrix(data = 0, nrow = dim(initg)[1], ncol = 1)
        tsamp <- samplemanybham(g = initg, x = x, a = a, ppp = p, dt = dt, 
                            skip = skip , samplesize = samplesize,
                            pr = pr, pc = pc, pf = pf, pb = pb, pinv = pinv,
                            nmodes = nmodes, modes = modes, 
                            sizelargestep = sizelargestep,
                            mtemp = mtemp, seed = seed)
        
        colnames(tsamp) <- c("links","mutual")
    } else if (model == "3params") {
        p  <- c(1,2,3)
        dt <- matrix(data = c(1,1,1,1,1,1), nrow = 3, ncol = 2)
        x  <- matrix(data = 0, nrow = dim(initg)[1], ncol = 1)
        allmus <- bhamidi3params(a[1],a[2],a[3]*dim(initg)[1])
        modes <- allmus[1,allmus[2,]<1]
        nmodes <- length(modes)
        tsamp <- samplemanybham(g = initg, x = x, a = a, ppp = p, dt = dt, 
                            skip = skip , samplesize = samplesize,
                            pr = pr, pc = pc, pf = pf, pb = pb, pinv = pinv,
                            nmodes = nmodes, modes = modes, 
                            sizelargestep = sizelargestep,
                            mtemp = mtemp, seed = seed)
        
        colnames(tsamp) <- c("links","mutual","indirect")
        
    } else if (model == "edgetriangle") {
        p  <- c(1,1,2)
        dt <- matrix(data = c(1,1,1,100), nrow = 2, ncol = 2)
        x  <- matrix(data = 0, nrow = dim(initg)[1], ncol = 1)
        tsamp <- samplemanybham(g = initg, x = x, a = a, ppp = p, dt = dt, 
                                skip = skip , samplesize = samplesize,
                                pr = pr, pc = pc, pf = pf, pb = pb, pinv = pinv,
                                nmodes = nmodes, modes = modes, 
                                sizelargestep = sizelargestep,
                                mtemp = mtemp, seed = seed)
        
        colnames(tsamp) <- c("links","triangles")
        
    } else if (model == "eit") {
        p  <- c(1,2,4)
        dt <- matrix(data = c(1,1,1,1,1,1,1,100), nrow = 4, ncol = 2)
        x  <- matrix(data = 0, nrow = dim(initg)[1], ncol = 1)
        tsamp <- samplemanybham(g = initg, x = x, a = a, ppp = p, dt = dt, 
                                skip = skip , samplesize = samplesize,
                                pr = pr, pc = pc, pf = pf, pb = pb, pinv = pinv,
                                nmodes = nmodes, modes = modes, 
                                sizelargestep = sizelargestep,
                                mtemp = mtemp, seed = seed)
        
        colnames(tsamp) <- c("links","mutual","indirect","triangles")
    }
    else {
        cat("\n sorry, this model is not implemented yet
            \n \n try 'erdos-renyi' or 'recip'")
        tsamp <- c()
    }
    return(tsamp)
}