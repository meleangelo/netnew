#' simulate a sample from directed erdos-renyi model
#' 
#' This function simulates a network from the erdos-renyi model with 
#' probability of link \code{p} and \code{n} nodes. The result is 
#' an adjacency matrix
#' @param p probability of a link
#' @param n number of nodes of the network
#' 
#' @author Angelo Mele 
#' @examples 
#' g <- erdos(p = 0.5, n = 10)
#' g
#' 
#' 

erdos <- function(p,n) {
    g <- matrix( data = rbinom(n = n^2, size = 1, prob = p) , ncol = n, nrow = n)
    diag(g) <- 0
    return(g)
}