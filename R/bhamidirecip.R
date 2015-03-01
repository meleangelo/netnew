#' computes solutions to equation in Bhamidi et al 2011 to check for high/low temp phase (reciprocity model)
#' 
#' this function computes solutions to equation in Bhamidi et al 2011 for the directed network.
#' Returns a matrix with solution (row 1) and first derivative at the solution (row 2)
#' 
#' @param a parameter for edges term
#' @param b parameter for reciprocity term
#' 
#' @author Angelo Mele 
#' @examples 
#' bhamidirecip(a = -2, b = 0.2)
#' 
#' 
bhamidirecip <- function(a, b) {    
    solution<-uniroot.all(phi_recip, c(0, 1), tol = 0.0001, a = a, b = b)
    derivative<-phi_prime_recip(solution, a, b)
    bhamidirecip <- rbind(solution,derivative)
    return(bhamidirecip)
}


