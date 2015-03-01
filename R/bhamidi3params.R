#' computes solutions to equation in Bhamidi et al 2011 to check for high/low temp phase (3 parameters model)
#' 
#' this function computes solutions to equation in Bhamidi et al 2011 for the directed network.
#' Returns a matrix with solution (row 1) and first derivative at the solution (row 2)
#' 
#' @param a parameter for edges term
#' @param b parameter for reciprocity term
#' 
#' @author Angelo Mele 
#' @examples 
#' bhamidi3params(a = -2, b = 0.2, c = 0.1)
#' 
#' 
bhamidi3params <- function(a, b, c) {    
    solution<-uniroot.all(phi_3params, c(0, 1), tol = 0.0001, a = a, b = b, c = c)
    derivative<-phi_prime_3params(solution, a, b, c)
    bhamidi3params <- rbind(solution,derivative)
    return(bhamidi3params)
}


