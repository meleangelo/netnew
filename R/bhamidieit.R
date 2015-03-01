#' computes solutions to equation in Bhamidi et al 2011 to check for high/low temp phase  (edges + indirect + triangles model)
#' 
#' this function computes solutions to equation in Bhamidi et al 2011 for the directed network.
#' Returns a matrix with solution (row 1) and first derivative at the solution (row 2)
#' 
#' @param a parameter for edges term
#' @param b parameter for reciprocity term
#' 
#' @author Angelo Mele 
#' @examples 
#' bhamidieit(a = -2, b = 0.2, c = 0.1)
#' 
#' 
bhamidieit <- function(a, b, c) {    
    solution<-uniroot.all(phi_eit, c(0, 1), tol = 0.000001, a = a, b = b, c = c)
    derivative<-phi_prime_eit(solution, a, b, c)
    bhamidieit <- rbind(solution,derivative)
    return(bhamidieit)
}


