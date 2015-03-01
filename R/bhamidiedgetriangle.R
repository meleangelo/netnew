#' computes solutions to equation in Bhamidi et al 2011 to check for high/low temp phase (edge - triangle model)
#' 
#' this function computes solutions to equation in Bhamidi et al 2011 for the directed network.
#' Returns a matrix with solution (row 1) and first derivative at the solution (row 2)
#' 
#' @param a parameter for edges term
#' @param b parameter for cyclic triangles term
#' 
#' @author Angelo Mele 
#' @examples 
#' bhamidiedgetriangle(a = -2, b = 0.2)
#' 
#' 
bhamidiedgetriangle <- function(a, b) {    
    solution<-uniroot.all(phi_edgetriangle, c(0, 1), tol = 0.0001, a = a, b = b)
    derivative<-phi_prime_edgetriangle(solution, a, b)
    bhamidiedgetriangle <- rbind(solution,derivative)
    return(bhamidiedgetriangle)
}


