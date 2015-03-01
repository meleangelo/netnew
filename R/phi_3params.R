#' function from Bhamidi et al 2011 to check for high/low temp phase (3 pararameters model)
#' 
#' this function is the same as Bhamidi et al 2011 for the directed network 
#' 
#' @param p prob of link (mean field approximation)
#' @param a parameter for edges term
#' @param b parameter for reciprocity term
#' @param c parameter for indirect links
#' 
#' @author Angelo Mele 
#' @examples 
#' phi_3params(p = 0.2, a = -2, b = 0.2, c = 0.1)
#' 
#' 

phi_3params <- function(p,a,b,c) {
    phi_3params <- exp(a+b*p+2*c*p)/( 1+ exp(a+b*p+2*c*p) ) - p
    return(phi_3params)
}
