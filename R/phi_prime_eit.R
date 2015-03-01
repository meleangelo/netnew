#' derivative of function from Bhamidi et al 2011 to check for high/low temp phase  (edges + indirect + triangles model)
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
#' phi_prime_eit(p = 0.2, a = -2, b = 0.2, c = 0.1)
#' 
#' 

phi_prime_eit<-function(p,a,b,c) {
    phi_prime_eit<-(2*b+2*c*p)*exp(a+2*b*p+c*p^2)/(( 1+ exp(a+2*b*p+c*p^2) )^2)
    return(phi_prime_eit)
}