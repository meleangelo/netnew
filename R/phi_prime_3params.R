#' derivative of function from Bhamidi et al 2011 to check for high/low temp phase (3 parameters model)
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
#' phi_prime_3params(p = 0.2, a = -2, b = 0.2, c = 0.1)
#' 
#' 

phi_prime_3params<-function(p,a,b,c) {
    phi_prime_3params<-(b+2*c)*exp(a+b*p+2*c*p)/(( 1+ exp(a+b*p+2*c*p) )^2)
    return(phi_prime_3params)
}