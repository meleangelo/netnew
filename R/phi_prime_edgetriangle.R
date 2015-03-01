#' derivative of function from Bhamidi et al 2011 to check for high/low temp phase (edge -triangle  model)
#' 
#' this function is the same as Bhamidi et al 2011 for the directed network 
#' 
#' @param p prob of link (mean field approximation)
#' @param a parameter for edges term
#' @param b parameter for cyclic triangle term
#' 
#' @author Angelo Mele 
#' @examples 
#' phi_prime_edgetriangle(p = 0.2, a = -2, b = 0.2)
#' 
#' 

phi_prime_edgetriangle<-function(p,a,b) {
    phi_prime_edgetriangle<-2*b*exp(a+b*p^2)/(( 1+ exp(a+b*p^2) )^2)
    return(phi_prime_edgetriangle)
}