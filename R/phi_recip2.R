#' function from Bhamidi et al 2011 to check for high/low temp phase (reciprocity model)
#' 
#' this function is the same as Bhamidi et al 2011 for the directed network 
#' 
#' @param p prob of link (mean field approximation)
#' @param a parameter for edges term
#' @param b parameter for reciprocity term
#' 
#' @author Angelo Mele 
#' @examples 
#' phi_recip(p = 0.2, a = -2, b = 0.2)
#' 
#' 

phi_recip2<-function(p,a,b,c) {
    phi_recip2<-(exp(a + 2*c*p) + exp(2*a+b + 2*c*p) )/( 1+ 2*exp(a + 2*c*p) + exp(2*a+b + 2*c*p) ) - p
    return(phi_recip2)
}
