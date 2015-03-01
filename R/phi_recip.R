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

phi_recip<-function(p,a,b) {
    phi_recip<-exp(a+b*p)/( 1+ exp(a+b*p) ) - p
    return(phi_recip)
}
