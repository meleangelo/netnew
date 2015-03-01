#' function from Bhamidi et al 2011 to check for high/low temp phase (edge - triangle model)
#' 
#' this function is the same as Bhamidi et al 2011 for the directed network 
#' 
#' @param p prob of link (mean field approximation)
#' @param a parameter for edges term
#' @param b parameter for cyclic triangle term
#' 
#' @author Angelo Mele 
#' @examples 
#' phi_edgetriangle(p = 0.2, a = -2, b = 0.2)
#' 
#' 

phi_edgetriangle<-function(p,a,b) {
    phi_edgetriangle<-exp(a+b*p^2)/( 1+ exp(a+b*p^2) ) - p
    return(phi_edgetriangle)
}
