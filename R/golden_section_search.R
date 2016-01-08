#' Golden section search
#' 
#' Uses golden section method to search given space for value that minimizes
#' function given
#' 
#' @keywords internal
#' @param f function
#' @param lower Lower bound to search in
#' @param upper Upper bound to search in
#' @param tol How exact the answer must be
#' @return Returns a value that minimizes the function input
#' @export
golden.section.search = function(f, lower, upper, tol = 1e-6){
  
  golden.ratio <- (sqrt(5)-1)/2
  
  c <- upper-(golden.ratio*(upper-lower))
  d <- lower+(golden.ratio*(upper-lower))
  
  while(abs(upper-lower)>tol){
    
    fc <- f(c)
    fd <- f(d)
    
    if(fc<fd){
      
      upper <- d
      d <- c
      c <- upper-golden.ratio*(upper-lower)
      
    }
    
    else{
      
      lower <- c
      c <- d
      d <- lower + golden.ratio*(upper-lower)
      
    }
    
  }
  
  return((lower+upper)/2)
  
}