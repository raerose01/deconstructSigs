
#' Calculates error
#' 
#' Finds error between calculated trinucleotide context fractions and inputted
#' ones
#' 
#' @keywords internal
#' @param tumor Actual trinucleotide context fractions
#' @param signatures Signatures matrix
#' @param w Weights matrix
#' @return Returns the sum squared error between calculated and actual
#'   trinucleotide context fractions
#' @export
getError = function(tumor, signatures, w){
  w_norm = w/sum(w)
  product = w_norm %*% signatures
  error = tumor - product
  penalized.error = error
  neg = which(penalized.error < 0)
  penalized.error[,neg] = 1*penalized.error[,neg]
  tot = sum(penalized.error * penalized.error)
  return(tot)
} 

#' Seeds weight matrix
#' 
#' Determines which single signature results in the lowest sum-squared error and
#' uses that as a seed to start from
#' 
#' @keywords internal
#' @param tumor Actual trinucleotide context fractions
#' @param signatures Signatures matrix
#' @return Returns the index corresponding to the signature that best describes
#'   the input data
#' @export
findSeed = function(tumor, signatures){
  w0 = vector()
  for(i in 1:nrow(signatures)){
    seed_matrix = matrix(0, nrow = nrow(tumor), ncol = nrow(signatures))
    seed_matrix[i] = 1
    tot = getError(tumor, signatures, seed_matrix)
    w0 = rbind(w0, tot)
  }
  return(which(w0 == min(w0))[1])
}

#' Updates the weights matrix
#' 
#' Determines what proportion of what signature, when added to the current
#' weights matrix, results in the lowest sum-squared error
#' 
#' @keywords internal
#' @param tumor Actual trinucleotide context fractions
#' @param signatures Signatures matrix
#' @param w Weights matrix
#' @param signatures.limit Number of signatures to limit the search to
#' @return Returns an updated weights matrix
#' @export
updateW_GR = function(tumor, signatures, w, signatures.limit, bound = 100){
  error_old = getError(tumor, signatures, w)
  boo = matrix(+Inf, nrow = 1, ncol = nrow(signatures))
  v = matrix(0, nrow = nrow(signatures), ncol = nrow(signatures))
  
  # Add another signature if needed
  if(length(w[w!="0"]) < signatures.limit){
    for(i in 1:nrow(signatures)){
      tmp = matrix(0, nrow = 1, ncol = nrow(signatures))
      toMinimize = function(x){
        tmp[1,i] = x
        w_new = w + tmp[1,]
        getError(tumor, signatures, w_new)
      }
      v[i,i] = golden.section.search(toMinimize, -w[i], bound, tol = 1e-6)
      w_new = w + v[i,]
      boo[1,i] = getError(tumor, signatures, w_new)
    }
  }
  
  # Only look at the signatures already present if there are too many being added
  if(length(w[w!="0"]) >= signatures.limit){
    signatures.available.to.change = colnames(w)[w!="0"]
    indices.for.signatures.available = which(rownames(signatures) %in% signatures.available.to.change)
    for(i in indices.for.signatures.available){
      #print(i)
      tmp = matrix(0, nrow = 1, ncol = nrow(signatures))
      toMinimize = function(x){
        tmp[1,i] = x
        w_new = w + tmp[1,]
        getError(tumor, signatures, w_new)
      }
      v[i,i] = golden.section.search(toMinimize, -w[i], bound, tol = 1e-6)
      w_new = w + v[i,]
      boo[1,i] = getError(tumor, signatures, w_new)
    }
  }
  
  ind = which(boo == min(boo))[1]
  #print(paste("old error: ", error_old, sep = ''))
  #print(paste("new error: ", min(boo), sep = ''))
  
  if(min(boo) < error_old){
    w[ind] = w[ind] + v[ind,ind]
  }
  
  return(w)
}

################################################