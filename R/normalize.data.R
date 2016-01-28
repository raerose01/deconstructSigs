#' Normalizes trinucleotide contexts
#' 
#' Normalizes the trinucleotide contexts
#' 
#' @keywords internal
#' @param col column names
#' @param trimer.counts count of the number of times each trimer is found in the area sequenced
#' @return Returns a normalized column based on the trimer counts
#' @export
norm.it <- function(col, trimer.ratio){
  trimer  <- paste(substr(colnames(col), 1, 1), substr(colnames(col), 3, 3), substr(colnames(col), 7, 7), sep = "")
  new.col <- col*trimer.ratio[trimer,]
  return(new.col)
}

#' Calculates trinucleotide context fraction
#' 
#' Determines trinucleotide context fraction from a mutation counts list. Given
#' a data frame or .txt file or mutation counts with columns as trincucleotide
#' contexts in the form A[C>A]A and the rows as sample id's and a data frame or
#' .txt file of the times each tri nucleotide context is seen in the sequencing
#' region with the rows in the format ACA, the function returns a data frame of
#' the mutation counts normalized by trinucleotide frequency
#' 
#' @keywords internal
#' @param mut.counts.ref data frame of counts of mutations in each trinucleotide
#'   context for each sample or location where .txt file is found
#' @param trimer.counts.ref data frame of counts of times each trinculeotide
#'   context is seen in sequencing area or location where the .txt file is found
#' @return Returns the trinucleotide context fraction
#' @export
getTriContextFraction <- function(mut.counts.ref, trimer.counts.method){
    
  if(exists("mut.counts.ref", mode = "list")){
    mut.counts <- mut.counts.ref
  } else {
    if(file.exists(mut.counts.ref)){
      mut.counts <- utils::read.table(mut.counts.ref, sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE)
    } else {
      print("mut.counts.ref is neither a file nor a loaded data frame")
    }
  }
  
  if(class(trimer.counts.method) == 'character'){
  
    # if default, return mut counts that sum to 1
    if(trimer.counts.method == 'default'){
      # make each row sum to 1
      norm.mut.counts           <- mut.counts/rowSums(mut.counts)
      return(norm.mut.counts)
    }
    
    # return mut counts divided by number of times that trinucleotide context is observed in the genome
    if(trimer.counts.method == 'genome'){
  
      tri.counts.wgs <- tri.counts.genome
      multiplicative.ratio      <- 1/tri.counts.wgs
      norm.mut.counts           <- sapply(colnames(mut.counts), function(x) {norm.it(mut.counts[,x,drop=F], trimer.ratio = multiplicative.ratio)})
      norm.mut.counts           <- data.frame(norm.mut.counts, row.names = rownames(mut.counts))
      colnames(norm.mut.counts) <- colnames(mut.counts)
      
      # make each row sum to 1
      norm.mut.counts           <- norm.mut.counts/rowSums(norm.mut.counts)
      return(norm.mut.counts)
   
    }
    
    # return mut counts divided by number of times that trinucleotide context is observed in the exome
    if(trimer.counts.method == 'exome'){
      
      tri.counts.wes <- tri.counts.exome
      multiplicative.ratio      <- 1/tri.counts.wes
      norm.mut.counts           <- sapply(colnames(mut.counts), function(x) {norm.it(mut.counts[,x,drop=F], trimer.ratio = multiplicative.ratio)})
      norm.mut.counts           <- data.frame(norm.mut.counts, row.names = rownames(mut.counts))
      colnames(norm.mut.counts) <- colnames(mut.counts)
      
      # make each row sum to 1
      norm.mut.counts           <- norm.mut.counts/rowSums(norm.mut.counts)
      return(norm.mut.counts)
    
    }
    
    # use the ratio of WGS/WES to normalize the input data
    if(trimer.counts.method == 'exome2genome'){
      
      tri.counts.wgs <- tri.counts.genome
      tri.counts.wes <- tri.counts.exome
  
      # multiply by WGS/WES tricontext ratio
      wgs.wes.ratio             <- tri.counts.wgs/tri.counts.wes
      norm.mut.counts           <- sapply(colnames(mut.counts), function(x) {norm.it(mut.counts[,x,drop=F], trimer.ratio = wgs.wes.ratio)})
      norm.mut.counts           <- data.frame(norm.mut.counts, row.names = rownames(mut.counts))
      colnames(norm.mut.counts) <- colnames(mut.counts)
      
      # make each row sum to 1
      norm.mut.counts           <- norm.mut.counts/rowSums(norm.mut.counts)
      return(norm.mut.counts)
    
    }
    
    if(!trimer.counts.method %in% c('default', 'exome', 'genome', 'exome2genome')){
      
      stop(paste(trimer.counts.method, ' is not set to one of the available options (\'default\', \'exome\', \'genome\', \'exome2genome\') and is not a data frame.', sep = ''))
      
    }
    
  }
  
  if(class(trimer.counts.method) %in% c('data.frame', 'matrix')){
    
    tri.counts.ratio <- trimer.counts.method

    norm.mut.counts           <- sapply(colnames(mut.counts), function(x) {norm.it(mut.counts[,x,drop=F], trimer.ratio = tri.counts.ratio)})
    norm.mut.counts           <- data.frame(norm.mut.counts, row.names = rownames(mut.counts))
    colnames(norm.mut.counts) <- colnames(mut.counts)
    
    # make each row sum to 1
    norm.mut.counts           <- norm.mut.counts/rowSums(norm.mut.counts)
    return(norm.mut.counts)
    
  }
  
}

# From pipeline version to match signatures.txt file
#' Re-formats column names
#' 
#' Changes the way the column names are formatted from nxn.y to the n[x>y]n format
#' 
#' @keywords internal
#' @param col Vector of column names
#' @return Returns a newly formatted vector of column names
#' @export
changeColumnNames <- function(col){
  new.col <- paste(substr(col, 1, 1), "[", substr(col, 2, 2), ">", substr(col, 8, 8), "]", substr(col, 3, 3), sep = "")
  return(new.col)
}

################################################