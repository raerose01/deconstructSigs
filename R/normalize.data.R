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
getTriContextFraction <- function(mut.counts.ref, trimer.counts.genome, trimer.counts.exome){
    
  if(exists("mut.counts.ref", mode = "list")){
    mut.counts <- mut.counts.ref
  } else {
    if(file.exists(mut.counts.ref)){
      mut.counts <- utils::read.table(mut.counts.ref, sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE)
    } else {
      print("mut.counts.ref is neither a file nor a loaded data frame")
    }
  }
  
  if(exists("trimer.counts.genome", mode = "list")){
    tri.counts.wgs <- trimer.counts.genome
  } else {
    if(file.exists(trimer.counts.genome)){
      tri.counts.wgs <- utils::read.table(trimer.counts.genome, sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE)
    } else {
      print("tri.counts.genome is neither a file nor a loaded data frame")
    }
  }
  
  if(exists("trimer.counts.exome", mode = "list")){
    tri.counts.wes <- trimer.counts.exome
  } else {
    if(file.exists(trimer.counts.exome)){
      tri.counts.wes <- utils::read.table(trimer.counts.exome, sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE)
    } else {
      print("tri.counts.exome is neither a file nor a loaded data frame")
    }
  }
  
  # multiply by WGS/WES tricontext ratio
  wgs.wes.ratio             <- tri.counts.wgs/tri.counts.wes
  norm.mut.counts           <- sapply(colnames(mut.counts), function(x) {norm.it(mut.counts[,x,drop=F], trimer.ratio = wgs.wes.ratio)})
  norm.mut.counts           <- data.frame(norm.mut.counts, row.names = rownames(mut.counts))
  colnames(norm.mut.counts) <- colnames(mut.counts)
  
  # make each row sum to 1
  #norm.mut.counts           <- apply(norm.mut.counts, 1, function(x) {x/sum(x)})
  #norm.mut.counts           <- t(norm.mut.counts)
  norm.mut.counts           <- norm.mut.counts/rowSums(norm.mut.counts)

  
  return(norm.mut.counts)
  
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