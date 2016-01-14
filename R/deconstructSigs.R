#' deconstructSigs
#'
#' Takes sample information in the form of the fraction of mutations
#' in each of 96 trinucleotide contexts and identifies the weighted combination 
#' of published signatures that, when summed, most closely reconstructs the 
#' mutational profile.
#'
#' @section Main functions:
#' \itemize{ 
#'   \item \code{\link{whichSignatures}} 
#'   \item \code{\link{mut.to.sigs.input}} 
#'   \item \code{\link{plotSignatures}} 
#'   \item \code{\link{makePie}}
#'   }
#'   
#' @docType package
#' @author Rachel Rosenthal rachel.rosenthal.14@ucl.ac.uk
#' @name deconstructSigs
NULL