#' Values for 100 randomly generated tumor samples
#' 
#' A dataset containing the value for the 96-trinucleotide contexts of
#' 100 randomly generated tumors.
#' 
#' @docType data
#' @keywords datasets
#' @name randomly.generated.tumors
#' @format A data frame of 100 rows and 96 columns
NULL

#' Published Signatures from Alexandrov et al 2013
#' 
#' A dataset containing the published signatures from Alexandrov
#' et al read into R as a data frame. Data found: \url{ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl}
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name signatures
#' @format A data frame of 27 rows and 96 columns
NULL

#' Published Signatures from Sanger COSMIC
#' 
#' A dataset containing the additional signatures identified
#' read into R as a data frame. Data found: \url{http://cancer.sanger.ac.uk/cosmic/signatures}
#' 
#' 
#' @docType data
#' @keywords datasets
#' @name updated.signatures
#' @format A data frame of 30 rows and 96 columns
NULL

#' The counts of every trinuclotide frequency in an exome
#' 
#' A datset containing the number of times each trinucleotide (ex: ACA)
#' is found in the exome region captured by sequencing.
#' 
#' @docType data
#' @keywords datasets
#' @name tri.counts.exome
#' @format A data frame of 32 rows and 1 column that contains the counts
NULL

#' Example output of whichSignatures()
#' 
#' A list that was generated from running whichSignatures().  It contains
#' the following:
#' 
#' \itemize{
#'   \item weights - weight of each signature found in the sample
#'   \item tumor - tricontext fractions that were used as input
#'   \item product - product of the weights by the reference signatures used
#'   \item diff - difference between the tumor tricontexts fractions used as input and those generated as the product
#'   \item unknown - the fraction of the tumor profile that could not be assigned to a signature
#' }
#' 
#' @docType data
#' @keywords datasets
#' @name example.output
#' @format A list of four elements
NULL

#' Example input to mut.to.sigs.input()
#' 
#' A data frame containing example mutations that can be used as input to mut.to.sigs.input().  Contains the following columns:
#' 
#' \itemize{
#'   \item Sample - sample name
#'   \item chr - chromosome number
#'   \item pos - chromosome position
#'   \item ref - reference base
#'   \item alt - alternate base
#'}
#' 
#' @docType data
#' @keywords datasets
#' @name sample.mut.ref
#' @format A data frame of 11 rows and 5 column that contains example mutations which could be used in mut.to.sigs.input()
NULL



