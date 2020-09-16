#' Which signatures are present
#' 
#' Determines how much of each signature is present in the sample given
#' 
#' @param tumor.ref Either a data frame or location of input text file, where 
#'   rows are samples, columns are trinucleotide contexts
#' @param sample.id Name of sample -- should be rowname of tumor.ref. Optional
#'   if the tumor.ref contains one single sample
#' @param signatures.ref Either a data frame or location of signature text file,
#'   where rows are signatures, columns are trinucleotide contexts
#' @param associated Vector of associated signatures. If given, will narrow the 
#'   signatures tested to only the ones listed.
#' @param signatures.limit Number of signatures to limit the search to
#' @param signature.cutoff Discard any signature contributions with a weight 
#'   less than this amount
#' @param contexts.needed FALSE if tumor.file is a context file, TRUE if it is 
#'   only mutation counts
#' @param tri.counts.method Set to either:
#' \itemize{
#'  \item 'default' -- no further normalization \item 'exome' -- normalized by
#'   number of times each trinucleotide context is observed in the exome \item
#'   'genome' -- normalized by number of times each trinucleotide context is
#'   observed in the genome \item 'exome2genome' -- multiplied by a ratio of that
#'   trinucleotide's occurence in the genome to the trinucleotide's occurence in
#'   the exome \item 'genome2exome' -- multiplied by a ratio of that
#'   trinucleotide's occurence in the exome to the trinucleotide's occurence in
#'   the genome \item data frame containing user defined scaling factor -- count
#'   data for each trinucleotide context is multiplied by the corresponding value
#'   given in the data frame }
#'   @param iterations Number of iterations to perform
#' @return A list of the weights for each signatures, the product when those are
#'   multiplied on the signatures, the difference between the tumor sample and 
#'   product, the tumor sample tricontext distribution given, and the unknown 
#'   weight.
#' @export
#' @section Normalization: If the input data frame only contains the counts of
#'   the mutations observed in each context, then the data frame must be
#'   normalized. In these cases, the value of `contexts.needed` should be TRUE.
#'   \cr The parameter, `tri.counts.method`, determines any additional
#'   normalization performed. Any user provided data frames should match the
#'   format of `tri.counts.exome` and `tri.counts.genome`. \cr The method of
#'   normalization chosen should match how the input signatures were normalized.
#'   For exome data, the 'exome2genome' method is appropriate for the signatures
#'   included in this package. For whole genome data, use the 'default' method
#'   to obtain consistent results.
#' @examples
#' test = bootstrap_whichSignatures(tumor.ref = randomly.generated.tumors,
#'                        sample.id = "2", 
#'                        contexts.needed = FALSE)


## Bootstrapping function for mmsig

bootstrap_whichSignatures <- function(mut.ref,
                                      sample.id,
                                      run.mut.to.sigs = FALSE,
                                      sample.id.col = "Sample",
                                      chr = "chr",
                                      pos = "pos",
                                      ref = "ref",
                                      alt = "alt",
                                      bsg = NULL,
                                      sig.type = "SBS",
                                      dbs_table = dbs_possible,
                                      signatures.ref = signatures.cosmic,
                                      associated = c(),
                                      signatures.limit = NA,
                                      signature.cutoff = 0.06,
                                      contexts.needed = TRUE,
                                      tri.counts.method = "default",
                                      iterations = 1000)
  {
  
  # will need to run mut.to.sigs.input
  # or need to provide counts matrix, not percentages
  if(run.mut.to.sigs){
    tumor.ref <- mut.to.sigs.input(mut.ref,
                                   sample.id = sample.id.col,
                                   chr = chr,
                                   pos = pos,
                                   ref = ref,
                                   alt = alt,
                                   bsg = bsg,
                                   sig.type = sig.type,
                                   dbs_table = dbs_table)
  }

  if(!run.mut.to.sigs){
    tumor.ref <- mut.ref
  }
  
  # run whichSignatures once
  simple_out <-  deconstructSigs:::whichSignatures(tumor.ref = tumor.ref, 
                                 sample.id = sample.id, 
                                 signatures.ref = signatures.ref, 
                                 associated = associated, 
                                 signatures.limit = signatures.limit, 
                                 signature.cutoff = signature.cutoff, 
                                 contexts.needed = contexts.needed, 
                                 tri.counts.method = tri.counts.method)
  signature_winners <- colnames(simple_out$weights)[which(simple_out$weights != 0)]
  
  unknown_val <- simple_out$unknown
  
  # begin bootstrapping
  
  # mutation types
  classes <- colnames(tumor.ref)
  
  sub <- as.integer(tumor.ref[sample.id,])
  total <- sum(sub)
  
  # sample new 96-classes profiles from the multinomial distribution
  bootMat <- data.frame(t(rmultinom(n = iterations, size = total, prob = sub/total)))
  colnames(bootMat) <- classes
  
  # run whichSignaturs on the bootstrapped matrix
  bootstrap_out <- lapply(rownames(bootMat), FUN = function(i){
    tmp <- whichSignatures(tumor.ref = bootMat, 
                           sample.id = as.numeric(i), 
                           signatures.ref = signatures.ref, 
                           associated = associated, 
                           signatures.limit = signatures.limit, 
                           signature.cutoff = signature.cutoff, 
                           contexts.needed = TRUE, 
                           tri.counts.method = "default")
    return(tmp)
  })
  
  # cheat and any signature not in simple_out is shunted into unknown
  bootstrap_out_clean <- bootstrap_out
  # bootstrap_out_clean <- lapply(1:length(bootstrap_out), FUN = function(x){
  #   tmp <- bootstrap_out[[x]]
  #   added_unknown <- sum(tmp$weights[,-which(colnames(tmp$weights) %in% signature_winners)])
  #   tmp$weights[,-which(colnames(tmp$weights) %in% signature_winners)] <- 0
  #   tmp$unknown <- tmp$unknown + added_unknown
  #   return(tmp)
  # })
  
  # summary statistics
  my_summary <- function(x){
    c(mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))
  }
  
  # get summary of weights
  all_weights <- do.call(rbind, lapply(bootstrap_out_clean, FUN = function(x) {return(x$weights)}))
  all_weights_summary <- data.frame(t(sapply(all_weights, my_summary)))
  colnames(all_weights_summary) <- c('mean', 'CI025', 'CI975')
  
  # get summary of mutation probabilities, by mutation class
  all_mut_probs <- do.call(plyr::rbind.fill, lapply(bootstrap_out_clean, FUN = function(x) {return(x$mutation_probability)}))
  all_mut_probs[is.na(all_mut_probs)] <- 0
  boop <- lapply(unique(all_mut_probs$MutationTypes), FUN = function(x){
    tmp <- all_mut_probs[which(all_mut_probs$MutationTypes == x),]
    out <- data.frame(t(sapply(tmp[!names(tmp) %in% c("Sample.Names", "MutationTypes")], my_summary)))
    names(out) <- c('mean', 'CI025', 'CI975')
    out$MutationTypes <- x
    out$signature <- rownames(out)
    return(out)
  })
  all_mut_probs_summary <- do.call(rbind, boop)
  rownames(all_mut_probs_summary) <- NULL
  
  tmp_out <- list(all_weights_summary, all_mut_probs_summary)
  names(tmp_out) <- c("mutSigsSummary", "mutProbsSummary")

  final_out <- c(simple_out, tmp_out)
  return(final_out)
}
