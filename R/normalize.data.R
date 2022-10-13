#' Calculates trinucleotide context fraction
#'
#' Determines trinucleotide context fraction from a mutation counts list. Given
#' a data frame or .txt file or mutation counts with columns as trincucleotide
#' contexts in the form A[C>A]A and the rows as sample id's and a data frame or
#' .txt file of the times each tri nucleotide context is seen in the sequencing
#' region with the rows in the format ACA, the function returns a data frame of
#' the mutation counts normalized by trinucleotide frequency
#'
#' @param mut.counts data frame of counts of mutations in each trinucleotide
#'   context for each sample
#' @return Returns the trinucleotide context fraction
#' @keywords internal
getTriContextFraction <- function(mut.counts, trimer.counts.method,
                                  genome.ref = NULL, chr.list = NULL,
                                  exome.range = NULL) {
  # if default, return mut counts that sum to 1
  if (identical(trimer.counts.method, "default")) {
    # make each row sum to 1
    norm.mut.counts <- mut.counts / rowSums(mut.counts)
    return(norm.mut.counts)
  }

  if (grepl("genome", trimer.counts.method)) {
    if (is.null(genome.ref)) {
      genome.ref <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
      tri.counts.wgs <- tri.counts.genome
    } else if (identical(genome.ref, BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)) {
      tri.counts.wgs <- tri.counts.genome
    } else {
      tri.counts.wgs <- Biostrings::trinucleotideFrequency(
        BSgenome::getSeq(
          genome.ref,
          intersect(GenomeInfoDb::seqnames(genome.ref), chr.list)
        )
      )
      tri.counts.wgs <- colSums(tri.counts.wgs)[make_n_nucleotide(3L)]
    }
  }

  if (grepl("exome", trimer.counts.method)) {
    if (is.null(exome.range)) {
      tri.counts.wes <- tri.counts.exome
    } else {
      exome.range <- GenomicRanges::reduce(exome.range)
      exome.range <- GenomeInfoDb::keepSeqlevels(
        exome.range, chr.list,
        pruning.mode = "tidy"
      )
      tri.counts.wes <- Biostrings::trinucleotideFrequency(
        BSgenome::getSeq(genome.ref, exome.range)
      )
      tri.counts.wes <- colSums(tri.counts.wes)[make_n_nucleotide(3L)]
    }
  }

  trimer.ratio <- switch(trimer.counts.method,
    # return mut counts divided by number of times that trinucleotide context is observed in the genome
    genome = 1L / tri.counts.wgs,

    # return mut counts divided by number of times that trinucleotide context is
    # observed in the exome
    exome = 1L / tri.counts.wes,

    # use the ratio of WGS/WES to normalize the input data
    exome2genome = tri.counts.wgs / tri.counts.wes,

    # use the ratio of WES/WGS to normalize the input data
    genome2exome = tri.counts.wes / tri.counts.wgs
  )

  norm.mut.counts <- sapply(colnames(mut.counts), function(x) {
    trimer <- strsplit(x, split = "")[[1L]]
    trimer <- paste0(trimer[c(1L, 3L, 7L)], collapse = "")
    mut.counts[[x]] * trimer.ratio[[trimer]]
  })
  norm.mut.counts <- data.frame(
    norm.mut.counts,
    row.names = rownames(mut.counts)
  )
  colnames(norm.mut.counts) <- colnames(mut.counts)

  # make each row sum to 1
  norm.mut.counts / rowSums(norm.mut.counts)
}

# From pipeline version to match signatures.txt file
#' Re-formats column names
#'
#' Changes the way the column names are formatted from nxn.y to the n[x>y]n format
#'
#' @param col Vector of column names
#' @return Returns a newly formatted vector of column names
#' @keywords internal
changeColumnNames <- function(col) {
  new.col <- paste(substr(col, 1, 1), "[", substr(col, 2, 2), ">", substr(col, 8, 8), "]", substr(col, 3, 3), sep = "")
  return(new.col)
}

################################################
