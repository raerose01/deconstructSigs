# Signatures: "Signature.1A", "Signature.1B", "Signature.2",  "Signature.3",  "Signature.4",  "Signature.5",  "Signature.6",  "Signature.7",  "Signature.8",  "Signature.9",  "Signature.10", "Signature.11", "Signature.12", "Signature.13", "Signature.14", "Signature.15", "Signature.16", "Signature.17", "Signature.18", "Signature.19", "Signature.20", "Signature.21", "Signature.R1", "Signature.R2", "Signature.R3", "Signature.U1", "Signature.U2"

#' Which signatures are present
#'
#' Determines how much of each signature is present in the sample given
#'
#' @param tumor.ref Either a data frame or location of input text file, where
#'   rows are samples, columns are trinucleotide contexts
#' @param sample.id Name of sample -- should be the rowname of tumor.ref.
#'   Default: `NULL`, means run deconstructSigs for all samples (rows)
#' @param signatures.ref Either a data frame or location of signature text file,
#'   where rows are signatures, columns are trinucleotide contexts. Default:
#'   NULL, means `signatures.nature2013`.
#' @param contexts.needed FALSE if tumor.file is a context file, TRUE if it is
#'   only mutation counts. Default: `TRUE`.
#' @param tri.counts.method Set to either:
#' \itemize{
#'  \item 'default' -- no further normalization
#'  \item 'exome' -- normalized by number of times each trinucleotide context is
#'   observed in the exome \item 'genome' -- normalized by number of times each
#'   trinucleotide context is observed in the genome
#' \item 'exome2genome' -- multiplied by a ratio of that trinucleotide's
#'   occurence in the genome to the trinucleotide's occurence in the exome
#' \item 'genome2exome' -- multiplied by a ratio of that trinucleotide's
#'   occurence in the exome to the trinucleotide's occurence in the genome}
#' @param genome.ref a reference genome [BSgenome] object to define
#' trinucleotide's occurence.
#' @param chr.list the chromose the use in the analysis. Default: `paste0("chr",
#' c(1:22, "X", "Y"))`.
#' @param exome.range a [GenomicRanges] object define the exome ranges.
#' @param associated Vector of associated signatures. If given, will narrow the
#'   signatures tested to only the ones listed.
#' @param signatures.limit Number of signatures to limit the search to
#' @param signature.cutoff Discard any signature contributions with a weight
#'   less than this amount
#' @return A list for every samples, the element of the list is the weights for
#'   each signatures, the product when those are multiplied on the signatures,
#'   the difference between the tumor sample and product, the tumor sample
#'   tricontext distribution given, and the unknown weight.
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
#' randomly.generated.tumors <- readRDS(
#'     system.file("extdata", "randomly.generated.tumors.rds",
#'                 package = "deconstructSigs")
#' )
#' test <- whichSignatures(
#'   tumor.ref = randomly.generated.tumors,
#'   sample.id = "2",
#'   contexts.needed = FALSE
#' )
#' @export
whichSignatures <- function(tumor.ref = NA, sample.id = NULL,
                            signatures.ref = NULL,
                            contexts.needed = TRUE,
                            tri.counts.method = "default",
                            genome.ref = NULL, chr.list = NULL,
                            exome.range = NULL,
                            associated = NULL,
                            signatures.limit = NA,
                            signature.cutoff = 0.06) {
  # check arguments --------------------------------------------------
  tri.counts.method <- match.arg(
    tri.counts.method,
    c("default", "genome", "exome", "exome2genome", "genome2exome")
  )
  if (is.null(signatures.ref)) {
    signatures.ref <- read_data("signatures.nature2013")
  }
  if (ncol(signatures.ref) == 78L && tri.counts.method != "default") {
    warning("Using default normalization for DBS signatures")
    tri.counts.method <- "default"
  }
  if (!identical(tri.counts.method, "default") && isTRUE(contexts.needed)) {
    if (!is.null(genome.ref) && !methods::is(genome.ref, "BSgenome")) {
      stop("genome.ref should be a `BSgenome` object or `NULL`")
    }
    if (is.null(chr.list)) chr.list <- paste0("chr", c(1:22, "X", "Y"))
    if (!is.null(exome.range) && !methods::is(exome.range, "GenomicRanges")) {
      stop("exome.range shoudl be a `GenomicRanges` object or `NULL`")
    }
  }
  # prepare data --------------------------------------------------------
  if (exists("tumor.ref", mode = "list")) {
    tumor <- tumor.ref
  } else if (is.character(tumor.ref)) {
    if (file.exists(tumor.ref)) {
      tumor <- utils::read.table(
        tumor.ref,
        sep = "\t", header = TRUE,
        as.is = TRUE, check.names = FALSE
      )
    } else {
      stop("tumor.ref is neither a file nor a loaded data frame", .call = FALSE)
    }
  } else {
    stop("Input tumor.ref needs to be a data frame or location of input text file", call. = FALSE)
  }

  if (isTRUE(contexts.needed)) {
    tumor <- getTriContextFraction(
      mut.counts = tumor,
      trimer.counts.method = tri.counts.method,
      genome.ref = genome.ref, chr.list = chr.list,
      exome.range = exome.range
    )
  }
  if (is.null(sample.id)) {
    sample.id <- rownames(tumor)
  } else {
    sample.id <- intersect(rownames(tumor), as.character(sample.id))
  }

  if (length(sample.id)) {
    # Take patient id given
    tumor <- as.matrix(tumor)
    res <- lapply(sample.id, function(x) {
      message("Processing sample ", x, appendLF = TRUE)
      sample_tumor <- tumor[x, , drop = FALSE]
      if (abs(rowSums(sample_tumor) - 1L) >= sqrt(.Machine$double.eps)) {
        stop(paste("Sample: ", x, " is not normalized\n", 'Consider using "contexts.needed = TRUE"', sep = " "))
      }
      deconstruct_sig_core(
        tumor = sample_tumor,
        signatures.ref = signatures.ref,
        associated = associated,
        signatures.limit = signatures.limit,
        signature.cutoff = signature.cutoff
      )
    })
    names(res) <- sample.id
    res
  } else {
    warning("No sample to run", call. = FALSE)
  }
}

deconstruct_sig_core <- function(tumor, signatures.ref,
                                 associated = NULL,
                                 signatures.limit = NA,
                                 signature.cutoff = 0.06) {
  # Read in Stratton signatures file
  if (exists("signatures.ref", mode = "list")) {
    signatures <- signatures.ref
  } else {
    if (file.exists(signatures.ref)) {
      signatures <- utils::read.csv(signatures.ref, sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE)
    } else {
      print("signatures.ref is neither a file nor a loaded data frame")
    }
  }

  signatures <- as.matrix(signatures)
  original.sigs <- signatures

  # Check column names are formatted correctly
  if (length(colnames(tumor)[colnames(tumor) %in% colnames(signatures)]) < length(colnames(signatures))) {
    colnames(tumor) <- changeColumnNames(colnames(tumor))
    if (length(colnames(tumor)[colnames(tumor) %in% colnames(signatures)]) < length(colnames(signatures))) {
      stop("Check column names on input file")
    }
  }

  # Ensure that columns in tumor match the order of those in signatures
  tumor <- tumor[, colnames(signatures), drop = FALSE]

  # Take a subset of the signatures
  if (is.character(associated)) {
    signatures <- signatures[
      rownames(signatures) %in% associated, ,
      drop = FALSE
    ]
  }

  if (is.na(signatures.limit)) {
    signatures.limit <- nrow(signatures)
  }

  # Remove signatures from possibilities if they have a "strong" peak not seen in the tumor sample
  zero.contexts <- colnames(tumor)[tumor < 0.01]
  corr.sigs <- which(
    signatures[, zero.contexts, drop = FALSE] >= 0.2,
    arr.ind = TRUE
  )
  signatures <- signatures[
    which(!rownames(signatures) %in% rownames(corr.sigs)), ,
    drop = FALSE
  ]

  # Set the weights matrix to 0
  weights <- matrix(0,
    nrow = nrow(tumor), ncol = nrow(signatures),
    dimnames = list(rownames(tumor), rownames(signatures))
  )

  seed <- findSeed(tumor, signatures)
  weights[seed] <- 1
  w <- weights * 10

  error_diff <- Inf
  error_threshold <- 1e-3

  num <- 0L
  while (error_diff > error_threshold) {
    num <- num + 1L
    # print(num)
    error_pre <- getError(tumor, signatures, w)
    if (error_pre == 0L) {
      break
    }
    # print(w)
    w <- updateW_GR(tumor, signatures, w, signatures.limit = signatures.limit)
    error_post <- getError(tumor, signatures, w)
    # print(paste("old_error = ", error_pre, sep = ""))
    # print(paste("new_error = ", error_post, sep = ""))
    # print(w)
    error_diff <- (error_pre - error_post) / error_pre
  }

  weights <- w / sum(w)
  unknown <- 0L

  ## filtering on a given threshold value (0.06 default)
  weights[weights < signature.cutoff] <- 0L
  unknown <- 1L - sum(weights)

  product <- weights %*% signatures
  diff <- tumor - product

  x <- matrix(
    data = 0L, nrow = 1L, ncol = nrow(original.sigs),
    dimnames = list(rownames(weights), rownames(original.sigs))
  )
  x <- data.frame(x)
  x[colnames(weights)] <- weights
  weights <- x

  out <- list(weights, tumor, product, diff, unknown)
  names(out) <- c("weights", "tumor", "product", "diff", "unknown")
  return(out)
}
