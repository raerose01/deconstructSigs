#' Converts a VCF file to correct input format
#' 
#' Given a VCF file, outputs a data frame with counts of how frequently a
#' mutation is found within each trinucleotide context per sample ID.  Output
#' can be used as input into getTriContextFraction.
#' 
#' The context sequence is taken from the BSgenome.Hsapiens.UCSC.hg19::Hsapiens
#' object, therefore the coordinates must correspond to the human hg19 assembly,
#' the UCSC version of the GRCh37 Homo sapiens assembly. This method will to
#' its best to translate chromosome names from other versions of the assembly
#' like NCBI or Ensembl. For instance, the following transformation will be
#' done: "1" -> "chr1"; "MT" -> "chrM"; "GL000245.1" -> "chrUn_gl000245"; etc.
#' 
#' This method relies on the VariantAnnotation package to read the VCF file.
#' 
#' @param vcf Location of the VCF file that is to be converted
#' @param bsg Only set if another genome build is required. Must be a BSgenome
#'   object.
#' @return A data frame that contains sample IDs for the rows and trinucleotide
#'   contexts for the columns. Each entry is the count of how many times a
#'   mutation with that trinucleotide context is seen in the sample.
#' @examples
#' \dontrun{
#' sigs.input = vcf.to.sigs.input(vcf = "variants.vcf")
#'}
#' @export
vcf.to.sigs.input <- function(vcf, bsg = NULL) {
  # Check dependency. Note that BiocGenerics and IRanges are dependencies of VariantAnnotation,
  # so no need to check for these.
  if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    stop("Package VariantAnnotation needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  gt <- VariantAnnotation::readGT(vcf, nucleotides = TRUE)

  # rownames are in the format chr:position_REF/ALT
  chr <- sub(":.+", "", rownames(gt))
  pos <- as.numeric(sub("_.+", "", sub(".+:", "", rownames(gt))))
  ref <- sub("[/|].+", "", sub(".+_", "", rownames(gt)))

  mut <- data.frame()
  for (sample in colnames(gt)) {
    a1 <- sub("[/|].+", "", gt[, sample])
    alt1 <- which(ref != a1)
    if (length(alt1) > 0) {
      mut <- rbind(mut, data.frame(sample = sample,
                                   chr = chr[alt1],
                                   pos = pos[alt1],
                                   ref = ref[alt1],
                                   alt = a1[alt1]))
    }  
    a2 <- sub(".+[/|]", "", gt[, sample])
    alt2 <- which(ref != a2 & a1 != a2)
    if (length(alt2) > 0) {
      mut <- rbind(mut, data.frame(sample = sample,
                                   chr = chr[alt2],
                                   pos = pos[alt2],
                                   ref = ref[alt2],
                                   alt = a2[alt2]))
    }
  }
  
  sigs <- mut.to.sigs.input(mut.ref = mut, sample.id = "sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg = bsg)

  return(sigs)
}
