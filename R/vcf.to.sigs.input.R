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
#' @return A data frame that contains sample IDs for the rows and trinucleotide
#'   contexts for the columns. Each entry is the count of how many times a
#'   mutation with that trinucleotide context is seen in the sample.
#' @examples
#' \dontrun{
#' sigs.input = vcf.to.sigs.input(vcf = "variants.vcf")
#'}
#' @export
vcf.to.sigs.input <- function(vcf) {
  # Check dependency. Note that BiocGenerics and IRanges are dependencies of VariantAnnotation,
  # so no need to check for these.
  if (!requireNamespace("VariantAnnotation", quietly = TRUE)) {
    stop("Package VariantAnnotation needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  genome <- "hg19"
  vcf.data <- VariantAnnotation::readVcf(vcf, genome)
  
  mut <- data.frame()
  for (sample in VariantAnnotation::samples(VariantAnnotation::header(vcf.data))) {
    alleles <- sort(unique(unlist(strsplit(as.matrix(VariantAnnotation::geno(vcf.data)$GT[, sample]), "[|/]"))))
    alleles <- as.numeric(grep("0", alleles, value = T, invert = T))
    for (allele in alleles) {
      has.alt.allele = grep(allele, as.matrix(VariantAnnotation::geno(vcf.data)$GT[, sample]))
      n <- length(has.alt.allele)
      if (n == 0) {
        next;
      }
      alt.array <- sapply(VariantAnnotation::alt(vcf.data)[has.alt.allele],
                          FUN = function(x) { return(as.character(x[[allele]])) } )
      mut <- rbind(mut, data.frame(sample = sample,
                                   chr = GenomeInfoDb::seqnames(vcf.data)[has.alt.allele],
                                   pos = BiocGenerics::start(IRanges::ranges(vcf.data))[has.alt.allele],
                                   ref = as.character(VariantAnnotation::ref(vcf.data))[has.alt.allele],
                                   alt = alt.array
      ))
    }
  }
  print(mut)
  return(mut.to.sigs.input(mut.ref = mut, sample.id = "sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt"))
}
