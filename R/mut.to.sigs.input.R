#' Converts mutation list to correct input format
#' 
#' Given a mutation list, outputs a data frame with counts of how frequently a 
#' mutation is found within each trinucleotide context per sample ID.  Output 
#' can be used as input into getTriContextFraction.
#' 
#' The context sequence is taken from the BSgenome.Hsapiens.UCSC.hgX::Hsapiens 
#' object. Therefore the coordinates must correspond to the human hgX assembly. 
#' Default is set to the UCSC hg19 assembly, which corresponds to the GRCh37 
#' assembly. If another assembly is required, it must already be present in the 
#' R workspace and fed as a parameter. This method will translate chromosome 
#' names from other versions of the assembly like NCBI or Ensembl. For instance,
#' the following transformation will be done: "1" -> "chr1"; "MT" -> "chrM"; 
#' "GL000245.1" -> "chrUn_gl000245"; etc.
#' 
#' @param mut.ref Location of the mutation file that is to be converted or name 
#'   of data frame in environment
#' @param sample.id Column name in the mutation file corresponding to the Sample
#'   ID
#' @param chr Column name in the mutation file corresponding to the chromosome
#' @param pos Column name in the mutation file corresponding to the mutation 
#'   position
#' @param ref Column name in the mutation file corresponding to the reference 
#'   base
#' @param alt Column name in the mutation file corresponding to the alternate 
#'   base
#' @param bsg Only set if another genome build is required. Must be a BSgenome
#'   object.
#' @param sig.type Are SBS or DBS signatures being used?
#' @param dbs_table Possible DBS changes and their reverse complements
#' @return A data frame that contains sample IDs for the rows and trinucleotide 
#'   contexts for the columns. Each entry is the count of how many times a 
#'   mutation with that trinucleotide context is seen in the sample.
#' @examples
#' \dontrun{
#'sigs.input = mut.to.sigs.input(mut.ref = sample.mut.ref, sample.id = "Sample",
#'chr = "chr", pos = "pos", ref = "ref", alt = "alt", bsg = 
#'BSgenome.Hsapiens.UCSC.hg19, sig.type = 'SBS', dbs_table = dbs_possible)
#'}
#' @export
mut.to.sigs.input = function(mut.ref, sample.id = 'Sample', chr = 'chr', pos = 'pos', ref = 'ref', alt = 'alt', bsg = NULL, sig.type = 'SBS', dbs_table = dbs_possible){
  
  if(exists("mut.ref", mode = "list")){
    mut.full <- mut.ref
  } else {
    if(file.exists(mut.ref)){
      mut.full <- utils::read.table(mut.ref, sep = "\t", header = TRUE, as.is = FALSE, check.names = FALSE)
    } else {
      stop("mut.ref is neither a file nor a loaded data frame")
    }
  }
  
  mut <- mut.full[,c(sample.id, chr, pos, ref, alt)]
  
  # don't need trinucleotide context if looking for DBS
  if(sig.type == 'DBS'){
    mut[,ref] <- as.character(mut[,ref])
    mut[,alt] <- as.character(mut[,alt])
    mut <- mut[which(nchar(mut[, ref]) ==2 & nchar(mut[,alt]) == 2),]
    mut$dbs <- paste(mut[,ref], mut[,alt], sep = '>') 
    mut$dbs_condensed <- dbs_possible$dbs_condensed[match(mut$dbs, dbs_possible$dbs)]
    final.df <- as.data.frame.matrix(table(mut[,sample.id], factor(mut$dbs_condensed, levels = unique(dbs_possible$dbs_condensed))), stringsAsFactors = FALSE)
  }
  
  # And now look at SBS
  #mut.lengths <- with(mut, nchar(as.character(mut[,ref])))
  #mut.lengths <- with(mut, nchar(as.character(ref)))
  #mut <- mut[which(mut.lengths == 1),]
  #mut$mut.lengths <- nchar(as.character(mut[, ref]))
  if(sig.type == 'SBS'){
    mut <- mut[which(mut[, ref] %in% c('A', 'T', 'C', 'G') & mut[, alt] %in% c('A', 'T', 'C', 'G')),]
  
    # Fix the chromosome names (in case they come from Ensembl instead of UCSC)
    mut[, chr] <- factor(mut[, chr])
    levels(mut[, chr]) <- sub("^([0-9XY])", "chr\\1", levels(mut[, chr]))
    levels(mut[, chr]) <- sub("^MT", "chrM", levels(mut[, chr]))
    levels(mut[, chr]) <- sub("^(GL[0-9]+).[0-9]", "chrUn_\\L\\1", levels(mut[, chr]), perl = T)
  
    # Check the genome version the user wants to use
    # If set to default, carry on happily
    if(is.null(bsg)){
      # Remove any entry in chromosomes that do not exist in the BSgenome.Hsapiens.UCSC.hg19::Hsapiens object    
      unknown.regions <- levels(mut[, chr])[which(!(levels(mut[, chr]) %in% GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)))]
      if (length(unknown.regions) > 0) {
        unknown.regions <- paste(unknown.regions, collapse = ',\ ')
        warning(paste('Check chr names -- not all match BSgenome.Hsapiens.UCSC.hg19::Hsapiens object:\n', unknown.regions, sep = ' '))      
        mut <- mut[mut[, chr] %in% GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens), ]
      }
      # Add in context
      mut$context = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, mut[,chr], mut[,pos]-1, mut[,pos]+1, as.character = T)
    }
    
    # If set to another build, use that one 
    # bsg parameter should be BSgenome object already
    if(!is.null(bsg)){
      if(class(bsg) != 'BSgenome'){
        stop('The bsg parameter needs to either be set to default or a BSgenome object.')
      }
      # Remove any entry in chromosomes that do not exist in the BSgenome.Hsapiens.UCSC.hgX::Hsapiens object
      unknown.regions <- levels(mut[, chr])[which(!(levels(mut[, chr]) %in% GenomeInfoDb::seqnames(bsg)))]
      if (length(unknown.regions) > 0) {
        unknown.regions <- paste(unknown.regions, collapse = ',\ ')
        warning(paste('Check chr names -- not all match',attr(bsg, which = 'pkgname'),'object:\n', unknown.regions, sep = ' '))
        mut <- mut[mut[, chr] %in% GenomeInfoDb::seqnames(bsg), ]
      }
      # Add in context
      mut$context = BSgenome::getSeq(bsg, mut[,chr], mut[,pos]-1, mut[,pos]+1, as.character = T) 
    }
   
    mut$mutcat = paste(mut[,ref], ">", mut[,alt], sep = "")
    
    if(any(substr(mut[,ref], 1, 1) != substr(mut[,'context'], 2, 2))){
      bad = mut[ which(substr(mut[,ref], 1, 1) != substr(mut[,'context'], 2, 2)), ]
      bad = paste(bad[,sample.id], bad[,chr], bad[,pos], bad[,ref], bad[,alt], sep = ':')
      bad = paste(bad, collapse = ',\ ')
      warning(paste('Check ref bases -- not all match context:\n ', bad, sep = ' '))
    }
  
    # Reverse complement the G's and A's
    gind = grep("G",substr(mut$mutcat,1,1))
    tind = grep("A",substr(mut$mutcat,1,1))
    
    mut$std.mutcat = mut$mutcat
    mut$std.mutcat[c(gind, tind)] <- gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", mut$std.mutcat[c(gind, tind)])))) # to lowercase
    mut$std.mutcat[c(gind, tind)] <- gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", mut$std.mutcat[c(gind, tind)])))) # complement
  
    mut$std.context = mut$context
    mut$std.context[c(gind, tind)] <- gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", mut$std.context[c(gind, tind)])))) # to lowercase
    mut$std.context[c(gind, tind)] <- gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", mut$std.context[c(gind, tind)])))) # complement
    mut$std.context[c(gind, tind)] <- sapply(strsplit(mut$std.context[c(gind, tind)], split = ""), function(str) {paste(rev(str), collapse = "")}) # reverse
  
    # Make the tricontext
    mut$tricontext = paste(substr(mut$std.context, 1, 1), "[", mut$std.mutcat, "]", substr(mut$std.context, 3, 3), sep = "")
    
    # Generate all possible trinucleotide contexts
    all.tri = c()
    for(i in c("A", "C", "G", "T")){
      for(j in c("C", "T")){
        for(k in c("A", "C", "G", "T")){
          if(j != k){
            for(l in c("A", "C", "G", "T")){
              tmp = paste(i, "[", j, ">", k, "]", l, sep = "")
              all.tri = c(all.tri, tmp)
            }
          }
        }
      }
    }
    
    all.tri <- all.tri[order(substr(all.tri, 3, 5))]
    
    final.matrix = matrix(0, ncol = 96, nrow = length(unique(mut[,sample.id])))
    colnames(final.matrix) = all.tri
    rownames(final.matrix) = unique(mut[,sample.id])
    
    # print(paste("[", date(), "]", "Fill in the context matrix"))
    for(i in unique(mut[,sample.id])){
      tmp = mut[which(mut[,sample.id] == i),]
      beep = table(tmp$tricontext)
      for(l in 1:length(beep)){
        trimer = names(beep[l])
        if(trimer %in% all.tri){
          final.matrix[i, trimer] = beep[trimer]
        }
      }
    }
    
    final.df = data.frame(final.matrix, check.names = FALSE)
  }
  
  bad = names(which(rowSums(final.df) <= 50))
  if(length(bad) > 0){
    bad = paste(bad, collapse = ',\ ')
    warning(paste('Some samples have fewer than 50 mutations:\n ', bad, sep = ' '))
  }
  
  return(final.df)
  
}
