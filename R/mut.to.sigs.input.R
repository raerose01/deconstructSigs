#' Converts mutation list to correct input format
#' 
#' Given a mutation list, outputs a data frame with counts of how frequently a
#' mutation is found within each trinucleotide context per sample ID.  Output
#' can be used as input into getTriContextFraction.
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
#' @return A data frame that contains sample IDs for the rows and trinucleotide
#'   contexts for the columns. Each entry is the count of how many times a
#'   mutation with that trinucleotide context is seen in the sample.
#' @examples
#' \dontrun{
#' sigs.input = mut.to.sigs.input(mut.ref = sample.mut.ref, 
#'                                sample.id = "Sample", 
#'                                chr = "chr", 
#'                                pos = "pos", 
#'                                ref = "ref", 
#'                                alt = "alt")
#'}
#' @export
mut.to.sigs.input = function(mut.ref, sample.id = 'Sample', chr = 'chr', pos = 'pos', ref = 'ref', alt = 'alt'){
  
  # print(paste("[", date(), "]", "Reading data"))
  if(exists("mut.ref", mode = "list")){
    mut.full <- mut.ref
  } else {
    if(file.exists(mut.ref)){
      mut.full <- utils::read.table(mut.ref, sep = "\t", header = TRUE, as.is = TRUE, check.names = FALSE)
    } else {
      stop("mut.ref is neither a file nor a loaded data frame")
    }
  }
  
  mut <- mut.full[,c(sample.id, chr, pos, ref, alt)]
  
  # Only take SNPs for now
  #mut.lengths <- with(mut, nchar(as.character(mut[,ref])))
  #mut.lengths <- with(mut, nchar(as.character(ref)))
  #mut <- mut[which(mut.lengths == 1),]
  mut$mut.lengths <- nchar(as.character(mut[,ref]))
  mut             <- mut[which(mut[,ref] %in% c('A', 'T', 'C', 'G') & mut[,alt] %in% c('A', 'T', 'C', 'G')),]
  
  # print(paste("[", date(), "]", "Adding in context"))
  # Add in context
  mut$context = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, mut[,chr], mut[,pos]-1, mut[,pos]+1, as.character = T)
  mut$mutcat = paste(mut[,ref], ">", mut[,alt], sep = "")
  
  # print(paste("[", date(), "]", "Checking reference bases"))
  if(any(substr(mut[,ref], 1, 1) != substr(mut[,'context'], 2, 2))){
    bad = mut[ which(substr(mut[,ref], 1, 1) != substr(mut[,'context'], 2, 2)), ]
    bad = paste(bad[,sample.id], bad[,chr], bad[,pos], bad[,ref], bad[,alt], sep = ':')
    bad = paste(bad, collapse = ',\ ')
    warning(paste('Check ref bases -- not all match context:\n ', bad, sep = ' '))
  }

  # print(paste("[", date(), "]", "Reverse complement the G's and A's - New algorithm"))
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

  # print(paste("[", date(), "]", "Make the tricontext"))
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
    tmp = subset(mut, mut[,sample.id] == i)
    beep = table(tmp$tricontext)
    for(l in 1:length(beep)){
      trimer = names(beep[l])
      if(trimer %in% all.tri){
        final.matrix[i, trimer] = beep[trimer]
      }
    }
  }
  # print(paste("[", date(), "]", "Done"))
  
  final.df = data.frame(final.matrix, check.names = F)
  bad = names(which(rowSums(final.df) <= 50))
  if(length(bad) > 0){
    bad = paste(bad, collapse = ',\ ')
    warning(paste('Some samples have fewer than 50 mutations:\n ', bad, sep = ' '))
  }
  
  return(final.df)
  
}