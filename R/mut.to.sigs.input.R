#' Returns complement of a sequence
#' 
#' Given an input sequence, returns either the reverse complement (default) or 
#' complement of the sequnce (nr = TRUE)
#' 
#' @keywords internal
#' @param tri Character vector of "A", "G", "C", or "T"
#' @param nr If TRUE returns only the complement of the input sequence
#' @return Returns a character vector containing the complemented input sequence
#' @export
findComp = function(tri, nr = FALSE) {
  pos=c("A","G","C","T")
  neg=c("T","C","G","A")
  tm=strsplit(tri,"")[[1]]
  out=c()
  for(i in 1:length(tm)) {
    if(tm[i] == "A"| tm[i] == "C"| tm[i] == "G"| tm[i] == "T"){
      ind=grep(tm[i],pos)
      out=c(out,neg[ind])
    }
    else{out=c(out,tm[i])} 
  }
  if(nr == FALSE){
    out=paste(rev(out),sep="",collapse="")
  }
  if(nr == TRUE){
    out=paste(out,sep="",collapse="")
  }
  return(out)
}


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
  
  # Add in context
  mut$context = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, mut[,chr], mut[,pos]-1, mut[,pos]+1, as.character = T)
  mut$mutcat = paste(mut[,ref], ">", mut[,alt], sep = "")
  
  if(any(substr(mut[,ref], 1, 1) != substr(mut[,'context'], 2, 2))){
    bad = mut[ which(substr(mut[,ref], 1, 1) != substr(mut[,'context'], 2, 2)), ]
    bad = paste(bad[,sample.id], bad[,chr], bad[,pos], bad[,ref], bad[,alt], sep = ':')
    bad = paste(bad, collapse = ',\ ')
    warning(paste('Check ref bases -- not all match context:\n ', bad, sep = ' '))
  }
  
  # Reverse complement the G's and A's
  mut$revmutcat = mut$mutcat
  mut$revcontext = mut$context
  
  gind = grep("G",substr(mut$mutcat,1,1))
  tind = grep("A",substr(mut$mutcat,1,1))
  
  if(length(gind) != 0){
    for(i in 1:length(gind)){
      mut$revmutcat[gind[i]]=findComp(mut$mutcat[gind[i]], nr = TRUE)
      mut$revcontext[gind[i]]=findComp(mut$context[gind[i]])
    }
  }
  
  if(length(tind) != 0){
    for(i in 1:length(tind)){
      mut$revmutcat[tind[i]]=findComp(mut$mutcat[tind[i]], nr = TRUE) 
      mut$revcontext[tind[i]]=findComp(mut$context[tind[i]])
    }
  }
  
  # Make the tricontext
  mut$tricontext = paste(substr(mut$revcontext, 1, 1), "[", mut$revmutcat, "]", substr(mut$revcontext, 3, 3), sep = "")
  
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
  
  final.df = data.frame(final.matrix, check.names = F)
  bad = names(which(rowSums(final.df) <= 50))
  if(length(bad) > 0){
    bad = paste(bad, collapse = ',\ ')
    warning(paste('Some samples have fewer than 50 mutations:\n ', bad, sep = ' '))
  }
  
  return(final.df)
  
}