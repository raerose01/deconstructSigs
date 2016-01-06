#' Formats the data for use in plotSignatures()
#' 
#' Makes the contexts useable in plotSignatures()
#' 
#' @keywords internal
#' @param contexts One of the entries in the output list from whichSignatures()
#' @return Returns a data frame with sample.id, full_context, fraction, and mutation as column names
formatContexts = function(contexts){
  for_plotting = reshape2::melt(contexts)
  colnames(for_plotting)[1] = "sample.id"
  colnames(for_plotting)[2] = "full_context"
  colnames(for_plotting)[3] = "fraction"
  for_plotting$mutation = substr(for_plotting$full_context, 3,5)
  for_plotting$full_context = as.factor(for_plotting$full_context)
  for_plotting$full_context = factor(for_plotting$full_context, levels = for_plotting$full_context, ordered=T)
  return(for_plotting)
}

#' Plots the result from whichSignatures()
#' 
#' Uses the output from whichSignatures() and creates a plot of the given tumor 
#' mutational spectrum and the calculated one
#' 
#' @param sigs.output The list output from whichSignatures()
#' @param sub A character vector that specifies cancer subtype if wanted
#' @return Plots the trinucleotide frequency for the given tumor on the top 
#'   panel, the calculated one on the middle panel, and the difference between 
#'   the two on the bottom panel.
#' @examples
#' plotSignatures(example.output, sub = "example")
plotSignatures = function(sigs.output, sub = ""){
  
  tumor   <- sigs.output[["tumor"]]
  product <- sigs.output[["product"]]
  diff    <- sigs.output[["diff"]]
  weights <- sigs.output[["weights"]]
  
  y_limit        <- 1.2 * max(tumor, product)
  tumor_plotting <- formatContexts(tumor)
  name           <- unique(tumor_plotting$sample.id)
  subtype        <- sub
  
  #ggp <- ggplot2::ggplot(tumor_plotting, ggplot2::aes(x=full_context, y=fraction, fill=mutation))
  ggp <- ggplot2::ggplot(data = tumor_plotting, ggplot2::aes_string(x = 'full_context', y = 'fraction', fill = 'mutation'))
  #plot1 = ggp + ggplot2::geom_histogram(stat="identity") + ggplot2::ggtitle(paste(name, " -- ", subtype, sep = "")) + ggplot2::scale_fill_manual(values=c("#1EBBEB", "#000000", "#E0242A", "#A3A3A3", "#A3CC6E", "#EDC6C2")) + ggplot2::theme_bw() +  ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, vjust = 0.5), plot.title = ggplot2::element_text(lineheight=.8, face="bold")) + ggplot2::ylim(0, y_limit)
  plot1 <- ggp + ggplot2::geom_histogram(stat="identity") + ggplot2::ggtitle(paste(name, " -- ", subtype, sep = "")) + ggplot2::scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + ggplot2::theme_bw() +  ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, vjust = 0.5), plot.title = ggplot2::element_text(lineheight=.8, face="bold")) + ggplot2::ylim(0, y_limit)
  
  tmp <- which(weights != 0)
  c <- paste(colnames(weights)[tmp[1]], " : ", round(weights[tmp[1]], 3), sep = "")
  if(length(tmp)>1){
    for(i in tmp[2:length(tmp)]){
      c = paste(c, " & ", colnames(weights)[i], " : ", round(weights[i], 3), sep = "")
    }
  }
  
  product_plotting <- formatContexts(product)
  #ggp <- ggplot2::ggplot(product_plotting, ggplot2::aes(x=full_context, y=fraction, fill=mutation))
  ggp <- ggplot2::ggplot(product_plotting, ggplot2::aes_string(x='full_context', y='fraction', fill='mutation'))
  #plot2 = ggp + ggplot2::geom_histogram(stat="identity") + ggplot2::ggtitle(c) + ggplot2::scale_fill_manual(values=c("#1EBBEB", "#000000", "#E0242A", "#A3A3A3", "#A3CC6E", "#EDC6C2")) + ggplot2::theme_bw() +  ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, vjust = 0.5), plot.title = ggplot2::element_text(lineheight=.8, face="bold")) + ggplot2::ylim(0, y_limit)
  plot2 <- ggp + ggplot2::geom_histogram(stat="identity") + ggplot2::ggtitle(c) + ggplot2::scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + ggplot2::theme_bw() +  ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, vjust = 0.5), plot.title = ggplot2::element_text(lineheight=.8, face="bold")) + ggplot2::ylim(0, y_limit)
  
  error_summed <- sqrt(sum(diff*diff))
  diff_plotting <- formatContexts(diff)
  #ggp <- ggplot2::ggplot(diff_plotting, ggplot2::aes(x=full_context, y=fraction, fill=mutation))
  ggp <- ggplot2::ggplot(diff_plotting, ggplot2::aes_string(x='full_context', y='fraction', fill='mutation'))
  #plot3 = ggp + ggplot2::geom_histogram(stat = "identity", position="identity") + ggplot2::ggtitle(paste("error = ", error_summed, sep = "")) + ggplot2::scale_fill_manual(values=c("#1EBBEB", "#000000", "#E0242A", "#A3A3A3", "#A3CC6E", "#EDC6C2")) + ggplot2::theme_bw() +  ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, vjust = 0.5), plot.title = ggplot2::element_text(lineheight=.8, face="bold")) + ggplot2::ylim(-y_limit, y_limit)
  plot3 <- ggp + ggplot2::geom_histogram(stat = "identity", position="identity") + ggplot2::ggtitle(paste("error = ", error_summed, sep = "")) + ggplot2::scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")) + ggplot2::theme_bw() +  ggplot2::theme(axis.text.x=ggplot2::element_text(angle = 90, vjust = 0.5), plot.title = ggplot2::element_text(lineheight=.8, face="bold")) + ggplot2::ylim(-y_limit, y_limit)
  
  gA <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot1))
  gB <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot2))
  gC <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(plot3))
  #grid.newpage()
  gridExtra::grid.arrange(gA, gB, gC, nrow = 3, heights = c(1/3, 1/3, 1/3))
  
}

#' Plots the weights from whichSignatures()
#' 
#' Uses the output from whichSignatures() and creates a pie chart
#' of the weights outputted
#' 
#' @param sigs.output The list output from whichSignatures()
#' @param sub A character vector that specifies cancer subtype if wanted
#' @return Plots a pie chart of the weights calculated in the given
#' tumor sample
#' @examples
#' makePie(example.output)
makePie <- function(sigs.output, sub = ""){
  
  weights              <- data.frame(sigs.output[["weights"]])
  unknown              <- sigs.output[["unknown"]]
  
  weights$unknown      <- unknown
  
  # Set up all the colors possible
  #all.sigs             <- colnames(weights)
  all.sigs             <- c("Signature.1A", "Signature.1B", "Signature.2",  "Signature.3", "Signature.4",  
                            "Signature.5", "Signature.6", "Signature.7", "Signature.8", "Signature.9",  
                            "Signature.10", "Signature.11", "Signature.12", "Signature.13", "Signature.14", 
                            "Signature.15", "Signature.16", "Signature.17", "Signature.18", "Signature.19", 
                            "Signature.20", "Signature.21", "Signature.R1", "Signature.R2", "Signature.R3",
                            "Signature.U1", "Signature.U2", "unknown")

  all.colors           <- c("#023FA5","#7D87B9","#BEC1D4","#D6BCC0","#BB7784",
                            "#FFFFCC","#4A6FE3","#8595E1","#B5BBE3","#E6AFB9",
                            "#E07B91","#D33F6A","#11C638","#8DD593","#C6DEC7",
                            "#EAD3C6","#F0B98D","#EF9708","#0FCFC0","#9CDED6",
                            "#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4","#866097",
                            "#008941","#A30059","#706563")
  
  all.colors           <- cbind(as.data.frame(all.sigs), as.data.frame(all.colors))
  colnames(all.colors) <- c("signature", "color")
  
  # Only include weights that are not 0
  ind                  <- which(weights > 0)
  weights              <- weights[,ind]
  
  # Set up color palette
  sigs.present         <- colnames(weights)
  colors.sigs.present = all.colors$color[match(sigs.present, all.colors$signature)]
  palette(as.character(colors.sigs.present))  
  
  pie(t(weights), col = factor(colnames(weights), levels = unique(colnames(weights))), labels = colnames(weights), main = paste(rownames(weights), " -- ", sub, sep = ""))

}
################################################ 




