#' Formats the data for use in plotSignatures()
#' 
#' Makes the contexts useable in plotSignatures()
#' 
#' @keywords internal
#' @param contexts One of the entries in the output list from whichSignatures()
#' @return Returns a data frame with sample.id, full_context, fraction, and 
#'   mutation as column names

#' @export
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
#' @param sub A character vector that specifies cancer subtype for plot title,
#'   if wanted
#' @return Plots the trinucleotide frequency for the given tumor on the top 
#'   panel, the calculated one on the middle panel, and the difference between 
#'   the two on the bottom panel.
#' @export
#' @examples
#' plotSignatures(example.output, sub = "example")

plotSignatures = function(sigs.output, sub = ""){
  
  op <- graphics::par()
  
  tumor   <- sigs.output[["tumor"]]
  product <- sigs.output[["product"]]
  diff    <- sigs.output[["diff"]]
  weights <- sigs.output[["weights"]]
  
  y_limit        <- 1.2 * max(tumor, product)
  tumor_plotting <- formatContexts(tumor)
  product_plotting <- formatContexts(product)
  error_summed <- round(sqrt(sum(diff*diff)), digits = 3)
  diff_plotting <- formatContexts(diff)
  
  name           <- unique(tumor_plotting$sample.id)
  subtype        <- sub
  tmp <- which(weights != 0)
  c <- paste(colnames(weights)[tmp[1]], " : ", round(weights[tmp[1]], 3), sep = "")
  if(length(tmp)>1){
    for(i in tmp[2:length(tmp)]){
      c <- paste(c, " & ", colnames(weights)[i], " : ", round(weights[i], 3), sep = "")
    }
  }
  
  if(subtype == ''){
    top.title <- name
  }
  if(subtype != ''){
    top.title <- paste(name, " -- ", subtype, sep = "")
  }
  
  graphics::par(mfrow = c(3,1), xpd = FALSE, mar=c(5, 4, 4, 2), oma=c(0,0,0,4))
  grDevices::palette(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))
  
  graphics::barplot(tumor_plotting$fraction, names.arg = tumor_plotting$full_context, cex.names = 0.7, las = 2, col = NA, ylim = c(0, y_limit), border = NA, xaxt='n', ann=FALSE, yaxt = 'n', space = 0.3)
  graphics::box()
  x = graphics::par('usr')
  graphics::abline(h = seq(from = 0, to = y_limit, by = 0.01), col = '#d3d3d350', lty = 1)
  graphics::abline(v = seq(from = x[1], to = x[2], by = 1), col = '#d3d3d350', lty = 1)
  graphics::barplot(tumor_plotting$fraction, names.arg = tumor_plotting$full_context, cex.names = 0.7, las = 2, col = factor(tumor_plotting$mutation), ylim = c(0, y_limit), border = NA, space = 0.3, main = top.title, ylab = 'fraction', add = TRUE)
  
  graphics::barplot(product_plotting$fraction, names.arg = product_plotting$full_context, cex.names = 0.7, las = 2, col = NA, ylim = c(0, y_limit), border = NA, xaxt='n', ann=FALSE, yaxt = 'n', space = 0.3)
  graphics::box()
  x = graphics::par('usr')
  graphics::abline(h = seq(from = 0, to = y_limit, by = 0.01), col = '#d3d3d350', lty = 1)
  graphics::abline(v = seq(from = x[1], to = x[2], by = 1), col = '#d3d3d350', lty = 1)
  graphics::barplot(product_plotting$fraction, names.arg = product_plotting$full_context, cex.names = 0.7, las = 2, col = factor(product_plotting$mutation), ylim = c(0, y_limit), border = NA, space = 0.3, main = c, ylab = 'fraction', add = TRUE)
  
  graphics::barplot(diff_plotting$fraction, names.arg = diff_plotting$full_context, cex.names = 0.7, las = 2, col = NA, ylim = c(-y_limit, y_limit), border = NA, xaxt='n', ann=FALSE, yaxt = 'n', space = 0.3)
  graphics::box()
  x = graphics::par('usr')
  graphics::abline(h = seq(from = -y_limit, to = y_limit, by = 0.02), col = '#d3d3d350', lty = 1)
  graphics::abline(v = seq(from = x[1], to = x[2], by = 1), col = '#d3d3d350', lty = 1)
  graphics::barplot(diff_plotting$fraction, names.arg = diff_plotting$full_context, cex.names = 0.7, las = 2, col = factor(diff_plotting$mutation), ylim = c(-y_limit, y_limit), border = 'black', space = 0.3, main = paste("error = ", error_summed, sep = ""), ylab = 'fraction', add = TRUE)
  
  graphics::par(fig=c(0,1,0,1), oma = c(1, 1, 1, 1), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  graphics::legend('right', legend = unique(tumor_plotting$mutation), col = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"), bty = 'n', ncol = 1, inset=c(-0,0), pch = 15, xpd = TRUE, pt.cex = 2.5)
  
  on.exit(suppressWarnings(graphics::par(op)))
  
}

#' Plots a tumor profile
#' 
#' Uses the output from whichSignatures() and creates a plot of the given tumor 
#' mutational spectrum and the calculated one
#' 
#' @param tumor A tumor matrix
#' @param sub A character vector that specifies cancer subtype for plot title,
#'   if wanted
#' @return Plots the trinucleotide frequency for the given tumor
#' @export
#' @examples
#' plotTumor(example.output[['tumor']], sub = "example")

plotTumor = function(tumor, sub = ""){
  
  op <- graphics::par()
  
  tumor          <- as.matrix(tumor)
  
  y_limit        <- 1.2 * max(tumor)
  tumor_plotting <- formatContexts(tumor)
  
  name           <- unique(tumor_plotting$sample.id)
  subtype        <- sub
  
  if(subtype == ''){
    top.title <- name
  }
  if(subtype != ''){
    top.title <- paste(name, " -- ", subtype, sep = "")
  }
  
  grDevices::palette(c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))
  
  graphics::barplot(tumor_plotting$fraction, names.arg = tumor_plotting$full_context, cex.names = 0.7, las = 2, col = NA, ylim = c(0, y_limit), border = NA, xaxt='n', ann=FALSE, yaxt = 'n', space = 0.3)
  graphics::box()
  x = graphics::par('usr')
  graphics::abline(h = seq(from = 0, to = y_limit, by = 0.01), col = '#d3d3d350', lty = 1)
  graphics::abline(v = seq(from = x[1], to = x[2], by = 1), col = '#d3d3d350', lty = 1)
  graphics::barplot(tumor_plotting$fraction, names.arg = tumor_plotting$full_context, cex.names = 0.7, las = 2, col = factor(tumor_plotting$mutation), ylim = c(0, y_limit), border = NA, space = 0.3, main = top.title, ylab = 'fraction', add = TRUE)
  
  graphics::legend('topright', legend = unique(tumor_plotting$mutation), col = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"), bty = 'n', ncol = 6, inset=c(-0,0), pch = 15, xpd = TRUE, pt.cex = 2.5)
  
  #graphics::par(fig=c(0,1,0,1), oma = c(1, 1, 1, 1), mar = c(0, 0, 0, 0), new = TRUE)
  #graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  #graphics::legend('topright', legend = unique(tumor_plotting$mutation), col = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"), bty = 'n', ncol = 1, inset=c(-0,0), pch = 15, xpd = TRUE, pt.cex = 2.5)
  
  on.exit(suppressWarnings(graphics::par(op)))

}

#' Plots the weights from whichSignatures()
#' 
#' Uses the output from whichSignatures() and creates a pie chart of the weights
#' outputted
#' 
#' @param sigs.output The list output from whichSignatures()
#' @param sub A character vector that specifies cancer subtype for plot title, 
#'   if wanted
#' @param add.color Optional character vector to assign additional colors for
#'   novel signatures
#' @return Plots a pie chart of the weights calculated in the given tumor sample
#' @export
#' @examples
#' makePie(example.output)
makePie <- function(sigs.output, sub = "", add.color = NULL){
  
  weights              <- data.frame(sigs.output[["weights"]])
  unknown              <- sigs.output[["unknown"]]
  
  weights$unknown      <- unknown
  
  # Set up all the colors possible
  all.sigs             <- c("Signature.1", "Signature.1A", "Signature.1B", "Signature.2",  "Signature.3", "Signature.4",  
                            "Signature.5", "Signature.6", "Signature.7", "Signature.8", "Signature.9",  
                            "Signature.10", "Signature.11", "Signature.12", "Signature.13", "Signature.14", 
                            "Signature.15", "Signature.16", "Signature.17", "Signature.18", "Signature.19", 
                            "Signature.20", "Signature.21", "Signature.R1", "Signature.R2", "Signature.R3",
                            "Signature.U1", "Signature.U2","Signature.22", "Signature.23", "Signature.24",
                            "Signature.25", "Signature.26", "Signature.27", "Signature.28", "Signature.29", "Signature.30",
                            "unknown")

  all.colors           <- c("#023FA5", "#023FA5","#7D87B9","#BEC1D4","#D6BCC0","#BB7784",
                            "gold","#4A6FE3","#8595E1","#B5BBE3","#E6AFB9",
                            "#E07B91","#D33F6A","#11C638","#8DD593","#C6DEC7",
                            "#EAD3C6","#F0B98D","#EF9708","#0FCFC0","#9CDED6",
                            "#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4","#866097",
                            "#008941","#A30059","#F6C4E1","#F79CD4","#866097",
                            "#008941","#A30059","#008080","#8B0000","#F4A460","#663399",
                            "#706563")

  
  all.colors           <- cbind(as.data.frame(all.sigs), as.data.frame(all.colors))
  colnames(all.colors) <- c("signature", "color")
  
  # Only include weights that are not 0
  ind                  <- which(weights > 0)
  weights              <- weights[,ind]
  
  missing.colors       <- colnames(weights)[!colnames(weights) %in% all.colors$signature]
  if(!is.null(add.color)){
    tmp <- data.frame(missing.colors, add.color)
    colnames(tmp) = c('signature', 'color')
    all.colors <- rbind(all.colors, tmp)
  }
  
  # Set up color palette
  sigs.present         <- colnames(weights)
  missing.colors       <- sigs.present[!sigs.present %in% all.colors$signature]
  if(length(missing.colors) > 0){
    warning(paste('No color assigned for: \n',  paste(missing.colors, collapse = ',\ '), '.\nTo assign one, use add.color parameter.', sep = ''))
  }
  colors.sigs.present = all.colors$color[match(sigs.present, all.colors$signature)]
  grDevices::palette(as.character(colors.sigs.present))
  
  if(sub == ''){
    top.title <- rownames(weights)
  }
  if(sub != ''){
    top.title <- paste(rownames(weights), " -- ", sub, sep = "")
  }
  
  graphics::pie(t(weights), col = factor(colnames(weights), levels = unique(colnames(weights))), labels = colnames(weights), main = top.title)

}
################################################ 




