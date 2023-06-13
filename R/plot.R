
#' Plot
#'
#' @param dataset gwasglue2 DataSet object
#' @param type Type of plot (Only available "manhattan" plots at the moment)
#' @param title Main title for the plot
#' @export
#'
#' @return A plot
plot_gwasglue <- function(dataset, type, title){
  
  if(type == "manhattan"){
    ntraits <- getLength(dataset)
    nb_rows <- ceiling(ntraits/2)
     # Add main title
    
    
    par(mfrow=c(nb_rows, 2))
    
    
      for (i in 1:ntraits){
        plot(dataset@summary_sets[[i]]@ss$position, -log10(dataset@summary_sets[[i]]@ss$p), main = "", xlab = "position", ylab = "-log10(p-value)", cex=0.8, pch=20)
        mtext(dataset@summary_sets[[i]]@metadata$trait, side = 3, line = 0.5)
      }
    mtext(as.expression(bquote(bold(.(title)))),
          side = 3,
          line = - 2.5,
          outer = TRUE,
          cex = 1.3)
  }
}


