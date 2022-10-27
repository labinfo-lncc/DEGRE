#' Create a volcano plot to visualize the proportion of downregulated and upregulated genes by applying a log2FC cutoff.
#'
#' @param results a data.frame object. It receives the output of the DEGRE function, filtered or not, as input.
#' @param log2FC_cutoff it stores the cutoff of the log2FoldChange. The default is 1.
#' @param padj it stores the cutoff of the P-adjusted value (Q-value). The default is 0.05.
#' @param font.x the font size of the x axis. The default is 10.
#' @param font.y the font size of the y axis. The default is 10.
#' @param font.tickslab the font size of the ticks lab. The default is 10.
#' @param downregulated_color the colors of the downregulated genes. The default is "coral2".
#' @param upregulated_color the colors of the upregulated genes. The default is "cornflowerblue".
#' @param xlab the x lab text. The default is "log2Foldchange".
#' @param ylab the y lab text. The default is "-log10(P-value)".
#' @param legend_position you need to specify here the position of the legend. The default is "right".
#' @param legend.title the title of the legend. The default is "Regulation".
#'
#' @examples
#' # Reading the count matrix and the design matrix for an example:
#' dir <- system.file("data", package = "DEGRE")
#' tab <- read.csv(file.path(dir,"count_matrix_for_example.csv"))
#' row.names(tab) <- tab[,1]; tab <- tab[,-1]
#' des <- read.csv(file.path(dir,"design_matrix_for_example.csv"))
#' # Running DEGRE function:
#' results <- DEGRE(count_matrix = tab,
#'                p_value_adjustment = "BH",
#'                design_matrix = des,
#'                formula = "condition + (1|sex)")
#' # Running the VolcanoDEGRE function
#' VolcanoDEGRE(results = results,
#'           log2FC_cutoff = 1,
#'           padj = 0.05,
#'           font.x = 10,
#'           font.y = 10,
#'           font.tickslab = 10,
#'           downregulated_color = "coral2",
#'           upregulated_color = "cornflowerblue",
#'           xlab = "log2Foldchange",
#'           ylab = "-log10(P-value)",
#'           legend_position = "right",
#'           legend.title = "Regulation")
#'
#' @export
VolcanoDEGRE <- function (results,
                          log2FC_cutoff = 1,
                          padj = 0.05,
                          font.x = 10,
                          font.y = 10,
                          font.tickslab = 10,
                          downregulated_color = "coral2",
                          upregulated_color = "cornflowerblue",
                          xlab = "log2Foldchange",
                          ylab = "-log10(P-value)",
                          legend_position = "right",
                          legend.title = "Regulation")
{

  if(missing(results))
    stop("You need to enter with the output of the DEGRE function.")

  results$diffexpressed <- "NO"
  results$diffexpressed[results$log2FC > log2FC_cutoff & results$`Q-value` < padj] <- "UP"
  results$diffexpressed[results$log2FC < -log2FC_cutoff & results$`Q-value` < padj] <- "DOWN"

  graph <- ggplot(data=results, aes(x=results$log2FC, y=-log10(results$`Q-value`), col=results$diffexpressed)) +
    geom_point() +
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values=c(paste0(downregulated_color), "grey", paste0(upregulated_color))) +
    geom_vline(xintercept=c(-log2FC_cutoff, log2FC_cutoff), col="red") +
    geom_hline(yintercept=-log10(padj), col="red")

  ggpar(graph,
        font.x = c(paste(font.x),"black"),
        font.y = c(paste(font.y),"black"),
        font.tickslab = c(paste(font.tickslab),"black"),
        ylab=ylab,
        xlab=xlab,
        legend = legend_position,
        legend.title = legend.title)
}
