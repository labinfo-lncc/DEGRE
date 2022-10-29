#' Create a bar plot showing the number of downregulated and upregulated genes.
#'
#' @param results a data.frame object. It receives the output of the DEGRE function, filtered or not, as input.
#' @param log2FC_cutoff it stores the cutoff of the log2FoldChange. The default is 1.
#' @param downregulated_color the bar color related to the number of downregulated genes. The default is "coral2".
#' @param upregulated_color the bar color related to the number of upregulated genes. The default is "cornflowerblue".
#' @param xlab the x lab text. The default is "Regulation".
#' @param ylab the y lab text. The default is "Number of genes".
#' @param font.x the font size of the x axis. The default is 10.
#' @param font.y the font size of the y axis. The default is 10.
#' @param font.tickslab the font size of the ticks lab. The default is 10.
#' @param legend_position you need to specify here the position of the legend. The default is "right".
#' @param legend.title the title of the legend. The default is "Regulation".
#'
#' \value{No return value, called for side effects}
#'
#' @example
#' # Reading the count matrix and the design matrix for an example:
#' dir <- system.file("data", package = "DEGRE")
#' tab <- read.csv(file.path(dir,"count_matrix_for_example.csv"))
#' row.names(tab) <- tab[,1]; tab <- tab[,-1]
#' des <- read.csv(file.path(dir,"design_matrix_for_example.csv"))
#' # Running DEGRE function:
#' results <- DEGRE(count_matrix = tab,
#'              p_value_adjustment = "BH",
#'              design_matrix = des,
#'              formula = "condition + (1|sex)")
#' # Running the BarGraphDEGRE function
#' BarGraphDEGRE(results = results,
#'            log2FC_cutoff = 1,
#'            downregulated_color = "coral2",
#'            upregulated_color = "cornflowerblue",
#'            xlab = "Regulation",
#'            ylab = "Number of genes",
#'            font.x = 10,
#'            font.y = 10,
#'            font.tickslab = 10,
#'            legend_position = "right",
#'            legend.title = "Regulation")
#'
#' @export
BarGraphDEGRE <- function(results,
                          log2FC_cutoff = 1,
                          downregulated_color = "coral2",
                          upregulated_color = "cornflowerblue",
                          xlab = "Regulation",
                          ylab = "Number of genes",
                          font.x = 10,
                          font.y = 10,
                          font.tickslab = 10,
                          legend_position = "right",
                          legend.title = "Regulation"){

  if(missing(results))
    stop("You need to enter with the output of the DEGRE function.")

  log2FC_cutoff <- as.numeric(log2FC_cutoff)
  font.x <- as.numeric(font.x)
  font.y <- as.numeric(font.y)
  font.tickslab <- as.numeric(font.tickslab)

  colnames(results) <- c("ID","log2FC","P_value","Q_value")

  filter_upregulated <- results[results$log2FC >= log2FC_cutoff,]
  if(dim(filter_upregulated)[1] == 0){
    filter_upregulated <- data.frame(ID = NULL,
                                     log2FC = NULL,
                                     P_value = NULL,
                                     Q_value = NULL,
                                     reg = NULL)
  } else{
    filter_upregulated$reg <- "Upregulated"
  }

  filter_downregulated <- results[results$log2FC <= -log2FC_cutoff,]
  if(dim(filter_downregulated)[1] == 0){
    filter_downregulated <- data.frame(ID = NULL,
                                       log2FC = NULL,
                                       P_value = NULL,
                                       Q_value = NULL,
                                       reg = NULL)
  } else{
    filter_downregulated$reg <- "Downregulated"
  }

  ds_bargraph <- rbind(filter_upregulated, filter_downregulated)

  if(dim(ds_bargraph)[1] == 0){
    stop("The log2FC cutoff you enter does not bring any results.")
  }

  graph <- ggplot(data = ds_bargraph, aes(x = reg, fill = reg))+
    geom_bar()+
    scale_fill_manual(values = c(downregulated_color, upregulated_color))

  ggpar(graph,
        font.x = c(paste(font.x),"black"),
        font.y = c(paste(font.y),"black"),
        font.tickslab = c(paste(font.tickslab),"black"),
        ylab=ylab,
        xlab=xlab,
        legend = legend_position,
        legend.title = legend.title)
}
