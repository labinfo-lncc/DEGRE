
## Rodar a função 'DEGRE()' antes das etapas abaixo.


teste_count_matrix <- read.table("Área de Trabalho/count_matrix_for_example.csv",
                                 header = TRUE, sep = ",")
teste_design_matrix <- read.table("Área de Trabalho/design_matrix_for_example.csv",
                                  header = TRUE, sep = ",")

results <- DEGRE(GSE = teste_count_matrix,
                 num_reps = 2,
                 p_value_adjustment = "BH",
                 design_matrix = teste_design_matrix,
                 formula = "condition + (1|sex)")

## Cutoff of the significance level
results_significance <- results[results$`Q-value` < 0.05,]

library(ggplot2)
library(ggpubr)
## To run the graph below, you must have filtered with the cutoff of the significance level.
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

BarGraphDEGRE(results = results_significance)
