# Inferring differentially expressed genes using generalized linear mixed models in the DEGRE R package

[![version](https://www.r-pkg.org/badges/version/DEGRE)](https://www.r-pkg.org/badges/version/DEGRE)
[![lastmonth](https://cranlogs.r-pkg.org/badges/last-month/DEGRE?color=blue)](https://r-pkg.org/pkg/DEGRE)
[![lastweek](https://cranlogs.r-pkg.org/badges/last-week/DEGRE)](https://cranlogs.r-pkg.org/badges/last-week/DEGRE)
[![lastday](https://cranlogs.r-pkg.org/badges/last-day/DEGRE)](https://cranlogs.r-pkg.org/badges/last-day/DEGRE)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/DEGRE)](https://cranlogs.r-pkg.org/badges/grand-total/DEGRE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7636564.svg)](https://doi.org/10.5281/zenodo.7636564)


## Scope

![alt text](https://github.com/labinfo-lncc/DEGRE/blob/main/man/figures/DEGRE%20logo2.png?raw=true)

The DEGRE is an R package that aims to identify Differentially Expressed Genes (DEGs) in a pairwise manner and considers the insertion of the individuals' random effects in the experimental design. This package has a preprocessing step responsible for filtering genes that could impair DEGs’ inference. For DEGs' inference, DEGRE uses Generalized Linear Mixed Model (GLMM) with the negative binomial distribution.



The repository for DEGRE is at GitHub on the https://github.com/labinfo-lncc/DEGRE. In this website you can report a bug and get help.



## Citation

Terra Machado D, Bernardes Brustolini OJ, Côrtes Martins Y, Grivet Mattoso Maia MA, Ribeiro de Vasconcelos AT. 2023. Inference of differentially expressed genes using generalized linear mixed models in a pairwise fashion. PeerJ 11:e15145 https://doi.org/10.7717/peerj.15145



## Installation of the DEGRE package

To install the DEGRE package, you must run:

```R
install.packages("DEGRE")
```



### Needed packages

The DEGRE package has several functions from the following packages:

- parglm, version 0.1.7
- glmmTMB, version 1.1.3
- doParallel, version 1.0.17
- foreach, version 1.5.2
- tibble, version 3.1.6
- dplyr, version 1.0.9
- ggplot2, version 3.3.5
- ggpubr, version  0.4.0
- ggrepel, version 0.9.1

Attention: if you are having trouble with the glmmTMB package, you can try removing and reinstalling it.



## Quick start

To use this package, you must have a matrix containing the read counts for each gene in each biological replicate and a design matrix containing the sample information. In the count matrix, the rows represent genes, while the columns represent the biological replicates of both experimental conditions. In the design matrix, the columns represent the relevant information of the samples.


To start using the package, you need to load the library:

```R
library("DEGRE")
```

Then you can apply the DEGRE's quick start function:

```R
# Reading the count matrix and the design matrix for an example:
dir <- system.file("data", package = "DEGRE")
tab <- read.csv(file.path(dir,"count_matrix_for_example.csv"))
row.names(tab) <- tab[,1]; tab <- tab[,-1]
des <- read.csv(file.path(dir,"design_matrix_for_example.csv"))

# Running DEGRE function:
results <- DEGRE(count_matrix = tab,
                 p_value_adjustment = "BH",
                 design_matrix = des,
                 formula = "condition + (1|sex)")
colnames(results) <- c("ID","log2FC","P.value","Q.value","avaragelogCPM")

```


## The output of the DEGRE function

Its output (represented as the `results` above) is a `data.frame` object. This object has the columns:

- ID - the gene IDs that the user inserted in the `count_matrix`. 
- log2FC - the log2 fold-change.
- *P*-value - the *P*-value computed for each gene in the Wald test.
- *Q*-value - the corrected *P*-values using Benjamini Hochberg (BH) or Bonferroni (BON) methods.
- averagelogCPM - the log2CPM for each gene.



## More about the DEGRE function

This function has the following arguments:

```R
results <- DEGRE(count_matrix, 
                 design_matrix,
                 formula,
                 p_value_adjustment)
```

As you can see, the DEGRE function has important parameters, which are described below:

- `count_matrix` - a `data.frame` object. It receives the matrix as input.
- `design_matrix` - a `data.frame` object. It receives the experimental design matrix. The sample names must be identified in the first column. This matrix can also have more columns with information for the fixed and the random effects for the samples.
- `formula` - it receives fixed and random effects descriptions.
- `p_value_adjustment` - All the *P*-values computed must be corrected and the DEGRE package offers two possibilities: "BH" (Benjamini-Hochberg) correction (default) and "BON" (Bonferroni) correction.

Obs.: The first column of `design_matrix` must be labeled "sample".

## What are random effects?
In RNA-seq, random effects are used to account for variation in gene expression levels that is due to sources other than the experimental factors of interest. This variation can come from a variety of sources, such as differences in the physiology of the individuals from whom the samples were collected, differences in the preparation of the samples, or differences in the sequencing process itself.

One example of random effects in RNA-seq is when you have multiple samples from different patients or individuals, as the expression level can differ among individuals due to the genetic variations, environment factors or other non-experimental factors. By including random effects that reflect the individual level, you can account for the variation among individuals, which can lead to more accurate and powerful inferences about the effects of the experimental factors.

Let's suppose we are doing an RNA-Seq experiment to identify DEGs in whole blood between a group of people that receives a specific treatment and a group that receive only water instead of medicine. The illustration of this comparison can be seen below:

<img src="https://github.com/labinfo-lncc/DEGRE/blob/main/man/figures/Github - fixed and random effects_FIGtreatcontrol.png" width="650">

Notice that we know that all people in group A received the treatment. In contrast, all people in group B received only water. This kind of comparison in which we can control the effect between both groups comprehends the fixed effects.

Suppose that in both groups, we have some heterogeneity related to the biological sex, in which there are men and women in different proportions. The illustration of this can be seen below:

<img src="https://github.com/labinfo-lncc/DEGRE/blob/main/man/figures/Github - fixed and random effects_FIGtreatcontrol_withbiologicalsex.png" width="650">

This kind of effect can be handle as a random effect to eliminate any potential bias related with variations between males and females¹²³. Other examples of random effects that can be harder to control may involve the different concentrations of biochemical molecules inside both groups implicating in the gene expression.

1. Thoral, Elisa, Quentin Queiros, Damien Roussel, Gilbert Dutto, Eric Gasset, David J. McKenzie, Caroline Romestaing, Jean-Marc Fromentin, Claire Saraux, and Loïc Teulier. 2021. “Changes in Foraging Mode Caused by a Decline in Prey Size Have Major Bioenergetic Consequences for a Small Pelagic Fish.” The Journal of Animal Ecology 90 (10): 2289–2301.

2. Kleyheeg, Erik, Bart A. Nolet, Sandra Otero-Ojea, and Merel B. Soons. 2018. “A Mechanistic Assessment of the Relationship between Gut Morphology and Endozoochorous Seed Dispersal by Waterfowl.” Ecology and Evolution 8 (22): 10857–67.

3. Donovan, Margaret K. R., Agnieszka D’Antonio-Chronowska, Matteo D’Antonio, and Kelly A. Frazer. 2020. “Cellular Deconvolution of GTEx Tissues Powers Discovery of Disease and Cell-Type Associated Regulatory Variants.” Nature Communications 11 (1): 955.

It is important to highlight that DEGRE is suitable for non-longitudinal experimental designs and future analyses including these types of experiments can be conducted. It is interesting to note that DEGRE can be applied in scenarios that have random effects that cause disturbances that are more difficult to be considered in programs that only accept fixed effects. Thus, it is also interesting to consider that increasing the number of biological replicates increases the detection power of true DEGs by DEGRE. 

## Some applications you can do!

### Filtering the gene IDs from the final data frame

Here are some codes you can apply to filter specific gene IDs.

Suppose you want to filter the results related to the ENSMUSG00000000881 gene ID. Then, you can try:

```R
filtering_gene <- results[results$ID == "ENSMUSG00000000881",]
print(filtering_gene)
```

Applying this filtering step, you get the following result:

```R
#                 ID   log2FC              P.value      Q.value
# ENSMUSG00000000881 -3.71189 1.83797593742305e-64 7.765448e-64
```



However, suppose you have specific gene IDs in a vector, for example:

```R
gene_IDs_to_filter <- c("ENSMUSG00000000881",
                        "ENSMUSG00000002820",
                        "ENSMUSG00000005610",
                        "ENSMUSG00000015484")
```

We can use the merge function to get the gene results. First, consider converting it into data.frame object:

```R
gene_IDs_to_filter <- data.frame(ID = gene_IDs_to_filter)
print(gene_IDs_to_filter)
```

```R
#                 ID
# ENSMUSG00000000881
# ENSMUSG00000002820
# ENSMUSG00000005610
# ENSMUSG00000015484
```

Then, you can use the merge function:

```R
filtering_gene <- merge(results, gene_IDs_to_filter, by = "ID")
print(filtering_gene)
```

```R
#                 ID     log2FC               P.value       Q.value
# ENSMUSG00000000881 -3.7118904  1.83797593742305e-64  7.765448e-64
# ENSMUSG00000002820 -4.4200312 5.10465155023083e-105 2.974780e-104
# ENSMUSG00000005610  1.4922973 2.41432633104465e-115 1.522467e-114
# ENSMUSG00000015484 -0.6455689  7.54077969856174e-06  1.098614e-05
```



### Getting the results from a *Q*-value cutoff

To filter the *Q*-value from a specific cutoff, you can do the following:

```R
results_q_value_cutoff <- results[results$`Q.value` < 0.05,]
```

Here we filtered based on the 5% of significance as a cutoff, but you can change it.

To check if it works properly, you can check the maximum *Q*-value from `results` and `results_q_value_cutoff`:

```R
max(results_q_value_cutoff$`Q.value`)
```

```R
# 0.0493556
```

Comparing it with the `results`, i.e., the data frame without the filtering, you get:

```R
max(results$`Q.value`)
```

```R
# 0.9519065
```

In which it shows you that the filter works properly.



### Getting the gene annotation from biomaRt package

If you only have the gene IDs but want the gene names, you can use the biomaRt package to get ensembl information.

First, if you don't have the biomaRt installed, please install it following the recommendations [here](https://bioconductor.org/packages/release/bioc/html/biomaRt.html).

Then, you must do the library:

```R
library("biomaRt")
```

Using the following code, you can get the gene names from the *Mus musculus* species:

```R
ensembl <- useMart("ensembl", 
                   dataset="mmusculus_gene_ensembl")

gene <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), 
              filters = 'ensembl_gene_id', 
              values = results$ID, 
              mart = ensembl)
```

You can use the merge function to bring the DEGRE output information of the genes. Attention: you need to rename the columns of the gene data frame in a way that the column with the gene IDs will have the same name in both data frames: `results` from the DEGRE function (DEGRE package) and `gene` from the getBM function (biomaRt package).

```R
colnames(gene)[1] <- "ID"
DEGREresults_gene_name <- merge(results, gene, by = "ID")
```

The data frame `DEGREresults_gene_name` now has the gene names from ensembl and all the information of the DEGRE function output.



## Plots that the DEGRE package offers to you

### Volcano plot

The DEGRE package has the VolcanoDEGRE function to visualize the proportion of downregulated and upregulated genes by applying a log2FC cutoff.

Here follows the quick start of the function with default parameters:

```R
VolcanoDEGRE(results = results,
             log2FC_cutoff = 1,
             padj = 0.05,
             delabel = "",
             font.x = 10,
             font.y = 10,
             font.tickslab = 10,
             downregulated_color = "coral2",
             upregulated_color = "cornflowerblue",
             xlab = "log2Foldchange",
             ylab = "-log10(P-value)",
             legend_position = "right",
             legend.title = "Regulation")
```

In the `results` argument, you must specify the data frame you want to plot. The following result is shown below.

<img src="https://github.com/labinfo-lncc/DEGRE/blob/main/man/figures/Volcano_DEGRE_example.png" width="500">

We recommend you to save it using the png function with the following parameters:

```R
png(filename = "Volcano_DEGRE_example.png", 
    width = 1700, 
    height = 1200, 
    res = 350)

# Then, you paste the code to do the bar plot.
Paste here the code to do the bar plot and run it.

# And so, you close the plotting:
dev.off()
```

Be free to change any of the arguments on the VolcanoDEGRE function.

The explained parameters are described below:

- `results` - a `data.frame` object. It receives the output of the DEGRE function, filtered or not, as input.

- `log2FC_cutoff` - it stores the cutoff of the log2FoldChange.

- `padj` - it stores the cutoff of the *P*-adjusted value (*Q*-value).

- `font.x` - the font size of the x axis.

- `font.y` - the font size of the y axis.

- `font.tickslab` - the font size of the ticks lab.

- `downregulated_color` - the colors of the downregulated genes. The default is "coral2".

- `upregulated_color` - the colors of the upregulated genes. The default is "cornflowerblue".

- `xlab` - the x lab text. The default is "log2Foldchange".

- `ylab` - the y lab text. The default is "-log10(P-value)".

- `legend_position` - you need to specify here the position of the legend. The default is "right".

- `legend.title` - the title of the legend.


To label the up and downregulated genes in `results` `data.frame`, you can run:

```R
log2FC_cutoff <- 1   # You can change this value.
padj <- 0.05   # You can change the 5% for the adjusted P-value cutoff.
results$diffexpressed <- "NO"
results$diffexpressed[results$log2FC > log2FC_cutoff & results$`Q.value` < padj] <- "UP"
results$diffexpressed[results$log2FC < -log2FC_cutoff & results$`Q.value` < padj] <- "DOWN"
```

### Bar plot

The DEGRE package has the BarGraphDEGRE function to get a bar plot showing the number of downregulated and upregulated genes.

Here follows the quick start of the function with default parameters:

```R
BarGraphDEGRE(results = results,
                     log2FC_cutoff = 1, 
                     downregulated_color = "coral2",
                     upregulated_color = "cornflowerblue",
                     xlab = "Regulation",
                     ylab = "Number of genes",
                     font.x = 10,
                     font.y = 10,
                     font.tickslab = 10,
                     legend_position = "right",
                     legend.title = "Regulation")
```

At the `results` argument, you must specify the data frame you want to plot. The following result is shown below.

<img src="https://github.com/labinfo-lncc/DEGRE/blob/main/man/figures/BarPlot_DEGRE_example.png" width="500">

We recommend you to save it using the png function with the following parameters:

```R
png(filename = "BarPlot_DEGRE_example.png", 
    width = 1500, 
    height = 1200, 
    res = 350)

# Then, you paste the code to do the bar plot.
Paste here the code to do the bar plot and run it.

# And so, you close the plotting:
dev.off()
```

Be free to change any of the arguments on the BarGraphDEGRE function.

The explained parameters are described below:

- `results` - a `data.frame` object. It receives the output of the DEGRE function, filtered or not, as input.
- `log2FC_cutoff` - it stores the cutoff of the log2FoldChange.
- `downregulated_color` - the bar color related to the number of downregulated genes. The default is "coral2".
- `upregulated_color` - the bar color related to the number of upregulated genes. The default is "cornflowerblue".
- `xlab` - the x lab text. The default is "Regulation".
- `ylab` - the y lab text. The default is "Number of genes".
- `font.x` - the font size of the x axis.
- `font.y` - the font size of the y axis.
- `font.tickslab` - the font size of the ticks lab.
- `legend_position` - you need to specify here the position of the legend. The default is "right".
- `legend.title` - the title of the legend.
