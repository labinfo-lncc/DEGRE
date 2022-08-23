# Inferring Differential Expression Genes with Random Effects using the DEGRE R package

## Scope

![alt text](https://github.com/labinfo-lncc/DEGRE/blob/main/DEGRE%20logo2.png?raw=true)

The DEGRE is an R package that aims the identification of Differential Expression Genes (DEGs) among two groups considering the insertion of the random effects in the experimental design matrix. These effects are identified previously in the experimental design. The package has the implementation of pre-processing steps to clean the count matrix of the gene reads, and it uses the Generalized Linear Mixed Model (GLMM) with the negative binomial distribution to infer the DEGs.



The repository for DEGRE is at GitHub on the www.websitehere.com. In this website you can report a bug and get help.



## Citation

The DEGRE package is under publication. If you use it in your research, please cite it when available.

You can cite it as here below:

```R
citation("DEGRE")
```

```R
# To cite package ‘DEGRE’ in publications use:
#   Douglas Terra Machado (2022). DEGRE: Identification of Differentially Expressed Genes
#   with Random Effects. R package version 1.0.
# A BibTeX entry for LaTeX users is
#   @Manual{,
#     title = {DEGRE: Identification of Differentially Expressed Genes with Random
# Effects},
#     author = {Douglas Terra Machado},
#     year = {2022},
#     note = {R package version 1.0},
#   }
#  ATTENTION: This citation information has been auto-generated from the package
#  DESCRIPTION file and may need manual editing, see ‘help("citation")’.
```



## Instalation of the DEGRE package

To install the DEGRE package, you must run:

```R
library(devtools)
install_github("labinfo-lncc/DEGRE",
                ref="main", auth_token = "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")
```



### Needed packages

The DEGRE package has functions that uses the following packages:

- parglm, version 0.1.7
- glmmTMB, version 1.1.3
- doParallel, version 1.0.17
- foreach, version 1.5.2
- car, version 3.0.12
- tibble, version 3.1.6
- gridExtra, version 2.3
- dplyr, version 1.0.9
- ggplot2, version 3.3.5
- ggpubr, version  0.4.0
- ggrepel, version 0.9.1

Attention: if you are having troubles with the glmmTMB package, you can try to remove it and re-install, it may work properly.



## Quick start

The DEGs are inferred according to its read counts among the replicates between two groups and the quantitative differences of the read counts are inferred to contribute in differential expression.

To use this package you must have a count matrix that contains the read counts for each gene in each of the sample replicates and also a design matrix that contains sample informations. In the count matrix, the genes are represented by the rows and the replicates are represented by the columns. In the design matrix, the columns represent the relevant information of the samples. 

To start using the package, you need to load the library:

```R
library("DEGRE")
```

Then you can apply the DEGRE function:

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

```



## Output of the DEGRE function

The output of the DEGRE function (represented as `results` above) is a `data.frame` object. This object has the columns:

- ID - the gene IDs that the user entered in the `count_matrix`. 
- log2FC - the log2 fold-change.
- *P*-value - the *P*-value computed for each gene in the Wald test (confirmar o nome certinho).
- *Q*-value - the corrected *P*-values using BH or BON methods.



## More about the DEGRE function

The DEGRE function has the following arguments:

```R
results <- DEGRE(count_matrix, 
                 num_reps,
                 design_matrix,
                 formula,
                 p_value_adjustment)
```

As you can see, the DEGRE function has important parameter in which are described below:

- `count_matrix` - a `data.frame` object. It receives the matrix with the reads for each gene as input. The rows are the gene IDs and the columns are the two groups with the identification of the replicates.
- `num_reps` - it is the argument that stores the integer number of replicates for each group and it must be the same.
- `design_matrix` - a `data.frame` object. It receives the experimental design matrix. The sample names must be identified in the first column. This matrix can also have more columns with informations for the fixed and the random effects for the samples.
- `formula` - it receives the fixed and random effects as a quotation mark argument.
- `p_value_adjustment` - All the *P*-values computed must be corrected and the DEGRE package offers two possibilities for this *P*-value correction: "BH" (Benjamini-Hochberg) correction, in which it is the default for the package, and "BON" (Bonferroni) correction.


## What are the fixed and random effects?
The fixed effects can be understood as variables which are easier to control in the experimental design. Let's suppose we are doing an RNA-Seq experiment to identify differentially expressed genes in whole blood between a group of people that receives a specific treatment compared to a group of people that receive only water instead of the medicine. The illustration of this comparison can be seen below:

<img src="https://github.com/labinfo-lncc/DEGRE/blob/main/Github - fixed and random effects_FIGtreatcontrol.png" width="650">

Notice that we know that all people in group A received the treatment while all people of group B received only water, this kind of comparison in which we can control the effect between both groups comprehends the fixed effects.

Now suppose that in the both groups we have some heterogeneity related to the biological sex, in which there are men and women in different proportions in the both groups. The illustration of this can be seen below:

<img src="https://github.com/labinfo-lncc/DEGRE/blob/main/Github - fixed and random effects_FIGtreatcontrol_withbiologicalsex.png" width="650">

This kind of effect in which we did not control is related to the random effect. Other examples of random effects that can be harder to control may involve the different concentration of biochemical molecules inside of the both groups and the concentration of these molecules involved in the gene expression.


## Some applications you can do!

### Filtering the gene IDs from the final data frame

Here follows some codes you can apply to filter specific gene IDs from the final data frame.

Suppose you want to filter the results related to the ENSMUSG00000000881 gene ID. Then, you can try:

```R
filtering_gene <- results[results$ID == "ENSMUSG00000000881",]
print(filtering_gene)
```

Applying this filtering step you will get the following result:

```R
#                 ID   log2FC              P-value      Q-value
# ENSMUSG00000000881 -3.71189 1.83797593742305e-64 7.765448e-64
```



However, suppose you have specific gene IDs in a vector, for example:

```R
gene_IDs_to_filter <- c("ENSMUSG00000000881",
                        "ENSMUSG00000002820",
                        "ENSMUSG00000005610",
                        "ENSMUSG00000015484")
```

We can use the merge function to get the genes results. First consider to convert it into a `data.frame` object:

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
#                 ID     log2FC               P-value       Q-value
# ENSMUSG00000000881 -3.7118904  1.83797593742305e-64  7.765448e-64
# ENSMUSG00000002820 -4.4200312 5.10465155023083e-105 2.974780e-104
# ENSMUSG00000005610  1.4922973 2.41432633104465e-115 1.522467e-114
# ENSMUSG00000015484 -0.6455689  7.54077969856174e-06  1.098614e-05
```



### Getting the results from a *Q*-value cutoff

To filter the *Q*-value from a specific cutoff, you can try the following:

```R
results_q_value_cutoff <- results[results$`Q-value` < 0.05,]
```

Here we filtered based on the 5% of significance as a cutoff, but you can change it in your way. 

To check if it works properly, you can check both the dimensions of `results` and `results_q_value_cutoff` data frames, but you can also try:

```R
max(results_q_value_cutoff$`Q-value`)
```

```R
# 0.0493556
```

Comparing it with the `results`, i.e., the data frame without the filtering, you get:

```R
max(results$`Q-value`)
```

```R
# 0.9519065
```

In which it shows to you that the filter works properly.



### Getting the gene annotation from biomaRt package

If you only have the gene IDs but want to get the gene names, you can do it using the biomaRt package to get ensembl informations.

First, if you don't have the biomaRt installed, please install it following the recommendations [here](https://bioconductor.org/packages/release/bioc/html/biomaRt.html).

Then, you must do the library:

```R
library("biomaRt")
```

Using the following code you can get the gene names from the *Mus musculus* specie:

```R
ensembl <- useMart("ensembl", 
                   dataset="mmusculus_gene_ensembl")

gene <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), 
              filters = 'ensembl_gene_id', 
              values = results$ID, 
              mart = ensembl)
```

To bring the DEGRE output informations, you can use the merge function. Attention: you need to rename the columns of the `gene` data frame in a way that the columns with the gene IDs will have the same name in both data frames: `results`, from the DEGRE function (DEGRE package), and `gene`, from the getBM function (biomaRt package). 

```R
colnames(gene)[1] <- "ID"
DEGREresults_gene_name <- merge(results, gene, by = "ID")
```

The data frame `DEGREresults_gene_name` now have the gene names from ensembl and all the informations of the DEGRE function output.



## Plots that the DEGRE package offers to you

### Volcano plot

The DEGRE package has the VolcanoDEGRE function to get a visualization of the proportion of downregulated and upregulated genes applying a log2FC cutoff.

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

At the `results` argument, you have to specify the data frame you want to plot. The following result is shown below.

<img src="https://github.com/labinfo-lncc/DEGRE/blob/main/Volcano_DEGRE_example.png" width="500">

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

The explained parameters is described below:

- `results` - a `data.frame` object. It receives the output of the DEGRE function, filtered by you or not, as input.

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

At the `results` argument, you have to specify the data frame you want to plot. The following result is shown below.

<img src="https://github.com/labinfo-lncc/DEGRE/blob/main/BarPlot_DEGRE_example.png" width="500">

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

The explained parameters is described below:

- `results` - a `data.frame` object. It receives the output of the DEGRE function, filtered by you or not, as input.
- `log2FC_cutoff` - it stores the cutoff of the log2FoldChange.
- `downregulated_color` - the color of the bar related to the number of downregulated genes. The default is "coral2".
- `upregulated_color` - the color of the bar related to the number of upregulated genes. The default is "cornflowerblue".
- `xlab` - the x lab text. The default is "Regulation".
- `ylab` - the y lab text. The default is "Number of genes".
- `font.x` - the font size of the x axis..
- `font.y` - the font size of the y axis..
- `font.tickslab` - the font size of the ticks lab.
- `legend_position` - you need to specify here the position of the legend. The default is "right".
- `legend.title` - the title of the legend.
