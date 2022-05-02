# Inferring Differential Expression Genes with Random Effects using the DEGRE R package

[TOC]

## Scope

![](/home/douglast/Área de Trabalho/DEGRE logo2.png)

The DEGRE is an R package to identify Differential Expression Genes (DEGs) among two groups considering the insertion of the random effects in the experimental design matrix. These effects are identified previously in the experimental design matrix. The package has the implementation of pre-processing steps to clean the count matrix of the gene reads, and it uses the Generalized Linear Mixed Model (GLMM) with the negative binomial distribution to infer the DEGs.



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

To install the DEGRE package, you must use:

xxxxxxxxxxxxxx



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

  

## Quick start

The DEGs are inferred according to its read counts among the replicates between two groups and the quantitative differences of the read counts are inferred to contribute in its differential expression (Figure 1).

To use this package you must have a count matrix that contains the read counts for each gene in each of the biological replicates and also a design matrix that contains the information related to the samples.  In the count matrix, the genes are represented by the rows and the replicates are represented by the columns. In the design matrix, the columns represent the relevant information of the samples. 

To start using the package, you need to load the library:

```R
library("DEGRE")
```

Then you can apply the DEGRE function:
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
- `formula` - it receives the fixed and random effects as a XXXX (ver GLM como é falado) argument. All the text in here must be into quotation marks, for example: "your formula here".  
- `p_value_adjustment` - All the *P*-values computed must be corrected and the DEGRE package offers two possibilities for this *P*-value correction: "BH" (Benjamini-Hochberg) correction, in which it is the default for the package, and "BON" (Bonferroni) correction.



## Output of the DEGRE function

The output of the DEGRE function (represented as `results` above) is a `data.frame` object. This object has the columns:

- ID - the gene IDs that the user entered in the `count_matrix`. 
- log2FC - the log2 fold-change.
- *P*-value - the *P*-value computed for each gene in the Wald test (confirmar o nome certinho).
- *Q*-value - the corrected *P*-values using BH or BON methods.

xxx