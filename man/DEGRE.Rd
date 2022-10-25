\name{DEGRE}
\alias{DEGRE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifying differentially expressed genes with DEGRE
}
\description{
Genes that are differentially expressed (DEGs) between two or more experimental conditions can be found in investigations using RNA-Seq data. The gene reads matrix is a count matrix that is used to infer DEGs, and fixed effects are often utilized in this inference. A high biological variability may impact the discovery of DEGs once it may be divergent between the fixed effects. This variability can be covered by the random effects. The DEGRE R package was designed to identify DEGs in RNA-Seq experiments considering fixed and the random effects on individuals. These effects are identified earlier in the experimental design matrix. DEGRE has the implementation of preprocessing procedures to clean the count matrix of the gene reads, and it fits a regression for each gene using the Generalized Linear Mixed Model with the negative binomial distribution,followed by a Wald test to test to assess the regression coefficients.
}
\usage{
DEGRE(count_matrix, p_value_adjustment = "BH", formula, design_matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count_matrix}{
The count matrix is a data frame that have the gene reads by samples. The rows have the IDs of the genes and the columns have the IDs of the biological replicates of the samples.
}
  \item{p_value_adjustment}{
The p-value adjustment involves the type of p-value correction the user wants to apply in the p-values multiple testing. The default on here is the Benjamini-Hochberg "BH", but can be alteres to Bonferroni "BON".
}
  \item{formula}{
The formula must be the fixed and the random factors that the user wants to apply in the inference of DEGs. All the levels from the factors must be present in the design and the count matrix.
}
  \item{design_matrix}{
The count matrix is a data frame in which the first column has the IDs of the biological replicates and the other columns have the fixed and random factors and the levels for each.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Manuscript is being written.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----

##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
