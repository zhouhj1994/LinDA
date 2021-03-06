% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linda.R
\name{linda}
\alias{linda}
\title{Perform LinDA (Linear model for differential abundance analysis)}
\usage{
linda(
  otu.tab,
  meta,
  formula,
  imputation = FALSE,
  pseudo.cnt = 0.5,
  p.adj.method = "BH",
  alpha = 0.05,
  prev.cut = 0,
  lib.cut = 1,
  winsor.quan = NULL,
  n.cores = 1
)
}
\arguments{
\item{otu.tab}{data frame or matrix representing observed OTU table. Row: taxa; column: samples.
NAs are not expected in OTU tables so are not allowed in function \code{linda}.}

\item{meta}{data frame of covariates. The rows of \code{meta} correspond to the columns of \code{otu.tab}.
NAs are allowed. If there are NAs, the corresponding samples will be removed in the analysis.}

\item{formula}{character. For example: \code{formula = '~x1*x2+x3+(1|id)'}. At least one fixed effect is required.}

\item{imputation}{TRUE or FALSE. If TRUE, zeros in \code{otu.tab} will be imputed;
otherwise, we add \code{pseudo.cnt} to each value in \code{otu.tab}.}

\item{pseudo.cnt}{a positive real value. Default is 0.5.}

\item{p.adj.method}{character; p-value adjusting approach. See R function \code{p.adjust}. Default is 'BH'.}

\item{alpha}{a real value between 0 and 1; significance level of differential abundance. Default is 0.05.}

\item{prev.cut}{a real value between 0 and 1; taxa with prevalence (percentage of nonzeros)
less than prev.cut are excluded. Default is 0 (no taxa will be excluded).}

\item{lib.cut}{a non-negative real value; samples with less than \code{lib.cut} read counts are excluded.
Default is 1 (no samples will be excluded).}

\item{winsor.quan}{a real value between 0 and 1; winsorization cutoff for \code{otu.tab}, e.g., 0.97. If not given,
winsorization process will not be conducted.}

\item{n.cores}{a positive integer. If \code{n.cores > 1} and formula is in a form of mixed-effect model,
\code{n.cores} parallels will be conducted. Default is 1.}
}
\value{
A list with the elements
\item{variables}{A vector of variable names of all fixed effects in \code{formula}. For example: \code{formula = '~x1*x2+x3+(1|id)'}.
Suppose \code{x1} and \code{x2} are numerical, and \code{x3} is a categorical variable of three levels: a, b and c.
Then the elements of \code{variables} would be \code{('x1', 'x2', 'x3b', 'x3c', 'x1:x2')}.}
\item{bias}{numeric vector; each element corresponds to one variable in \code{variables};
the estimated bias of the regression coefficients due to compositional effect.}
\item{output}{a list of data frames with columns 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj', 'reject',
 'df'; \code{names(output)} is equal to \code{variables}; the rows of the data frame are taxa.
 Note: if there are taxa being excluded due to \code{prev.cut}, the number of rows of the output data frame
 is not equal to the number of rows of \code{otu.tab}. Taxa are identified by the rownames.
 If the rownames of \code{otu.tab} are NULL, then \code{1 : nrow(otu.tab)} is set as the rownames of \code{otu.tab}.
 \itemize{
 \item{baseMean:}{ 2 to the power of intercept coefficients (normalized by one million)}
 \item{log2FoldChange:}{ bias-corrected coefficients of the variable of interest}
 \item{lfcSE:}{ standard errors of the coefficients of the variable of interest}
 \item{stat:}{ \code{log2FoldChange / lfcSE}}
 \item{pvalue:}{ \code{2 * pt(-abs(stat), df)}}
 \item{padj:}{ \code{p.adjust(pvalue, method = p.adj.method)}}
 \item{reject:}{ \code{padj <= alpha}}
 \item{df:}{ degree of freedom. Number of samples minus number of variables for fixed-effect model; estimates
 from R package \code{lmerTest} with Satterthwaite's method of approximation for linear mixed-effect model.}
 }}
\item{otu.tab.use}{the OTU table used for analysis (the \code{otu.tab} after the preprocessing:
samples that have NAs in the variables in \code{formula} or have less than \code{lib.cut} read counts are removed;
taxa with prevalence less than \code{prev.cut} are removed and data is winsorized if \code{!is.null(winsor.quan)};
and zeros are treated, i.e., imputed or pseudo-count added).}
\item{meta.use}{the meta data used for analysis (only variables in \code{formula} are stored; samples that have NAs
or have less than \code{lib.cut} read counts are removed; numerical variables are scaled).}
}
\description{
The function implements a simple, robust and highly scalable approach to tackle
the compositional effects in differential abundance analysis. It fits linear regression models
on the centered log2-ratio transformed data, identifies a bias term due to the transformation
and compositional effect, and corrects the bias using the mode of the regression coefficients.
It supports the mixed-effect models.
}
\examples{

#install package "phyloseq" for importing "smokers" dataset
ind <- smokers$meta$AIRWAYSITE == 'Throat'
otu.tab <- as.data.frame(smokers$otu[, ind])
meta <- cbind.data.frame(Smoke = factor(smokers$meta$SMOKER[ind]),
                         Sex = factor(smokers$meta$SEX[ind]),
                         Site = factor(smokers$meta$SIDEOFBODY[ind]),
                         SubjectID = factor(smokers$meta$HOST_SUBJECT_ID[ind]))
ind1 <- which(meta$Site == 'Left')
res.left <- linda(otu.tab[, ind1], meta[ind1, ], formula = '~Smoke+Sex', alpha = 0.1,
                  prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
ind2 <- which(meta$Site == 'Right')
res.right <- linda(otu.tab[, ind2], meta[ind2, ], formula = '~Smoke+Sex', alpha = 0.1,
                  prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
rownames(res.left$output[[1]])[which(res.left$output[[1]]$reject)]
rownames(res.right$output[[1]])[which(res.right$output[[1]]$reject)]

linda.obj <- linda(otu.tab, meta, formula = '~Smoke+Sex+(1|SubjectID)', alpha = 0.1,
                   prev.cut = 0.1, lib.cut = 1000, winsor.quan = 0.97)
linda.plot(linda.obj, c('Smokey', 'Sexmale'), alpha = 0.1, lfc.cut = 1,
           legend = TRUE, directory = NULL, width = 11, height = 8)

}
\references{
Huijuan Zhou, Xianyang Zhang, Kejun He & Jun Chen. LinDA: Linear Models for Differential Abundance
Analysis of Microbiome Compositional Data.
}
\author{
Huijuan Zhou
}
