% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refineMIBC.R
\name{refineMIBC}
\alias{refineMIBC}
\title{Refine muscle-invasive bladder cancer}
\usage{
refineMIBC(
  expr = NULL,
  scaleFlag = TRUE,
  centerFlag = TRUE,
  distance = "cosine",
  seed = 123456,
  doPlot = FALSE,
  fig.path = getwd(),
  fig.name = "NTP_PREDICTED_HEATMAP",
  res.path = getwd(),
  res.name = "PREDICTED_MIBC_SUBTYPE"
)
}
\arguments{
\item{expr}{A numeric expression matrix with row features and sample columns.}

\item{scaleFlag}{A logic value to indicate if the expression data should be further scaled. TRUE by default.}

\item{centerFlag}{A logic value to indicate if the expression data should be further centered. TRUE by default.}

\item{distance}{A string value to indicate the distance measurement. Allowed values contain c('cosine', 'pearson', 'spearman', 'kendall'); 'cosine' by default.}

\item{seed}{An integer value for p-value reproducibility; 123456 by default.}

\item{doPlot}{A logic value to indicate whether to produce prediction heatmap; FALSE by default.}

\item{fig.path}{A string value to indicate the output path for storing the nearest template prediction heatmap.}

\item{fig.name}{A string value to indicate the name of the nearest template prediction heatmap.}

\item{res.path}{A string value to indicate the path for saving the prediction result.}

\item{res.name}{A string value to indicate the name of the output prediction table.}
}
\value{
A figure of predictive heatmap by NTP (.pdf) and a data.frame storing the results of nearest template prediction.
}
\description{
Refine muscle-invasive bladder cancer into four subtypes based on transcriptome profile by nearest template prediction algorithm.
}
\examples{
library(refineMIBC)
load(system.file("extdata", "demo.RData", package = "refineMIBC", mustWork = TRUE)) # load example data
expr <- demo$MIBC.expr
subt <- demo$MIBC.subt
iCS  <- refineMIBC(expr       = expr,
                   scaleFlag  = TRUE,
                   centerFlag = TRUE,
                   doPlot     = FALSE)
head(iCS)
}
\references{
Hoshida, Y. (2010). Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE 5, e15543.

            Eide, P.W., Bruun, J., Lothe, R.A. et al. CMScaller: an R package for consensus molecular subtyping of colorectal cancer pre-clinical models. Sci Rep 7, 16618 (2017).
}
