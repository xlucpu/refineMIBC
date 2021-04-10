#' @name refineMIBC
#' @title Refine muscle-invasive bladder cancer
#' @description Refine muscle-invasive bladder cancer into four subtypes based on transcriptome profile by nearest template prediction algorithm.
#'
#' @param expr A numeric expression matrix with row features and sample columns.
#' @param scaleFlag A logic value to indicate if the expression data should be further scaled. TRUE by default.
#' @param centerFlag A logic value to indicate if the expression data should be further centered. TRUE by default.
#' @param distance A string value to indicate the distance measurement. Allowed values contain c('cosine', 'pearson', 'spearman', 'kendall'); 'cosine' by default.
#' @param seed An integer value for p-value reproducibility; 123456 by default.
#' @param doPlot A logic value to indicate whether to produce prediction heatmap; FALSE by default.
#' @param fig.path A string value to indicate the output path for storing the nearest template prediction heatmap.
#' @param fig.name A string value to indicate the name of the nearest template prediction heatmap.
#' @param res.path A string value to indicate the path for saving the prediction result.
#' @param res.name A string value to indicate the name of the output prediction table.
#'
#' @return A figure of predictive heatmap by NTP (.pdf) and a data.frame storing the results of nearest template prediction.
#' @export
#' @importFrom CMScaller ntp subHeatmap
#' @importFrom grDevices pdf dev.off
#' @references Hoshida, Y. (2010). Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE 5, e15543.
#'
#'             Eide, P.W., Bruun, J., Lothe, R.A. et al. CMScaller: an R package for consensus molecular subtyping of colorectal cancer pre-clinical models. Sci Rep 7, 16618 (2017).
#' @examples
refineMIBC <- function(expr       = NULL,
                       scaleFlag  = TRUE,
                       centerFlag = TRUE,
                       distance   = "cosine",
                       seed       = 123456,
                       doPlot     = FALSE,
                       fig.path   = getwd(),
                       fig.name   = "NTP_PREDICTED_HEATMAP",
                       res.path   = getwd(),
                       res.name   = "PREDICTED_MIBC_SUBTYPE") {

  # load internal templates for classifying MIBC
  data("templates.RData")

  # initial checking
  if(!is.element(distance, c("cosine", "pearson", "spearman", "kendall"))) {
    stop("the argument of distance should be one of cosine, pearson, spearman, or kendall.")
  }
  if(max(expr) < 25 | (max(expr) >= 25 & min(expr) < 0)) {
    message("--expression profile seems to have veen standardised (z-score or log transformation), no more action will be performed.")
    gset <- expr
  }
  if(max(expr) >= 25 & min(expr) >= 0){
    message("--log2 transformation done for expression data.")
    gset <- log2(expr + 1)
  }

  com_feat <- intersect(rownames(gset), templates$probe)
  message(paste0("--original template has ",nrow(templates), " biomarkers and ", length(com_feat)," are matched in external expression profile."))
  gset <- gset[com_feat, , drop = FALSE]
  templates <- templates[which(templates$probe %in% com_feat), , drop = FALSE]

  if(is.element(0,as.numeric(table(templates$class)))) {
    stop("at least one class has no probes/genes matched in template file!")
  }

  # scale data
  emat <- t(scale(t(gset), scale = scaleFlag, center = centerFlag))

  # perform NTP
  if(doPlot) {
    outFig <- paste0(fig.name,".pdf")
    pdf(file = file.path(fig.path, outFig), width = 5, height = 5)
    ntp.res <- ntp(emat      = emat,
                   templates = templates,
                   doPlot    = doPlot,
                   nPerm     = 1000,
                   distance  = distance,
                   nCores    = 1,
                   seed      = seed,
                   verbose   = TRUE)
    invisible(dev.off())

  } else {
    ntp.res <- ntp(emat      = emat,
                   templates = templates,
                   doPlot    = doPlot,
                   nPerm     = 1000,
                   distance  = distance,
                   nCores    = 1,
                   seed      = seed,
                   verbose   = TRUE)
  }

  return(ntp.res)
}
