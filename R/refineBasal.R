#' @name refineBasal
#' @title Refine basal-like muscle invasive bladder cancer
#' @description Refine pre-identified basal-like muscle invasive bladder cancer into either basal-inflamed or basal-noninflamed subtype based on a random forest classifier.
#'
#' @param expr A numeric expression matrix with row features and sample columns.
#' @param isBasal A logic vector to indicate if the sample is basal-like or not. NULL by default and all samples will be assumed as basal-like.
#' @param res.path A string value to indicate the path for saving the prediction result.
#' @param res.name A string value to indicate the name of the output prediction table.
#'
#' @return A data.frame storing the random forest prediction results, including `samID` for sample name, `prob` for probability as basal-inflamed, and `basal` for final classifications using 0.5 as cutoff.
#' @export
#'
#' @import varSelRF
#' @examples
refineBasal <- function(expr       = NULL,
                        isBasal    = NULL,
                        res.path   = getwd(),
                        res.name   = "REFINED_MIBC_BASAL") {

  # load internal random forest classifier for refining basal-like MIBC
  #data("rfClassifier.RData")

  # initial check
  if(!all(is.element(rfClassifier$selected.vars, rownames(expr)))) {
    stop(paste0("expression profile must have the following features: \n ", paste(rfClassifier$selected.vars, collapse = "\n "),"\n\nmissing required features: ", paste(setdiff(rfClassifier$selected.vars, rownames(expr)), collapse = "\n")))
  }
  if(max(expr) < 25 | (max(expr) >= 25 & min(expr) < 0)) {
    message("--expression profile seems to have veen standardised (z-score or log transformation), no more action will be performed.")
    gset <- expr
  }
  if(max(expr) >= 25 & min(expr) >= 0){
    message("--log2 transformation done for expression data.")
    gset <- log2(expr + 1)
  }

  # prediction using random forest classifier
  if(is.null(isBasal)) { # all samples were assumed to be basal-like
    pred <- predict(rfClassifier$rf.model,
                    newdata = subset(t(gset), select = rf$selected.vars),
                    type = "prob")

    rf.res <- data.frame(samID = colnames(gset),
                         prob = as.numeric(pred[,1]),
                         basal = ifelse(as.numeric(pred[,1]) > 0.5,"inflamed","noninflamed"),
                         stringsAsFactors = F)

  } else {
    pred <- predict(rfClassifier$rf.model,
                    newdata = subset(t(gset[,isBasal]), select = rf$selected.vars), # use basal-like samples only
                    type = "prob")

    rf.res <- data.frame(samID = colnames(gset[,isBasal]),
                         prob = as.numeric(pred[,1]),
                         basal = ifelse(as.numeric(pred[,1]) > 0.5,"inflamed","noninflamed"),
                         stringsAsFactors = F)
  }

  return(rf.res)
}
