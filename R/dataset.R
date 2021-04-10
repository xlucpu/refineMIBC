#' Template for classifying MIBC using NTP
#'
#' A data frame storing the subtype-specific markers for MIBC.
#'
#' @format A data frame with 120 rows (30 markers for each subtype) and 3 variables.
#' \describe{contains ‘templates’ the subtype-specific up-regulated markers for MIBC.
"templates"

#' Random forest classifier for refining basal-like MIBC
#'
#' An R object derived from varSelRF.
#'
#' @format An object of class varSelRF
#' \describe{contains ‘rfClassifier’ an object of class varSelRF and can be used for external prediction.
"rfClassifier"

