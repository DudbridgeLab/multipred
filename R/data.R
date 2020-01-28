#' Polygenic risk scores
#'
#' Simulated liabilities and polygenic risk scores
#' for six complex diseases, with corresponding disease outcomes
#' and risk predictions.
#'
#' @format A list with 7 elements.
#' Each element is a matrix with 10,000 rows and 6 columns
#' corresponding to simulated data on 10,000 individuals for 6 diseases.
#' Diseases are modelled on type-2 diabetes, coronary artery disease,
#' Crohn's disease, ulcerative colitis,
#' schizophrenia and rheumatoid arthritis.
#'
#' \describe{
#' \item{\code{PRS}}{Polygenic risk scores}
#' \item{\code{liability}}{Liabilities}
#' \item{\code{risk}}{Disease risks corresponding to \code{PRS}}
#' \item{\code{disease}}{Disease outcomes corresponding to \code{liability}}
#' \item{\code{VL}}{Variance-covariance matrix of \code{liability}}
#' \item{\code{VX}}{Variance-covariance matrix of \code{PRS}}
#' \item{\code{prevalence}}{Disease prevalences}
#' }
"PRSdata"
