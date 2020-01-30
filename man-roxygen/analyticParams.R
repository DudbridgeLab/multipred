#' @param VL Variance-covariance matrix of liability.  Must have 1 on diagonal.
#' @param VX Variance-covariance matrix of predictors.
#' @param VLX Cross-covariance matrix between liabilities and predictors.
#' Entry on row i, column j, is covariance between liability i and predictor j.
#' Diagonal entries are the liability variances explained for each trait.
#' @param thresh Vector of risk thresholds for predicting an event.
#' If NULL, which is the default, concordance is the only measure that can be calculated.
#' @param prev Vector of prevalences, ie population risks, for each trait.
#'
#' @return A list with the following components
#' @return \code{sens} Sensitivity
#' @return \code{spec} Specificity
#' @return \code{PPV} Positive predictive value
#' @return \code{NPV} Negative predictive value
#' @return \code{C} Concordance
#' @return \code{RU} Relative utility
