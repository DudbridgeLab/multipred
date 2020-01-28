#' @param x Matrix of predicted risks.  Each row corresponds to an individual, each column to an outcome.  Each entry should be a risk between 0 and 1.
#' @param y Matrix of outcomes.  Each row corresponds to an individual, each column to an outcome.  Must contain binary outcomes coded as 0 and 1.
#' @param thresh Vector of risk thresholds.  For each row of x, each outcome is predicted to occur for which the risk exceeds the corresponding element of thresh.  These predictions are then compared to the elements of y.
#' If NULL, which is the default, concordance is the only measure that can be calculated.
#'
#' @return A list with the following components
#' @return \code{sens} Sensitivity
#' @return \code{spec} Specificity
#' @return \code{PPV} Positive predictive value
#' @return \code{NPV} Negative predictive value
#' @return \code{C} Concordance
#' @return \code{RU} Relative utility
