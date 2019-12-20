#' @param x Matrix of predicted risks.  Each row corresponds to an individual, each column to a trait.  Each entry should be a risk between 0 and 1.
#' @param y Matrix of traits.  Each row corresponds to an individual, each column to a trait.  Must contain binary events coded as 0 and 1.
#' @param thresh Vector of risk thresholds.  For each row of x, an event is predicted for each trait that exceeds the corresponding element of thresh.  These predictions are then compared to the elements of y.
#' If NULL, which is the default, concordance is the only measure that can be calculated.
#'
#' @return A list with the following components
#' @return \code{sens} Sensitivity
#' @return \code{spec} Specificity
#' @return \code{PPV} Positive predictive value
#' @return \code{NPV} Negative predictive value
#' @return \code{C} Concordance
#' @return \code{RU} Relative utility
