#' multipred
#'
#' Package for calculating measures of accuracy for risk predictors of multiple outcomes.
#'
#' Accuracy can be evaluated in four senses: outcome-wise, joint, and family-wise (weak sense and strong sense).
#' For convenience the weak family-wise sense is also called "screening", and the strong family-wise sense
#' simply "family-wise".  In each sense, accuracy can be measured empirically, within data sets given as input, or theoretically,
#' given parameters of an underlying multivariate liability threshold model.
#'
#' Throughout the documentation, an "outcome" means one of several binary variables observed in an individual,
#' and that an outcome "occurs" when the variable has the positive state.
#'
#' Outcome-wise measures calculate standard univariate measures of accuracy over all outcomes and individuals.
#'
#' Joint measures consider the prediction of all outcomes occuring simultaneously within an individual.
#'
#' Screening measures consider the prediction of at least one outcome occuring within an individual.
#' It is not necessary that the predicted outcomes are the same ones that actually occur.
#'
#' Family-wise measures consider the prediction of at least one outcome occuring within an individual.
#' There must be at least one predicted outcome that actually occurs.
#'
#' @section Functions:
#' \code{\link{outcomeWise}}
#'
#' \code{\link{joint}}
#'
#' \code{\link{screening}}
#'
#' \code{\link{familyWise}}
#'
#' \code{\link{analyticOutcomeWise}}
#'
#' \code{\link{analyticJoint}}
#'
#' \code{\link{analyticScreening}}
#'
#' \code{\link{analyticFamilyWise}}
#'
#' @docType package
#' @name multipred
NULL

