#' @title Healthcare data set
#'
#' @description Data set containing the healthcare expense of 129,257 customers of a Brazilian healthcare
#' company between 2006 and 2009.
#' @details The expenses are in Reais (Brazilian currency) and were deflated to the January 2006 value. In order to
#' fit the models of the paper in the references it is necessary to turn every observation whose expense is less than
#' R$ 100 into zero.
#' @return \item{ID}{The ID of the customer.}
#' @return \item{sex}{The sex of the customer.}
#' @return \item{age}{The age of the customer on the considered year.}
#' @return \item{expense}{The healthcare expense of the customer on the considered year.}
#' @return \item{log_expense}{The logarithm of the healthcare expense of the customer on the considered year.}
#' @return \item{year}{The considered year.}
#' @return \item{previous_expense}{The healthcare expense of the customer on the previous year.}
#' @return \item{log_previous_expense}{The logarithm of the healthcare expense of the customer on the previous year.}
#' @references Marcondes, D.; Peixoto, C.; Maia, A. C.; A Survey of a Hurdle Model for Heavy-Tailed Data Based on the Generalized Lambda Distribution. (2017) \emph{arxiv1712.02183}
"healthcare"
