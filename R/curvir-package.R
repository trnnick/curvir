#' @keywords internal package
#' @aliases curvir-package
"_PACKAGE"

#' curvir:Reserve demand curve modelling toolbox
#'
#' The \pkg{curvir} package provides tools for building reserve demand curves for central bank operations.
#'
#' @section Parametric curve modelling:
#' \itemize{
#'  \item \code{\link{curve}} - fit a parametric curve
#'  \item \code{\link{curveopt}} - optimise curve parameters (called through \code{\link{curve}})
#'  \item \code{\link{curvepred}} - predict curve values (called through \code{\link{predict}})
#'  \item \code{\link{invcurve}} - produce predictions from the inverse curve
#' }
#'
#' @section Non-parametric curve modelling:
#' \itemize{
#'  \item \code{\link{npcurve}} - fit a non-parametric curve
#' }
#'
#' @section Curve specification:
#' \itemize{
#'  \item \code{\link{cvfit}} - Joint automatic specification of curve type and regressors
#'  \item \code{\link{varselect}} - Automatic specification of regressors
#'  \item \code{\link{cvnpcurve}} - Provide cross-validated errors for non-paramteric curves
#'  \item \code{\link{cvfitplot}} - Summary plot for \code{\link{cvfit}} and \code{\link{cvnpcurve}}
#' }
#'
#' @section General:
#' \itemize{
#'  \item \code{\link{predict}}  Predict parametric and non-parametric curves
#'  \item \code{\link{plot}}  Plot parametric and non-parametric curves
#' }
#'
#' @section Data:
#' \itemize{
#'  \item \code{\link{ecb}} ECB dataset from Chen et al. (2023)
#' }
#'
#' @inherit curve references
#' @docType package
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}, Zhuohui Chen, \email{zchen4@imf.org}, Romain R. Veyrune.
#'
#' @name curvir
#'
#' @importFrom pso psoptim
#' @importFrom cvTools cvFolds
#' @importFrom abind abind
#' @importFrom parallel detectCores clusterApplyLB makeCluster stopCluster clusterCall
#' @importFrom pbapply pblapply pboptions
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics axTicks axis lines points title legend par
#' @importFrom stats IQR cor median optim prcomp quantile sd as.formula formula predict qnorm
#' @importFrom utils tail
#' @importFrom randomForest randomForest
#' @importFrom quantregForest quantregForest
#' @importFrom mgcv gam
#' @importFrom qgam mqgam qdo
#'
NULL
