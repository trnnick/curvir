#' ECB excess reserve data
#'
#' The short-term interest rate used for the European Central Bank (ECB) is the volume-weighted Euro Overnight Index Average (EONIA) rate.
#' Various explanatory variables are provided, as listed in the Chen et al. (2023).
#' Data are collected from 1999 to 2019, resulting in 239 maintenance periods where all data is available.
#'
#' @format \code{ecb} a list containing
#' \itemize{
#'   \item \code{rate} the normalised EONIA rate.
#'   \item \code{x} a matrix including the ECB excess reserves and various regressors.
#' }
#' @inherit curve references
#'
#' @examples
#' plot(ecb$x[,1],ecb$rate,xlab="Excess Reserves",ylab="Rate")
#'
"ecb"
