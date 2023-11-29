#' Reserve demand non-parametric curve
#'
#' Fits a non-parametric reserve demand curve between excess reserves and normalised rates
#'
#' @param type The type of the reserve demand curve. This can be any of \code{rforecast} for random forecast or \code{spline} for spline regression.
#' @param ... Additional arguments (unused).
#' @inheritParams curve
#'
#' @return Returns a model of class \code{npcurvir}. This includes
#' \itemize{
#'  \item \code{type} the type of the curve.
#'  \item \code{fit} the non-parametric model for the mean.
#'  \item \code{fitQ} The non-parametric model for the quantiles.
#'  \item \code{data} a list including the \code{y}, \code{x}, and \code{dummy} used for the fitting of the curve.
#'  \item \code{q} the interval used in the fitting of the curve.
#' }
#'
#' @seealso \code{\link{predict.npcurvir}}, and \code{\link{plot.npcurvir}}.
#'
#' @inherit curve author references
#'
#' @examples
#' # Use ECB example data
#' rate <- ecb$rate
#' x <- ecb$x[,1,drop=FALSE]
#' npcurve(x,rate)
#'
#' @export npcurve

npcurve <- function(x,y,type=c("rforest","spline"),dummy=NULL,q=NULL,...){

  type <- match.arg(type,c("rforest","spline"))

  # Prepare data as a data.frame
  dat <- datPrepNP(y,x,dummy)

  # Optimise curve
  switch(type,
         "rforest" = {fit <- randomForest(dat[,-1,drop=FALSE],y)},
         "spline" = {
           # Create the estimation equation
           if (ncol(x)*10 >= nrow(dat)){
             kmax <- floor(nrow(dat)*.7/(ncol(x)+as.numeric(!is.null(dummy))))
             eq <- as.formula(paste0("y~",paste0(paste0("s(",colnames(x)),paste0(",k=",kmax,")"),collapse="+")))
           } else {
             eq <- as.formula(paste0("y~",paste0(paste0("s(",colnames(x)),")",collapse="+")))
           }
           # Fit the model
           fit <- gam(eq,data=dat)
           })

  # Estimate quantiles if needed
  if (!is.null(q)){
    switch(type,
           "rforest" = {fitQ <- quantregForest(dat[,-1,drop=FALSE],y,keep.inbag=TRUE)},
           "spline" = {
             p <- (1-q)/2
             fitQ <- mqgam(formula(fit),dat,qu=c(p,1-p))
           })
  } else {
    fitQ <- NULL
  }

  # Save curve model
  data = list(y=y,x=x,dummy=dummy)
  model <- list(type=type,fit=fit,fitQ=fitQ,data=data,q=q)
  model <- structure(model,class="npcurvir")

  return(model)

}

#' Predict method for npcurvir reserve demand models
#'
#' Predicted values based on npcurvir model object
#'
#' @param object A model fit with \code{\link{npcurve}}.
#' @param newdata New input data organised as the x matrix in \code{\link{curve}}. If \code{NULL} then the data used to fit the model is re-used.
#' @param newdummy New input dummy organised as the dummy vector in \code{\link{curve}}. If \code{NULL} then the dummy used to fit the model is re-used.
#' @param ... Further arguments (unused)
#'
#' @return Returns a matrix of predicted values. If the model has estimates for intervals then it will provide upper and lower intervals.
#'
#' @inherit curve references author
#' @seealso \code{\link{npcurve}}.
#'
#' @examples
#' # Use ECB example data
#' rate <- ecb$rate
#' x <- ecb$x[,1,drop=FALSE]
#' fit <- npcurve(x,rate)
#' predict(fit)
#'
#' @describeIn predict Predicted values for non-parametric curves
#' @export
#' @method predict npcurvir
predict.npcurvir <- function(object, newdata=NULL, newdummy=NULL,...){

  # Use old data if nothing new is provided
  if (is.null(newdata)){
    newdata <- object$data$x
  }

  if (is.null(newdummy)){
    newdummy <- object$data$dummy
  }

  if (!is.null(newdummy)){
    if (length(newdummy) != nrow(newdata)){
      stop("Length of newdummy must match the rows of newdata.")
    }
  }

  # Prepare data
  newX <- datPrepNP(y=NULL,newdata,newdummy)

  switch(object$type,
         "rforest" = {
            yhatMean <- matrix(predict(object$fit,newdata=newX),ncol=1)
          },
         "spline" = {
           yhatMean <- matrix(predict(object$fit,newdata=newX,se.fit=TRUE)$fit,ncol=1)
         })
  colnames(yhatMean) <- "mean"

  # Predict intervals
  if (!is.null(object$fitQ)){

    p <- (1-object$q)/2
    switch(object$type,
           "rforest" = {
             yhatQ <- predict(object$fitQ,newdata=newX,what=c(p,1-p))
           },
           "spline" = {
             yhatQ <- qdo(object$fitQ,qu=c(p,1-p),predict,newdata=newX)
             yhatQ <- abind(yhatQ,along=2)
           })
    colnames(yhatQ) <- c(paste0("lower ",object$q*100,"%"),
                         paste0("upper ",object$q*100,"%"))
  } else {
    yhatQ <- NULL
  }

  return(cbind(yhatMean,yhatQ))

}

#' Cross-validated errors for non-parametric curve
#'
#' Obtain cross-validated errors for a non-parametric curve with given sample splits.
#'
#' @inheritParams npcurve
#' @param cvIndex A matrix detailing how the sample was split for the cross-validation. Output from \code{\link{cvfit}}.
#'
#' @return Returns summary cross-validated errors, comparable with the output from  \code{\link{cvfit}}.
#'
#' @inherit curve references author
#' @seealso \code{\link{cvfit}}, and \code{\link{cvfitplot}}.
#'
#' @examples
#' \donttest{
#'   # Use ECB example data
#'   rate <- ecb$rate
#'   x <- ecb$x[,1:3,drop=FALSE]
#'   cvKeep <- cvfit(x,rate,folds=5,alltype=c("logistic","arctan"),parallel=TRUE)
#'   # Get non-parametric curve cross-validated errors
#'   cvRF <- cvnpcurve(x,rate,cvKeep$cvIndex)
#'   cvSP <- cvnpcurve(x,rate,cvKeep$cvIndex,type="spline")
#' }
#'
#' @export cvnpcurve

# Get cross-validated errors from npcurve
cvnpcurve <- function(x,y,cvIndex,type="rforest",dummy=NULL){

  type <- match.arg(type,c("rforest","spline"))

  folds <- max(cvIndex[,1])
  cvError <- array(NA,c(1,folds))
  for (i in 1:folds){
    # Subset
    yTrn <- y[cvIndex[cvIndex[,1]!=i,2]]
    xTrn <- x[cvIndex[cvIndex[,1]!=i,2],,drop=FALSE]
    yTst <- y[cvIndex[cvIndex[,1]==i,2]]
    xTst <- x[cvIndex[cvIndex[,1]==i,2],,drop=FALSE]
    if (!is.null(dummy)){
      dTrn <- dummy[cvIndex[cvIndex[,1]!=i,2]]
      dTst <- dummy[cvIndex[cvIndex[,1]==i,2]]
    } else {
      dTrn <- dTst <- NULL
    }

    # Train for this fold
    fit <- npcurve(xTrn,yTrn,type=type,dummy=dTrn)
    # Predict for this fold
    yhat <- predict(fit,newdata=xTst,newdummy=dTst)
    # Obtain error
    cvError[,i] <- mean((yTst - yhat)^2)
  }

  errors <- c(mean(cvError),median(cvError),sd(cvError))
  names(errors) <- c("Mean MSE","Median MSE","Sd MSE")
  return(errors)
}

## Internal functions

#' @export
#' @method print npcurvir
print.npcurvir <- function(x, ...){
  type <- x$type

  writeLines(paste0(toupper(substr(x$type,1,1)), tolower(substr(x$type,2,nchar(x$type))), " type reserve demand non-parametric curve"))
  if (ncol(x$data$x)-1>0){
    writeLines(paste0(ncol(x$data$x)-1, " additional regressors included."))
  }
  if (!is.null(x$dummy)){
    writeLines("Shift indicator included.")
  }
  if (!is.null(x$q)){
    writeLines("")
    writeLines(paste0("Model has estimates for ",round(x$q*100,0),"% intervals."))
  }

}

# Data preparation for non-parametric curve model fit
datPrepNP <- function(y,x,dummy){

  # Check column names of x
  clnmX <- colnames(x)
  if (is.null(clnmX)){
    colnames(x) <- paste0("x",1:ncol(x))
  }

  if (!is.null(y)){
    y <- cbind(y)
    colnames(y) <- "y"
  }

  if (!is.null(dummy)){
    dummy <- cbind(dummy)
    colnames(dummy) <- "dummy"
  }

  dat <- as.data.frame(cbind(y,x,dummy))

  return(dat)

}
