#' Automatically select explanatory variables for a curve type
#'
#' Using cross-validation automatically select explanatory variables for a reserve
#' demand curve type.
#'
#' @inheritParams curveopt
#' @param folds Folds to use for cross-validation.
#' @param search Search strategy for variable inclusion. With \code{forward} the model starts from no explanatory variables, and evaluates potential additions. With \code{backward} the model starts with all variables and evaluates for potential exclusions.
#'
#' @return Returns a list with the recommended variable selection choice:
#' \itemize{
#'  \item \code{keep} a logical vector with which variables to keep.
#'  \item \code{errors} statistics of the cross-validated MSE error.
#' }
#'
#' @inherit curve author references
#' @seealso \code{\link{curve}}, and \code{\link{cvfit}}.
#'
#' @examples
#' \donttest{
#'   # Use ECB example data
#'   rate <- ecb$rate
#'   x <- ecb$x[,1:3,drop=FALSE]
#'   varKeep <- varselect(x,rate,folds=5)
#'   # Print result
#'   print(varKeep)
#'   # Fit curve with the selected variables
#'   curve(x[,varKeep$keep,drop=FALSE],rate)
#' }
#'
#' @export varselect

# Find the best set of x's using backwards search and cross-validation
varselect <- function(x,y,type="logistic",folds=10,constant=c(TRUE,FALSE),sign=NULL,reps=3,
                      search=c("backward","forward"),wsel=c("select","combine"),
                      dummy=NULL){

  constant <- constant[1]
  search <- match.arg(search,c("backward","forward"))
  wsel <- match.arg(wsel,c("select","combine"))
  type <- match.arg(type,c("logistic","redLogistic","fixLogistic",
                           "doubleExp","exponential","fixExponential",
                           "arctan",
                           "linear"))

  # Find number of variables
  p <- ncol(x)

  # Check the style of the folds input
  if (length(folds)==1){
    # Create cross-validation folds
    cvIndex <- createFolds(y,folds)
  } else {
    cvIndex <- folds
  }

  # Store best results per run
  cvErrorTrace <- list()

  # Set initial model
  if (search == "forward"){
    incl <- c(TRUE,rep(FALSE,p-1))
  } else {
    incl <- rep(TRUE,p)
  }

  cvError <- cvCore(y=y,x=x,incl=incl,type=type,constant=constant,cvIndex=cvIndex,sign=sign,dummy=dummy)
  mdlTrace <- c(incl,mean(cvError),median(cvError),sd(cvError))
  cvErrorTrace[[1]] <- mdlTrace

  # Go into the variable selection if there is something to select
  if (p>1){
    if (search == "forward"){
      reject <- rep(TRUE,p-1)
    } else {
      reject <- rep(FALSE,p-1)
    }
    finish <- FALSE
  } else {
    finish <- TRUE
  }

  # We explore the variables in a backward manner
  # Check whether the CV-error improves when we remove any of the
  # variables
  while (!finish){

    # Create an index of what variables to explore
    if (search == "forward"){
      inclAll <- cbind(1, diag(p-1))
      inclAll <- inclAll == 1
      inclAll <- inclAll | matrix(rep(c(TRUE,!reject),p-1),nrow=p-1,byrow=TRUE)
      inclAll <- inclAll[reject,,drop=FALSE]
    } else {
      inclAll <- cbind(1,matrix(rep(0,sum(reject)*sum(!reject)),nrow=sum(!reject)),-1*(diag(sum(!reject))-1))
      inclAll <- inclAll[,c(1,which(reject)+1,setdiff(2:p,which(reject)+1)),drop=FALSE]
      inclAll <- inclAll == 1
    }

    # Save the errors and setup for this run
    cvErrorGroup <- list()
    for (j in 1:nrow(inclAll)){
      # Get CV error
      inclTemp <- inclAll[j,]
      cvError <- cvCore(y=y,x=x,incl=inclTemp,type=type,constant=constant,cvIndex=cvIndex,sign=sign,reps=reps,dummy=dummy)
      mdlTrace <- c(inclTemp,mean(cvError),median(cvError),sd(cvError))
      cvErrorGroup[[j]] <- mdlTrace
    }
    # Find the minimum error within the run
    cvErrorGroup <- t(abind(cvErrorGroup,along=2))
    minIdx <- which.min(cvErrorGroup[,p+1])
    cvErrorTrace[[length(cvErrorTrace)+1]] <- cvErrorGroup[minIdx,]

    # print(abind(cvErrorTrace,along=2))

    # Compare with previous run to see if we improved by rejecting
    # a variable
    cvErrorTemp <- abind(tail(cvErrorTrace,2),along=2)

    # Continue only if error keeps improving
    if (cvErrorTemp[p+1,1] <= cvErrorTemp[p+1,2]){
      finish <- TRUE
    } else {
      reject <- !cvErrorTemp[2:p,2]
    }

    # or stop when all variables are evaluates
    if (search == "forward"){
      if (all(!reject)){
        finish <- TRUE
      }
    } else {
      if (all(reject)){
        finish <- TRUE
      }
    }
  }

  # Get best model setup and cross-validated errors
  sumStats <- abind(cvErrorTrace,along=2)
  selIdx <- which.min(sumStats[p+1,])
  keep <- sumStats[1:p,selIdx] == 1
  errors <- sumStats[(p+1):(p+3),selIdx]
  names(errors) <- c("Mean MSE","Median MSE","Sd MSE")

  # Check if names exist
  clNames <- colnames(x)
  if (!is.null(clNames)){
    names(keep) <- clNames
  }

  return(list(keep=keep,errors=errors))

}

#' Automatically select curve type and explanatory variables
#'
#' Using cross-validation automatically select explanatory variables jointly with curve type.
#' When running \code{cvfit} there is no need to use \code{\link{varselect}} separately.
#'
#' @inheritParams varselect
#' @param parallel Initialise and use parallel processing for the cross-validation. Note that the maximum number of cores that will be used is equal to the alternative types of curves that are evaluated (\code{alltype}).
#' @param usepbapply A logical to indicate whether to use pbapply to report on the parallel calculation progress. Note that pbapply does not support load balancing.
#' @param alltype A vector of the curve types to consider in the selection process.
#'
#' @return Returns a list with the recommended variable selection choice:
#' \itemize{
#'  \item \code{type} the type of selected curve.
#'  \item \code{keep} a logical vector with which variables to keep.
#'  \item \code{varRes} the result from \code{\link{varselect}} for each evaluated curve type.
#'  \item \code{cvIndex} a matrix detailing how the sample was split for the cross-validation. First column is the fold number and second column is the index of the observation.
#' }
#' Use \code{\link{cvfitplot}} to visualise the output.
#'
#' @inherit curve author references
#' @seealso \code{\link{cvfitplot}}, \code{\link{curve}}, and \code{\link{varselect}}.
#'
#' @examples
#' \donttest{
#'   # Use ECB example data
#'   rate <- ecb$rate
#'   x <- ecb$x[,1:3,drop=FALSE]
#'   cvKeep <- cvfit(x,rate,folds=5,alltype=c("logistic","arctan"),parallel=TRUE)
#'   # Print result
#'   print(cvKeep)
#'   # Fit curve with the selected variables
#'   curve(x[,cvKeep$keep,drop=FALSE],rate,type=cvKeep$type)
#' }
#'
#' @export cvfit

# Identify best curve and variables to include
cvfit <- function(x,y,folds=10,constant=c(TRUE,FALSE),sign=NULL,reps=3,
                  parallel=c(FALSE,TRUE),usepbapply=c(FALSE,TRUE),
                  alltype=c("logistic","redLogistic","doubleExp","exponential","arctan","linear"),
                  search=c("backward","forward"),wsel=c("select","combine"),
                  dummy=NULL){

  constant <- constant[1]
  parallel <- parallel[1]
  usepbapply <- usepbapply[1]
  search <- match.arg(search,c("backward","forward"))
  wsel <- match.arg(wsel,c("select","combine"))

  # alltype <- c("logistic","redLogistic","doubleExp","exponential","linear")
  # alltype <- c("exponential","fixExponential")
  t <- length(alltype)

  # Create common cross-validation folds
  cvIndex <- createFolds(y,folds)

  # We will choose the curve approximation and the selection
  # of variables simultaneously
  if (parallel==FALSE){
    # Without parallel
    varRes <- list()
    for (j in 1:t){
      type <- alltype[j]
      # We use a common set of CV indices, so courses are comparable
      varRes[[j]] <- varselect(x,y,type,folds=cvIndex,constant=constant,sign=sign,reps=reps,search=search,wsel=wsel,dummy=dummy)
    }
  } else {
    # With parallel
    cl <- setParallel(TRUE,maxCore=t)
    # Prepare inputs
    parInputs <- list(x=x,y=y,alltype=alltype,folds=cvIndex,constant=constant,sign=sign,reps=reps,search=search,wsel=wsel,dummy=dummy)
    # Run parallel
    iMax <- 1:t
    if (usepbapply){
      pboptions(use_lb=TRUE)
      varRes <- pblapply(iMax,FUN=varselectParallel,parInputs,cl=cl)
    } else {
      varRes <- clusterApplyLB(cl,iMax,varselectParallel,parInputs=parInputs)
    }
    # Close cluster
    setParallel(FALSE,cl)
  }

  names(varRes) <- alltype

  # Get errors
  cvErrors <- abind(lapply(varRes,function(x){x$errors}),along=2)

  # Get model specifications
  keep <- abind(lapply(varRes,function(x){x$keep}),along=2)

  # Find best specification
  minErr <- which.min(cvErrors[1,])
  keepSel <- keep[,minErr]
  typeSel <- alltype[minErr]

  return(list(type=typeSel,keep=keepSel,varRes=varRes,cvIndex=cvIndex))

}

## Internal functions

# The core part of the cross-validation used for curve and variable selection
cvCore <- function(y,x,incl,type,constant,cvIndex,sign,reps=3,wsel=c("select","combine"),dummy=NULL){
  # This function performs one cross-validation
  # with the settings provided in the inputs

  wsel <- match.arg(wsel,c("select","combine"))

  folds <- max(cvIndex[,1])
  x <- x[,incl,drop=FALSE]
  m <- modelSize(x,type,constant,dummy)
  wOptCV <- array(NA,c(folds,m))
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
    wOptCV[i,] <- curveopt(x=xTrn,y=yTrn,type=type,constant=constant,sign=sign,reps=reps,dummy=dTrn)$w
    if (constant){
      yhat <- curvepred(cbind(1,xTst),wOptCV[i,],type,dummy=dTst)
    } else {
      yhat <- curvepred(xTst,wOptCV[i,],type,dummy=dTst)
    }
    cvError[,i] <- mean((yTst - yhat)^2)
  }
  return(cvError)
}

# Create common subsets for the cross-validation
createFolds <- function(y,folds=10){
  cvIndex <- cvFolds(length(y),K=folds)
  cvIndex <- cbind(cvIndex$which,cvIndex$subsets)
  colnames(cvIndex) <- c("Fold","Sample")
  return(cvIndex)
}

varselectParallel <- function(i,parInputs){

  # Spread variables
  x <- parInputs$x
  y <- parInputs$y
  type <- parInputs$alltype[i]
  folds <- parInputs$folds
  constant <- parInputs$constant
  sign <- parInputs$sign
  reps <- parInputs$reps
  search <- parInputs$search
  wsel <- parInputs$wsel
  dummy <- parInputs$dummy

  return(varselect(x=x,y=y,type=type,folds=folds,constant=constant,
                   sign=sign,reps=reps,search=search,wsel=wsel,
                   dummy=dummy))

}

# Start/stop parallel cluster
setParallel <- function(on=c(TRUE,FALSE),cl=NULL,maxCore=NULL){

  on <- on[1]

  if (on){
    # Initialise parallel
    crs <- detectCores()
    if (is.null(maxCore)){
      clSize <- crs-1
    } else {
      clSize <- min(maxCore,crs-1)
    }
    cl <- makeCluster(getOption("cl.cores", clSize))
    writeLines(paste("Running with", clSize, 'cores'))
    # Load packages to cluster
    invisible(clusterCall(cl, function(pkgs) {
      library(curvir)
    }))
  } else {
    stopCluster(cl)
    cl <- NULL
  }
  return(if(!is.null(cl)){cl})
}
