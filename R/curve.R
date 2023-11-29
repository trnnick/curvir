#' Reserve demand curve
#'
#' Fits the reserve demand curve between excess reserves and normalised rates
#'
#' @param y A vector of normalised interest rates.
#' @param x A matrix of explanatory variables. Excess reserve must be the first input.Additional regressor follow (optional).
#' @param type The type of the reserve demand curve. This can be any of \code{logistic}, \code{redLogistic}, \code{fixLogistic}, \code{doubleExp}, \code{exponential}, \code{fixExponential}, \code{arctan}, \code{linear}. See details in \code{\link{curve}}
#' @param dummy Optional input to signify a regime change (vertical shifts in the curve). Must be a vector of equal length to the rows of \code{x}. If not needed use \code{NULL}.
#' @param q Target interval. This is a scalar below 1, for example 0.9 is the 90\% interval. If \code{NULL} then no quantiles are estimated.
#' @param ... Additional arguments passed to optimiser \code{\link{curveopt}}.
#'
#' @note An additional column for the constant is automatically generated, unless requested otherwise.
#'
#' @return Returns a model of class \code{curvir}. This includes
#' \itemize{
#'  \item \code{type} the type of the curve.
#'  \item \code{constant} a logical indicating the use of a constant.
#'  \item \code{w} a list including: \code{mean} the curve parameters for the mean of the curve, \code{upper} and \code{lower} the parameters for the curve at the upper and lower intervals.
#'  \item \code{data} a list including the \code{y}, \code{x}, and \code{dummy} used for the fitting of the curve.
#'  \item \code{mse} the MSE from the fitting of the curve (the mean only).
#'  \item \code{q} the interval used in the fitting of the curve.
#' }
#'
#' @author Nikolaos Kourentzes, \email{nikolaos@kourentzes.com}
#'
#' @seealso \code{\link{predict.curvir}}, \code{\link{plot.curvir}}, and \code{\link{curveopt}}.
#'
#' @references Chen, Z., Kourentzes, N., & Veyrune, R. (2023). \href{https://www.imf.org/en/Publications/WP/Issues/2023/09/01/Modeling-the-Reserve-Demand-to-Facilitate-Central-Bank-Operations-538754}{Modeling the Reserve Demand to Facilitate Central Bank Operations.} IMF Working Papers, 2023(179).
#'
#' @details For a description of the parametric curves, see the provided reference. Below we list their functions:
#' \itemize{
#'  \item \code{logisitc} (Logistic) \deqn{r_i = \alpha + \kappa / (1 - \beta e^{g(\bm{C}_i)}) + \varepsilon_i}
#'  \item \code{redLogistic} (Reduced logistic) \deqn{r_i = \alpha + 1 / (1 - \beta e^{g(\bm{C}_i)}) + \varepsilon_i}
#'  \item \code{fixLogistic} (Fixed logistic) \deqn{r_i = \alpha + 1 / (1 - e^{g(\bm{C}_i)}) + \varepsilon_i}
#'  \item \code{doubleExp} (Double exponential) \deqn{r_i = \alpha + \beta e^{\rho e^{g(\bm{C}_i)}} + \varepsilon_i}
#'  \item \code{exponential} (Exponential) \deqn{r_i = \alpha + \beta e^{g(\bm{C}_i)} + \varepsilon_i}
#'  \item \code{fixExponential} (Fixed exponential) \deqn{r_i = \beta e^{g(\bm{C}_i)} + \varepsilon_i}
#'  \item \code{arctan} (Arctangent) \deqn{r_i = \alpha + \beta \arctan ( g(\bm{C}_i))  + \varepsilon_i}
#'  \item \code{linear} (Linear) \deqn{r_i = g(\bm{C}_i) + \varepsilon_i}
#' }
#' And \eqn{g(\bm{C}) = c + \bm{C} w_g}, where \eqn{\alpha}, \eqn{\beta}, \eqn{\kappa}, \eqn{\rho} are curve parameters,
#' \eqn{c} is a constant togglable by \code{constant}, \eqn{\bm{C}} are the regressors including the excess reserves. \eqn{w_g} their coefficients, and finally \eqn{\varepsilon_i} is the error term of the curve.
#'
#' @examples
#' \dontshow{
#'   rate <- head(ecb$rate,10)
#'   x <- ecb$x[1:10,1,drop=FALSE]
#'   curve(x,rate,rep=1,type="fixExponential")
#' }
#' \donttest{
#'   # Use ECB example data
#'   rate <- ecb$rate
#'   x <- ecb$x[,1,drop=FALSE]
#'   curve(x,rate)
#'
#'   # An arctangent curve
#'   curve(x,rate,type="arctan")
#'  }
#'
#' @export curve

curve <- function(x,y,type="logistic",dummy=NULL,q=NULL,...){

  type <- match.arg(type,c("logistic","redLogistic","fixLogistic",
                           "doubleExp","exponential","fixExponential",
                           "arctan",
                           "linear"))

  # Optimise curve
  wOpt <- curveopt(x=x,y=y,type=type,dummy=dummy,q=NULL,...)

  # Extract arguments from ellipsis and if missing use defaults from curveopt()
  elarguments <- list(...)
  if ((length(elarguments)>0) && any(grepl("constant",names(elarguments)))){
    constant <- elarguments$constant
  } else {
    constant <- eval(formals(curveopt)$constant)[1]
  }

  # Estimate quantiles
  if (!is.null(q)){
    # Get fitted values
    if (constant){
      xhat <- cbind(1,x)
    } else {
      xhat <- x
    }
    yhat <- curvepred(x=xhat,w=wOpt$w,type=type,dummy=dummy)
    # Optimise # Restricted to not cross yhat
    p <- (1-q)/2
    wU <- curveopt(x,y,type=type,q=1-p,winit=wOpt$w,yhat=yhat,dummy=dummy,constant=constant,...)$w
    wL <- curveopt(x,y,type=type,q=p,winit=wOpt$w,yhat=yhat,dummy=dummy,constant=constant,...)$w
  } else {
    wU <- wL <- NULL
  }

  # Save curve model
  data = list(y=y,x=x,dummy=dummy)
  w <- list(mean=wOpt$w,upper=wU,lower=wL)
  model <- list(type=type,constant=constant,w=w,data=data,mse=wOpt$mse,q=q)
  model <- structure(model,class="curvir")

  return(model)

}

#' Predict method for curvir reserve demand models
#'
#' Predicted values based on curvir model object
#'
#' @param object A model fit with \code{\link{curve}}.
#' @param newdata New input data organised as the x matrix in \code{\link{curve}}. If \code{NULL} then the data used to fit the model is re-used.
#' @param newdummy New input dummy organised as the dummy vector in \code{\link{curve}}. If \code{NULL} then the dummy used to fit the model is re-used.
#' @param ... Further arguments (unused)
#'
#' @return Returns a matrix of predicted values. If the model has estimates for intervals then it will provide upper and lower intervals.
#'
#' @inherit curve references author
#' @seealso \code{\link{curve}}.
#'
#' @examples
#' \dontshow{
#'   rate <- head(ecb$rate,10)
#'   x <- ecb$x[1:10,1,drop=FALSE]
#'   fit <- curve(x,rate,rep=1,type="fixExponential")
#'   predict(fit)
#' }
#' \donttest{
#'   # Use ECB example data
#'   rate <- ecb$rate
#'   x <- ecb$x[,1,drop=FALSE]
#'   fit <- curve(x,rate)
#'   predict(fit)
#'
#'   # An example with new data
#'   predict(fit,newdata=tail(x))
#'  }
#'
#' @describeIn predict Predicted values for parametric curves
#' @export
#' @method predict curvir
predict.curvir <- function(object, newdata=NULL, newdummy=NULL,...){

  # Use old data if nothing new is provided
  if (is.null(newdata)){
    newdata <- object$data$x
  }

  if (is.null(newdummy)){
    newdummy <- object$data$dummy
  }

  # Add constant as needed
  if (object$constant){
    newdata <- cbind(1,newdata)
  }

  if (!is.null(newdummy)){
    if (length(newdummy) != nrow(newdata)){
      stop("Length of newdummy must match the rows of newdata.")
    }
  }

  yhatMean <- curvepred(newdata,object$w$mean,object$type,dummy=newdummy)
  colnames(yhatMean) <- "mean"

  # Interval
  if (!is.null(object$w$upper)){
    yhatUpper <- curvepred(newdata,object$w$upper,object$type,dummy=newdummy)
    colnames(yhatUpper) <- paste0("upper ",object$q*100,"%")
  } else {
    yhatUpper <- NULL
  }
  if (!is.null(object$w$lowe)){
    yhatLower <- curvepred(newdata,object$w$lower,object$type,dummy=newdummy)
    colnames(yhatLower) <- paste0("lower ",object$q*100,"%")
  } else {
    yhatLower <- NULL
  }

  return(cbind(yhatMean,yhatLower,yhatUpper))

}

#' Reserve demand curve predicted values
#'
#' Provides the predicted values for the reserve demand curve of choice. For general use prefer the
#' predict() function, which handles the constant internally.
#'
#' @param x A matrix with the inputs. If there is a constant in the estimated curve, then the first column in \code{x} should be the constant. Next column should be the excess reserves. Additional regressor follow (optional).
#' @param w The vector of weights for the desired curve. Estimated using the \code{\link{curveopt}} function.
#' @inheritParams curve
#'
#' @return Returns a vector of the predicted values.
#'
#' @seealso \code{\link{curve}}, and \code{\link{curveopt}}.
#' @inherit curve references details author
#' @export curvepred

# Note: Constant must be the first column in x, and excess reserve the next
curvepred <- function(x,w,type="logistic",dummy=NULL){

  type <- match.arg(type,c("logistic","redLogistic","fixLogistic",
                           "doubleExp","exponential","fixExponential",
                           "arctan",
                           "linear"))

  p <- ncol(x)
  b <- tail(w,p)
  g <- gX(x,b)

  yhat <- switch(type,
                 logistic = logistic(w,g,dummy),
                 redLogistic = redLogistic(w,g,dummy),
                 fixLogistic = fixLogistic(w,g,dummy),
                 doubleExp = doubleExp(w,g,dummy),
                 exponential = exponential(w,g,dummy),
                 fixExponential = fixExponential(w,g,dummy),
                 arctan = arctan(w,g,dummy),
                 linear = linear(w,g,dummy)
  )

  return(yhat)

}

#' Optimise curve parameters
#'
#' Finds optimal curve parameters.
#'
#' @inheritParams curve
#' @param constant A logical (\code{TRUE} or \code{FALSE}) whether to include a constant or not.
#' @param reps Number of repetitions for the particle swarm optimisation.
#' @param sign A vector of equal length to the number of additional regressors in \code{x} (excluding the constant (if used) and the excess reserves) of positive and negative values (any) that will be used to obtain signs to restrict the values of the estimated parameters. Use \code{NULL} for no restrictions.
#' @param q The desired quantile to optimise for. Use \code{NULL} to get the conditional expectation.
#' @param winit A vector of initial values for the optimisation. This will also carry over to sign restrictions if \code{sameSign==TRUE}.
#' @param yhat Useful when estimating quantiles. Supply here the predicted values for the conditional expectation to add restrictions for the quantiles to not cross the conditional expectation. Use \code{NULL} to not add any restrictions.
#' @param wsel Use the minimum error set of parameters (\code{select}) or the combination of the pool parameters using the heuristic in Kourentzes et al., (2019) (\code{combine}).
#' @param sameSign Used if \code{winit != NULL} to take any sign restrictions from \code{winit}.
#'
#' @return Returns a list of
#' \itemize{
#'  \item \code{w} The optimal parameters
#'  \item \code{mse} The Mean Squared Error of the fitted curve.
#' }
#'
#' @seealso \code{\link{curve}}, and \code{\link{curvepred}}.
#' @inherit curve author
#' @references \itemize{
#'  \item Chen, Z., Kourentzes, N., & Veyrune, R. (2023). \href{https://www.imf.org/en/Publications/WP/Issues/2023/09/01/Modeling-the-Reserve-Demand-to-Facilitate-Central-Bank-Operations-538754}{Modeling the Reserve Demand to Facilitate Central Bank Operations.} IMF Working Papers, 2023(179).
#'  \item Kourentzes, N., Barrow, D., & Petropoulos, F. (2019). Another look at forecast selection and combination: Evidence from forecast pooling. International Journal of Production Economics, 209, 226-235.
#' }
#'
#' @examples
#' \dontshow{
#'   rate <- head(ecb$rate,10)
#'   x <- ecb$x[1:10,1,drop=FALSE]
#'   curveopt(x,rate,reps=1,type="fixExponential")
#' }
#' \donttest{
#'   # Use ECB example data
#'   rate <- ecb$rate
#'   x <- ecb$x[,1,drop=FALSE]
#'   curveopt(x,rate)
#' }
#'
#' @export curveopt

# Function to optimise curve parameters
curveopt <- function(x,y,type="logistic",constant=c(TRUE,FALSE),
                     reps=3,sign=NULL,q=NULL,winit=NULL,yhat=NULL,
                     wsel=c("select","combine"),dummy=NULL,
                     sameSign=c(TRUE,FALSE)){

  constant <- constant[1]
  sameSign <- sameSign[1]
  wsel <- match.arg(wsel,c("select","combine"))

  # Optimise coefficients

  # Clean data
  idx <- !(apply(x,1,function(x){any(is.na(x))}) | is.na(y))
  xClean <- x[idx,,drop=FALSE]
  xClean <- as.matrix(cbind(xClean))
  yClean <- y[idx]
  if (!is.null(dummy)){
    dummy <- dummy[idx]
  }

  # # Order data
  # idx <- order(xClean[,1])
  # xClean <- xClean[idx,,drop=FALSE]
  # yClean <- yClean[idx]

  # Normalise variance - needed for better optimisation
  sdx <- apply(xClean,2,sd)
  xCleanSc <- xClean / matrix(rep(sdx,nrow(xClean)),ncol=ncol(xClean),byrow=TRUE)
  # xClean[,1] <- xClean[,1] / sdx[1]

  # Introduce constant
  if (constant){
    X <- cbind(1,xCleanSc)
  }

  p <- ncol(X) # p includes the constant
  eta <- as.numeric(!is.null(dummy))

  if (is.null(winit)){
    winit <- switch(type,
                    logistic = c(0,1,0,rep(-0.1,p+eta)),
                    redLogistic = c(0,1,rep(-0.1,p+eta)),
                    fixLogistic = c(0,rep(-0.1,p+eta)),
                    doubleExp = c(0.2,0,-0.1,1,rep(0,p-1+eta)),
                    exponential = c(0.2,0,1,rep(0,p-1+eta)),
                    fixExponential = c(0.2,1,rep(0,p-1+eta)),
                    arctan = c(1,-1,rep(0,p+eta)),
                    linear = rep(0,p+eta)
    )
  } else {
    winit <- winit
    # If a winit is provided and sameSign is TRUE than retain
    # signs for the X's
    if (sameSign){
      sign <- tail(sign(winit),p-constant)
    }
  }

  # Restrict for signs of explanatory variables
  upper <- rep(20,length(winit))
  lower <- rep(-20,length(winit))
  if (!is.null(sign)){
    upper[length(winit) - (ncol(x) + constant + eta) + 1 + which(sign < 0)] <- 0
    lower[length(winit) - (ncol(x) + constant + eta) + 1 + which(sign > 0)] <- 0
  }

  # Pre-optmise using Nelder-Meade
  opt1 <- optim(winit,loss,y=yClean,x=X,type=type,q=q,yhat=yhat,dummy=dummy)

  # Override the solution of pre-optimisation if it violates the bounds
  opt1$par[opt1$par>upper] <- upper[opt1$par>upper]
  opt1$par[opt1$par<lower] <- lower[opt1$par<lower]

  # We will optimise multiple times to reduce the stochasticity of the results
  wOptAll <- array(NA,c(reps,length(winit)))
  fnAll <- vector("numeric",reps)
  for (r in 1:reps){
    # Optimise using particle swarm
    psopts <- list(maxit=4000,type="SPSO2007",fnscale=1,s=20)
    opts <- psoptim(opt1$par, loss, upper=upper, lower=lower, control=psopts,
                    y=yClean ,x=X ,type=type, q=q, yhat=yhat, dummy=dummy)

    wOptAll[r,] <- opts$par
    fnAll[r] <- opts$value
  }

  repsRnk <- order(fnAll)
  wOptAll <- wOptAll[repsRnk,,drop=FALSE]
  fnAll <- fnAll[repsRnk]

  if (wsel == "combine"){
    # Use islands to select what to combine
    if (reps != 1){

      # Get correlation from min MSE solution
      crit <- -1*(cor(t(wOptAll))[1,])
      # crit <- c(0,diff(fnAll))

      # Re-order based on critical value
      critRnk <- order(crit)
      crit <- crit[critRnk]
      wOptAll <- wOptAll[critRnk,,drop=FALSE]
      fnAll <- fnAll[critRnk]

      # Combination islands
      critD <- c(0,diff(crit))
      # Create a list of increasing elements
      crt <- lapply(1:length(critD),function(x){critD[1:x]})
      s <- 1.5 # Factor to scale IQR
      th <- sapply(crt,function(x){s*IQR(x)+quantile(x,0.75)})
      th[1] <- NA
      # Find where the threshold is crossed
      cKeep <- which(critD>th)
      if (length(cKeep)==0){
        cKeep <- reps
      } else {
        cKeep <- cKeep[1]-1
      }
      cKeep <- max(cKeep,1)

      # Override for low correlation
      cKeep <- min(cKeep,tail(which(cor(t(wOptAll))[1,]>=0.8),1))

      # plot(X[,1+constant],yClean)
      # cmpT <- colorRampPalette(RColorBrewer::brewer.pal(9,"Reds")[3:8])(reps)
      # for (i in 1:reps){
      #   ord <- order(X[,1+constant])
      #   temp <- curvepred(X,wOptAll[i,],type=type)
      #   lines(X[ord,1+constant],temp[ord],col=cmpT[i])
      # }
      # legend("topright",paste0("rep",1:reps),col=cmpT,bty="n",cex=0.7,lty=1)

      # plot(log(critD),type="o",pch=20)
      # lines(log(th),col="red")
      # abline(v=cKeep)

      # Do a weighted combination on fnAll (MSE)
      keepIdx <- c(rep(TRUE,cKeep),rep(FALSE,reps-cKeep))
      temp <- wOptAll[keepIdx,,drop=FALSE]
      # wOpt <- as.numeric(rbind((1/(crit[keepIdx]+1.1))/sum(1/(crit[keepIdx]+1.1))) %*% temp)
      wOpt <- as.numeric(rbind(1/fnAll[keepIdx]/sum((1/fnAll[keepIdx]))) %*% temp)
    } else {
      # cKeep <- 1
      wOpt <- wOptAll
    }
  } else {
    # Select only top performing set of coefficients
    wOpt <- wOptAll[which.min(fnAll),]
  }

  # Re-scale parameters
  wOpt[(length(wOpt)-ncol(xClean)+1):length(wOpt)] <- wOpt[(length(wOpt)-ncol(xClean)+1):length(wOpt)] / sdx
  # wOpt[(length(wOpt)-ncol(xClean)+1)] <- wOpt[(length(wOpt)-ncol(xClean)+1)] / sdx[1]

  if (constant){
    xhat <- cbind(1,xClean)
  }

  # ord <- order(xhat[,1+constant])
  # plot(xhat[ord,1+constant],yClean[ord])
  # lines(xhat[ord,1+constant],
  #       curvepred(xhat[ord,,drop=FALSE],wOpt,type),col="red")

  # Calculate in-sample MSE on original scale
  mse <- mean((yClean - curvepred(xhat[,,drop=FALSE],wOpt,type,dummy))^2)

  # Name parameters
  wOpt <- namePar(wOpt,type,xClean,constant,dummy)

  return(list(w=wOpt,mse=mse))

}

#' Calculate the inverse curve prediction
#'
#' Calculate the predicted reserves given some rate, i.e., calculate the prediction of the inverse curve.
#'
#' @param object A model fit with \code{\link{curve}}.
#' @param ynew The input rate. If \code{NULL} this corresponds to the values from \code{predict(object)}.
#' @param xnew The values for the additional regressors that were used in the curve fit. Must be a matrix, ordered (columns) as they were input in the fitting of the curve. The constant is dealt with automatically. Do not input the excess reserves. If \code{NULL} this is picked up from the data used to fit the curve.
#' @param dummynew The values for the indicator, if one was used in the fitting of the curve. If \code{NULL} then the data used in the fitting of the model are used.
#' @param warn A logical (\code{TRUE} or \code{FALSE}) to issue a warning if the resulting values are more than 10\% away from the min-max of the excess reserves used to estimate the curve.
#'
#' @return Returns a vector of values of the predicted reserves
#'
#' @seealso \code{\link{curve}}, and \code{\link{predict}}.
#' @inherit curve author references
#'
#' @examples
#' \dontshow{
#'   rate <- ecb$rate
#'   x <- ecb$x[,1,drop=FALSE]
#'   fit <- curve(x,rate,rep=1,type="fixExponential")
#'   invcurve(fit)
#' }
#' \donttest{
#'   # Use ECB example data
#'   rate <- ecb$rate
#'   x <- ecb$x[,1,drop=FALSE]
#'   fit <- curve(x,rate,type="logistic")
#'   invcurve(fit)
#'
#'   # Use a different input rate
#'   invcurve(fit,ynew=0.1)
#' }
#'
#' @export invcurve

invcurve <- function(object,ynew=NULL,xnew=NULL,dummynew=NULL,warn=c(TRUE,FALSE)){

  type <- object$type
  w <- object$w$mean
  constant <- object$constant
  warn <- warn[1]

  # Keep the mean prediction
  yhat <- predict(object)[,1,drop=FALSE]

  if (is.null(ynew)){
    # Use data in the curve object
    ynew <- yhat
    if (is.null(xnew)){
      xnew <- object$data$x
      # Remove the excess reserves
      xnew <- xnew[,-1,drop=FALSE]
    }
    if (is.null(ynew)){
      dummynew <- object$data$dummy
    }
  } else {
    # Data in ynew is external, make sure that the
    # correct inputs are provided
    p <- ncol(object$data$x)-1 # remove the excess reserves
    if (!is.null(xnew)){
      if (ncol(xnew) != p){
        stop("xnew must be a matrix with the same number of explanatory variables as the matrix used to fit the curve. Do not input the excess reserve here.")
      }
      if (length(ynew) != nrow(xnew)){
        stop("xnew must have equal number of rows as the values in ynew.")
      }
    }
    # Check the new dummy
    if (!is.null(object$data$dummy)){
      if (is.null(dummynew)){
        stop("The curve has been fitted with a dummy. Please input a dummy information.")
      }
      if (length(ynew) != length(dummynew)){
        stop("dummy must have the same length as ynew.")
      }
    }
  }

  invy <- switch(type,
                 logistic = invlogistic(w,ynew),
                 redLogistic = invredLogistic(w,ynew),
                 fixLogistic = invfixLogistic(w,ynew),
                 doubleExp = invdoubleExp(w,ynew),
                 exponential = invexponential(w,ynew),
                 fixExponential = invfixExponential(w,ynew),
                 arctan = invarctan(w,ynew),
                 linear = invlinear(w,ynew)
  )

  # Inverse the g(C) part for given explanatory values
  # First construct the value of all X's without the reserves
  # Obtain coefficients
  if (!is.null(xnew)){
    p <- ncol(xnew) + constant
  } else {
    p <- as.numeric(constant)
  }
  b <- tail(w,p+1) # include +1 for the excess reserves
  # Add constant if needed
  if (constant){
    xnew <- cbind(1,xnew)
  }
  # xalpha is the "constant" that goes in the inverse of gX
  # This incorporates the effect of all explanatories, the constant, and the dummy
  # xbeta is the coefficient for the excess reserves
  xalpha <- xnew %*% b[-(1+constant)]
  if (!is.null(object$data$dummy)){
    eta <- w[4]
    xalpha = eta * dummynew
  }
  xbeta <- b[(1+constant)]
  invyy <- invgX(invy,xalpha,xbeta)

  if (any(invyy < min(object$data$x[,1]) - 0.1*diff(range(object$data$x[,1]))) |
      any(invyy > max(object$data$x[,1]) + 0.1*diff(range(object$data$x[,1]))) &&
      warn){
    warning("Predicted value(s) is more than 10% away from the min/max of the data used for fitting the curve.")
  }

  return(invyy)

}

## Internal functions

#' @export
#' @method print curvir
print.curvir <- function(x, ...){
  type <- x$type

  writeLines(paste0(toupper(substr(x$type,1,1)), tolower(substr(x$type,2,nchar(x$type))), " type reserve demand curve"))
  if (ncol(x$data$x)-1>0){
    writeLines(paste0(ncol(x$data$x)-1, " additional regressors included."))
  }
  if (!is.null(x$dummy)){
    writeLines("Shift indicator included.")
  }
  writeLines("")
  writeLines("Curve parameters:")
  print(x$w$mean)
  if (!is.null(x$q)){
    writeLines("")
    writeLines(paste0("Model has estimates for ",round(x$q*100,0),"% intervals."))
  }

}

# Name curve parameters
namePar <- function(wOpt,type="logistic",x=NULL,constant=c(TRUE,FALSE),dummy=NULL){

  constant <- constant[1]
  type <- match.arg(type,c("logistic","redLogistic","fixLogistic",
                           "doubleExp","exponential","fixExponential",
                           "arctan",
                           "linear"))

  # Add curve parameters names
  wNm <- switch(type,
                logistic = c("alpha","beta","kappa"),
                redLogistic = c("alpha","beta"),
                fixLogistic = c("alpha"),
                doubleExp = c("alpha","beta","rho"),
                exponential = c("alpha","beta"),
                fixExponential = c("beta"),
                arctan = c("alpha","beta"),
                linear = NULL
  )
  # If there is a dummy name eta
  if (!is.null(dummy)){
    wNm <- c(wNm,"eta")
  }
  # If there is a constant, name it
  if (constant){
    wNm <- c(wNm,"constant")
  }
  # Names of X's
  if (!is.null(x)){
    wNm <- c(wNm,colnames(x))
  } else {
    xNm <- paste0("X",1:(length(wOpt)-length(wNm)))
    xNm[1] <- "Excess"
    wNm <- c(wNm,xNm)
  }

  names(wOpt) <- wNm
  return(wOpt)

}

# Loss function for optimising curve parameters
loss <- function(w,y,x,type,q=NULL,yhat=NULL,dummy=NULL){

  # w contains the parameters
  # y is the target
  # x is the explanatories as a matrix
  # q is used to optimise for a desired quantile

  yfit <- curvepred(x,w,type,dummy)
  e <- y - yfit

  if (is.null(q)){
    e <- sum(e*e)
  } else {
    b <- rep(q,length(y))
    b[y < yfit] <- 1 - q
    e <- sum(b * abs(y - yfit))
  }

  # Make quantiles to not cross expectation
  if (!is.null(q) && !is.null(yhat)){
    d <- yfit - yhat
    if (any(sign(q-0.5) * d < 0)){
      e <- (e+1)^5
    }
  }

  # e <- sum(abs(e)) + diff(range(e))

  # # Ensure monotonicity
  # if (!any(is.na(yfit)) & !any(is.infinite(yfit))){
  #   if (any(sign(diff(yfit))<=0)){
  #     e <- (e+1)^5
  #   }
  # }

  # if (is.na(e)){
  #   e <- 10^10
  # }

  return(e)

}


# Various curves follow, call them through curvepred()
logistic <- function(w,g,dummy){

  alpha <- w[1]
  beta <- w[2]
  kappa <- w[3]

  yhat <- alpha + kappa/(1-beta*exp(g))

  if (!is.null(dummy)){
    eta <- w[4]
    yhat <- yhat + eta*dummy
  }

  return(yhat)

}


# Estimate g(X) - the exponent part
gX <- function(x,b){
  # If there is a constant, add it to X before it is fed here
  g <- x %*% b
  return(g)
}

# Various curves follow, call them through curvepred()
logistic <- function(w,g,dummy){

  alpha <- w[1]
  beta <- w[2]
  kappa <- w[3]

  yhat <- alpha + kappa/(1-beta*exp(g))

  if (!is.null(dummy)){
    eta <- w[4]
    yhat <- yhat + eta*dummy
  }

  return(yhat)

}

redLogistic <- function(w,g,dummy){

  alpha <- w[1]
  beta <- w[2]

  yhat <- alpha + 1/(1-beta*exp(g))

  if (!is.null(dummy)){
    eta <- w[3]
    yhat <- yhat + eta*dummy
  }

  return(yhat)

}

fixLogistic <- function(w,g,dummy){

  alpha <- w[1]

  yhat <- alpha + 1/(1-exp(g))

  if (!is.null(dummy)){
    eta <- w[2]
    yhat <- yhat + eta*dummy
  }

  return(yhat)

}

doubleExp <- function(w,g,dummy){

  alpha <- w[1]
  beta <- w[2]
  rho <- w[3]

  yhat <- alpha + beta * exp(rho*exp(g))

  if (!is.null(dummy)){
    eta <- w[4]
    yhat <- yhat + eta*dummy
  }

  return(yhat)

}

exponential <- function(w,g,dummy){

  alpha <- w[1]
  beta <- w[2]

  yhat <- alpha + beta * exp(g)

  if (!is.null(dummy)){
    eta <- w[3]
    yhat <- yhat + eta*dummy
  }

  return(yhat)

}

fixExponential <- function(w,g,dummy){

  beta <- w[1]

  yhat <- beta*exp(g)

  if (!is.null(dummy)){
    eta <- w[2]
    yhat <- yhat + eta*dummy
  }

  return(yhat)

}

linear <- function(w,g,dummy){

  yhat <- g

  if (!is.null(dummy)){
    eta <- w[1]
    yhat <- yhat + eta*dummy
  }

  return(yhat)

}

arctan <- function(w,g,dummy){

  alpha <- w[1]
  beta <- w[2]

  yhat <- alpha + beta*atan(g)

  if (!is.null(dummy)){
    eta <- w[3]
    yhat <- yhat + eta*dummy
  }

  return(yhat)

}

invlogistic <- function(w,ynew){

  alpha <- w[1]
  beta <- w[2]
  kappa <- w[3]

  # Bound by alpha
  ynew[ynew<alpha] <- alpha + 1e-18
  # Deal with the asymptote
  dy <- (alpha - ynew)
  dy[dy == 0] <- -1e-18

  # Find inverse
  yinv <- log(-(-kappa/(dy) - 1)/beta)

  return(yinv)

}

# Estimate g(X) - the exponent part
invgX <- function(x,alpha,beta){
  # If there is a constant, add it to X before it is fed here
  g <- (x - alpha)/beta
  return(g)
}

invredLogistic <- function(w,ynew){

  alpha <- w[1]
  beta <- w[2]

  # Bound by alpha
  ynew[ynew<alpha] <- alpha + 1e-18
  # Deal with the asymptote
  dy <- (alpha - ynew)
  dy[dy == 0] <- -1e-18
  # Find inverse
  yinv <- log(-(-1/(dy) - 1)/beta)

  return(yinv)

}

invfixLogistic <- function(w,ynew){

  alpha <- w[1]

  # Bound by alpha
  ynew[ynew<alpha] <- alpha + 1e-18
  # Deal with the asymptote
  dy <- (alpha - ynew)
  dy[dy == 0] <- -1e-18
  # Find inverse
  yinv <- log(1/(dy) + 1)

  return(yinv)

}

invdoubleExp <- function(w,ynew){

  alpha <- w[1]
  beta <- w[2]
  rho <- w[3]

  # Deal with the asymptote
  # Deal with negatives
  dy <- (ynew-alpha)/beta
  dy[dy <= 0] <- 1e-18
  # Deal with negatives again
  dy <- log(dy)/rho
  dy[dy <= 0] <- 1e-18
  # Find inverse
  yinv <- log(dy)

  return(yinv)

}

invexponential <- function(w,ynew){

  alpha <- w[1]
  beta <- w[2]

  # Deal with the asymptote
  dy <- (ynew-alpha)/beta
  dy[dy <= 0] <- 1e-18
  # Find inverse
  yinv <- log(dy)

  return(yinv)

}

invfixExponential <- function(w,ynew){

  beta <- w[1]
  # Find inverse
  yinv <- log(ynew/beta)

  return(yinv)

}

invlinear <- function(w,ynew){

  yinv <- ynew

  return(yinv)

}

invarctan <- function(w,ynew){

  alpha <- w[1]
  beta <- w[2]

  yinv <- tan((ynew - alpha)/beta)
  return(yinv)

}

# Outputs the number of parameters in a curve
modelSize <- function(x,type,constant=c(TRUE,FALSE),dummy=NULL){

  constant <- constant[1]
  constant <- as.numeric(constant)
  eta <- as.numeric(!is.null(dummy))

  p <- switch(type,
              logistic = ncol(x)+3+constant+eta,
              redLogistic = ncol(x)+2+constant+eta,
              fixLogistic = ncol(x)+1+constant+eta,
              doubleExp = ncol(x)+3+constant+eta,
              exponential = ncol(x)+2+constant+eta,
              fixExponential = ncol(x)+1+constant+eta,
              arctan = ncol(x)+2+constant+eta,
              linear = ncol(x)+constant+eta
  )

  return(p)

}
