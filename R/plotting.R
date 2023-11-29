#' Plot method for curvir reserve demand models
#'
#' Plot a reserve demand curve estimated using \code{\link{curve}}.
#'
#' @param x A model fit with \code{\link{curve}}.
#' @param ploty Logical (\code{TRUE} or \code{FALSE}) to plot data used for fitting or not.
#' @param plotq Logical (\code{TRUE} or \code{FALSE}) to plot intervals.
#' @param usemean Logical (\code{TRUE} or \code{FALSE}) to replace all explanatory variables (except from excess reserves) with their mean. It cannot be used if a \code{dummy} is used in the fitting of the curve.
#' @param prcmp Logical (\code{TRUE} or \code{FALSE}) to use principal components to project all variables into two dimensions. Requires \code{usemean == FALSE}. This may be useful when multiple explanatory variables are used.
#' @param useline Logical (\code{TRUE} or \code{FALSE}) to force line plots when \code{usemean == FALSE} or \code{prcomp == FALSE} and there are many explanatory variables.
#' @param main Use to overide the default plot title.
#' @param pch Use to overide the dafault marker type for points.
#' @param ... Additional inputs (currently unused).
#'
#' @return No returned value. Produces a plot of the estimated curve.
#'
#' @inherit curve references author
#' @seealso \code{\link{curve}}.
#'
#' @examples
#' \dontshow{
#'   rate <- head(ecb$rate,10)
#'   x <- ecb$x[1:10,1,drop=FALSE]
#'   fit <- curve(x,rate,rep=1,type="fixExponential")
#'   plot(fit)
#' }
#' \donttest{
#'   # Use ECB example data
#'   rate <- ecb$rate
#'   x <- ecb$x[,1,drop=FALSE]
#'   fit <- curve(x,rate)
#'   plot(fit)
#' }
#'
#' @describeIn plot Plot parametric curves
#' @export
#' @method plot curvir
plot.curvir <- function(x,ploty=c(FALSE,TRUE),plotq=c(TRUE,FALSE),
                        usemean=c(FALSE,TRUE),prcmp=c(FALSE,TRUE),
                        useline=c(FALSE,TRUE),main=NULL,pch=20,...){

  # Replace all X's with their mean (for plotting).
  # This cannot be used when dummy is used.
  usemean <- usemean[1]
  # If usemean=FALSE and ncol(x)>1 then use principal components for plotting
  prcmp <- prcmp[1]
  # If there are many X's (or a dummy) and usemean or prcmp is not used, then
  # order points and connect with a line.
  useline <- useline[1]
  # Plot actuals
  ploty <- ploty[1]
  # Plot intervals
  plotq <- plotq[1]

  # Retrieve data from model
  xdat <- x$data$x
  dummy <- x$data$dummy
  type <- x$type
  constant <- x$constant

  if (ploty){
    y <- x$data$y
  } else {
    y <- NULL
  }

  # Predict output
  if (ncol(xdat)>1 | is.null(dummy)){
    if (usemean & is.null(dummy)){
      xhat <- seq(min(xdat[,1]),max(xdat[,1]),length.out=1000)
      xhat <- cbind(xhat,matrix(rep((colMeans(xdat))[-1],length(xhat)),nrow=length(xhat),byrow=TRUE))
      colnames(xhat) <- colnames(xdat)
    } else {
      xhat <- xdat
    }
    if (constant){
      xhat <- cbind(1,xhat)
    }
  } else {
    if (is.null(dummy)){
      xhat <- cbind(seq(min(xdat),max(xdat),length.out=1000))
    } else {
      xhat <- cbind(xdat)
    }

    if (constant){
      xhat <- cbind(1,xhat)
    }
  }

  yhat <- curvepred(x=xhat,w=x$w$mean,type=type,dummy=dummy)
  if (plotq && (!is.null(x$q))){
    yhatQU <- curvepred(x=xhat,w=x$w$upper,type=type,dummy=dummy)
    yhatQL <- curvepred(x=xhat,w=x$w$lower,type=type,dummy=dummy)
  } else {
    yhatQU <- yhatQL <- NULL
  }

  # Produce plot
  if (is.null(main)){
    mtext <- type
  } else {
    mtext <- main
  }
  plotCore(xdat,xhat,y,yhat,yhatQL,yhatQU,dummy,constant,main=mtext,pch,usemean,prcmp,useline)

}

#' Plot method for npcurvir reserve demand models
#'
#' Plot a non-parametric reserve demand curve estimated using \code{\link{npcurve}}.
#'
#' @inherit curve references author
#' @seealso \code{\link{npcurve}}.
#'
#' @examples
#' # Use ECB example data
#' rate <- ecb$rate
#' x <- ecb$x[,1,drop=FALSE]
#' fit <- npcurve(x,rate)
#' plot(fit)
#'
#' @describeIn plot Plot non-parametric curves
#' @export
#' @method plot npcurvir
plot.npcurvir <- function(x,ploty=c(FALSE,TRUE),plotq=c(TRUE,FALSE),
                        usemean=c(FALSE,TRUE),prcmp=c(FALSE,TRUE),
                        useline=c(FALSE,TRUE),main=NULL,pch=20,...){

  # Replace all X's with their mean (for plotting).
  # This cannot be used when dummy is used.
  usemean <- usemean[1]
  # If usemean=FALSE and ncol(x)>1 then use principal components for plotting
  prcmp <- prcmp[1]
  # If there are many X's (or a dummy) and usemean or prcmp is not used, then
  # order points and connect with a line.
  useline <- useline[1]
  # Plot actuals
  ploty <- ploty[1]
  # Plot intervals
  plotq <- plotq[1]

  # Retrieve data from model
  xdat <- x$data$x
  dummy <- x$data$dummy
  type <- x$type

  if (ploty){
    y <- x$data$y
  } else {
    y <- NULL
  }

  # Predict output
  if (ncol(xdat)>1 | is.null(dummy)){
    if (usemean & is.null(dummy)){
      xhat <- seq(min(xdat[,1]),max(xdat[,1]),length.out=1000)
      xhat <- cbind(xhat,matrix(rep((colMeans(xdat))[-1],length(xhat)),nrow=length(xhat),byrow=TRUE))
      colnames(xhat) <- colnames(xdat)
    } else {
      xhat <- xdat
    }
  } else {
    if (is.null(dummy)){
      xhat <- cbind(seq(min(xdat),max(xdat),length.out=1000))
    } else {
      xhat <- cbind(xdat)
    }
  }

  yhat <- predict(x,newdata=xhat,newdummy=dummy)

  if (plotq && (!is.null(x$q))){
    yhatQU <- yhat[,3]
    yhatQL <- yhat[,2]
  } else {
    yhatQU <- yhatQL <- NULL
  }
  yhat <- yhat[,1]

  # Produce plot
  if (is.null(main)){
    mtext <- type
  } else {
    mtext <- main
  }
  plotCore(xdat,xhat,y,yhat,yhatQL,yhatQU,dummy,constant=FALSE,main=mtext,pch,usemean,prcmp,useline)

}

#' Plot output of cvfit to facilitate comparison
#'
#' Plot summarised error of different curves specification from \code{\link{cvfit}}.
#' Assuming normal errors, plot the mean cross-validated error and the 95% interval around it for each curve and its best selection of variables.
#'
#' @param cvKeep The output of \code{\link{cvfit}}.
#' @param xlock Focus horizontal axis on \code{mean} or \code{median}, The latter provides a narrower range, facilitating comparison when there are outliers.
#' @param cvRF Include cross-validation results from random forecast non-parametric curves. Obtain these from \code{\link{cvnpcurve}} with \code{type="rforest"}. Use \code{NULL} to ignore.
#' @param cvSP Include cross-validation results from spline regression non-parametric curves. Obtain these from \code{\link{cvnpcurve}} with \code{type="spline"}. Use \code{NULL} to ignore.
#'
#' @return No returned value. Produces a summary plot of cross-validated errors.
#'
#' @inherit curve references author
#' @seealso \code{\link{cvfit}}.
#'
#' @examples
#' \donttest{
#'   # Use ECB example data
#'   rate <- ecb$rate
#'   x <- ecb$x[,1:3,drop=FALSE]
#'   cvKeep <- cvfit(x,rate,folds=5,alltype=c("logistic","arctan"),parallel=TRUE)
#'   cvfitplot(cvKeep)
#'   # Add results from non-parameteric curves
#'   cvRF <- cvnpcurve(x,rate,cvKeep$cvIndex)
#'   cvSP <- cvnpcurve(x,rate,cvKeep$cvIndex,type="spline")
#'   cvfitplot(cvKeep,cvRF=cvRF,cvSP=cvSP)
#' }
#'
#' @export cvfitplot

cvfitplot <- function(cvKeep,xlock=c("mean","median"),cvRF=NULL,cvSP=NULL){

  # Ensure graceful exit
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar),add = TRUE)

  # Check inputs
  xlock <- match.arg(xlock,c("mean","median"))

  # Increase horizontal space of plot to accommodate labels
  vmar <- oldpar$mar
  vmar[2] <- vmar[2] + 2
  vmar[4] <- vmar[4] - 1
  par(mar=vmar)

  pErr <- abind(lapply(cvKeep$varRes,function(x){x$errors}),along=2)
  if (!is.null(cvRF)){
    pErr <- cbind(pErr,cvRF)
    colnames(pErr) <- c(colnames(pErr)[1:(ncol(pErr)-1)],"rforest")
  }
  if (!is.null(cvSP)){
    pErr <- cbind(pErr,cvSP)
    colnames(pErr) <- c(colnames(pErr)[1:(ncol(pErr)-1)],"spline")
  }
  m <- ncol(pErr)
  f <- qnorm(0.975)
  s <- which.min(pErr[1,])
  eRng <- rbind(pErr[1,]-f*pErr[3,],pErr[1,],pErr[1,]+f*pErr[3,])
  eRng[eRng<0] <- 0
  if (xlock=="mean"){
    xx <- range(eRng)
  } else {
    xx <- c(min(eRng),max(pErr[2,]))
  }
  xx <- xx +c(-1,1)*0.04*diff(xx)
  xx[xx<0] <- 0
  plot(NA,NA,xlim=xx,ylim=c(0,m+1.5),xlab="",ylab="",yaxt="n")
  for (i in 1:m){
    if (i == s){
      cmp <- RColorBrewer::brewer.pal(3,"Set1")[2]
    } else {
      cmp <- "black"
    }
    lines(eRng[c(1,3),i],c(i,i),col=cmp)
    lines(rep(eRng[1,i],2),i+c(-1,1)*0.1,col=cmp)
    lines(rep(eRng[3,i],2),i+c(-1,1)*0.1,col=cmp)
    points(eRng[2,i],i,pch=21,bg=if(i==s){cmp}else{"white"},lwd=2,col=cmp)
    points(pErr[2,i],i,pch=4,lwd=2,col=cmp)
  }
  legend("topleft",c("Mean","Median"),pch=c(1,4),horiz=TRUE,bty="n",lwd=2,lty=NA)
  title(xlab="CV MSE",line=2.5)
  axis(2,at=1:m,labels=colnames(pErr),las=1)

}

## Internal functions

# Core function for plots of individual curves
plotCore <- function(x,xhat,y,yhat,yhatQL,yhatQU,dummy,constant,
                     main,pch,usemean,prcmp,useline){

  # Produce plot
  cmp <- brewer.pal(3,"Set1")
  yy <- range(c(yhat,yhatQL,yhatQU,y))
  yy <- yy + c(-1,1)*0.04*diff(yy)
  xx <- range(x[,1])
  xx <- xx + c(-1,1)*0.04*diff(xx)
  if (!prcmp | !(ncol(x)>1 && usemean==FALSE)){
    plot(NA,NA,xlim=xx,ylim=yy,xlab="",ylab="",xaxt="n")
    xTck <- axTicks(1)
    axis(1,at=xTck,labels=formatC(xTck,format="d"))
    if (!is.null(y)){
      points(x[,1],y,pch=20)
    }
  }

  if (ncol(x)>1 && usemean==FALSE){
    if (prcmp){
      if (!is.null(dummy)){
        pc <- prcomp(cbind(x[,],dummy))
      } else {
        pc <- prcomp(x[,])
      }
      idx <- order(pc$x[,1])
      if (!is.null(y)){
        plot(pc$x[idx,1],y[idx],type="p",pch=20,ylab="",xlab="")
      }
      lines(pc$x[idx,1],yhat[idx],col=cmp[1])
      if (!is.null(yhatQU)){
        lines(pc$x[idx,1],yhatQU[idx],col=cmp[2])
      }
      if (!is.null(yhatQL)){
        lines(pc$x[idx,1],yhatQL[idx],col=cmp[2])
      }
    } else {
      if (useline){
        idx <- order(xhat[,1+constant])
        lines(xhat[idx,1+constant],yhat[idx],col=cmp[1])
        if (!is.null(yhatQU)){
          lines(xhat[idx,1+constant],yhatQU[idx],col=cmp[2])
        }
        if (!is.null(yhatQL)){
          lines(xhat[idx,1+constant],yhatQL[idx],col=cmp[2])
        }
      } else {
        points(xhat[,1+constant],yhat,col=cmp[1],pch=pch)
        if (!is.null(yhatQU)){
          points(xhat[,1+constant],yhatQU,col=cmp[2],pch=pch)
        }
        if (!is.null(yhatQL)){
          points(xhat[,1+constant],yhatQL,col=cmp[2],pch=pch)
        }
      }
    }
  } else {
    if (is.null(dummy) | useline){
      idx <- order(xhat[,1+constant])
      lines(xhat[idx,1+constant],yhat[idx],col=cmp[1],lwd=1.5)
      if (!is.null(yhatQU)){
        lines(xhat[idx,1+constant],yhatQU[idx],col=cmp[2],lwd=1.5)
      }
      if (!is.null(yhatQL)){
        lines(xhat[idx,1+constant],yhatQL[idx],col=cmp[2],lwd=1.5)
      }
    } else {
      points(xhat[,1+constant],yhat,col=cmp[1],pch=pch)
      if (!is.null(yhatQU)){
        points(xhat[,1+constant],yhatQU,col=cmp[2],pch=pch)
      }
      if (!is.null(yhatQL)){
        points(xhat[,1+constant],yhatQL,col=cmp[2],pch=pch)
      }
    }
  }
  title(main=main,font.main=1,line=1)

}
