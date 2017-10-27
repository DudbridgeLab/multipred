#' Screening measures of predictive accuracy
#'
#' Calculates screening sensitivity, specificity, positive and negative predictive value, concordance and net benefit for a vector of predictors.
#'
#' Screening measures consider the prediction of at least one event among several, without regard to whether the correct events are predicted.

#' @param targetProb Probability of at least one event.
#' If NULL, which is the default, then targetProb is estimated from the data, ignoring ascertainment.
#' @inheritParams eventWise
#'
#' @export
screening=function(x,y,thresh=NULL,targetProb=NULL,nsample=0) {

  # coerce x and y to matrices
  x = as.matrix(x)
  y = as.matrix(y)

  check=checkMatrixDimensions(x,y)
  if (!is.null(check)) {
    print(paste("x and y",check))
    return(NULL)
  }

  check=checkVectorMatrixDimensions(thresh,y)
  if (!is.null(check)) {
    print(paste("thresh",check))
    return(NULL)
  }

  check=checkBinary(y)
  if (!is.null(check)) {
    print(paste("y",check))
    return(NULL)
  }

  # predicted binary traits
  predictedTrait=x
  for(i in 1:dim(x)[1]) predictedTrait[i,]=x[i,]>=thresh

  # sensitivity
  sens = mean(apply(as.matrix(predictedTrait[apply(y,1,max)==1,]),1,max)==1)

  # specificity
  spec = mean(apply(as.matrix(predictedTrait[apply(y,1,max)==0,]),1,max)==0)

  # positive predictive value
  PPV = mean(apply(as.matrix(y[apply(predictedTrait,1,max)==1,]),1,max)==1)

  # negative predictive value
  NPV = mean(apply(as.matrix(y[apply(predictedTrait,1,max)==0,]),1,max)==0)

  # concordance
  Cnumer=0
  Cdenom=0
  # complete enumeration
  if (nsample==0) {
    # individuals with at least one event
    for(i in which(apply(y,1,max)==1)) {
      xvector=NULL
      yvector=NULL
      # individuals with no events
      for(j in which(apply(y,1,max)==0)) {
        Cnumer = Cnumer + (sum(x[i,]>x[j,])>0)
        Cdenom =Cdenom + 1
      }
    }
    C=as.numeric(Cnumer/Cdenom)
  }
  # random sampling
  else {
    isample=sample(which(apply(y,1,max)==1),nsample,replace=TRUE)
    jsample=sample(which(apply(y,1,max)==0),nsample,replace=TRUE)
    Cnumer=0
    Cdenom=0
    for(i in 1:nsample) { # only consider pairs where one score vector is systematically <= than the other
      if (sum(x[isample[i],]>x[jsample[i],])==0 |
          sum(x[isample[i],]<x[jsample[i],])==0) {
        Cnumer = Cnumer + (sum(x[isample[i],]>x[jsample[i],])>0)
        Cdenom = Cdenom +1
      }
    }
    C = as.numeric(Cnumer/Cdenom)
    print(Cdenom)
  }

  # net benefit
  if (is.null(targetProb)) targetProb=mean(apply(y,1,max))
  NB = sens - (1-spec) * min(thresh/(1-thresh)) * (1-targetProb)/targetProb

  list(sens=sens,
       spec=spec,
       PPV=PPV,
       NPV=NPV,
       C=C,
       NB=NB
  )
}
