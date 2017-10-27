#' Calculate multi-trait measures of predictive accuracy
#'
#' Calculates several measures of accuracy for a vector of predictors.
#'
#' Event-wise accuracy.
#' Individual-wise accuracy: screening, joint and family-wise

#' @param x Matrix of predictions.  Each row corresponds to an individual, each column to a trait.  Each entry is a continuous value.
#' @param y Matrix of traits.  Each row corresponds to an individual, each column to a trait.  May contain binary events (0/1) or continuous values.
#' @param thresh Vector of thresholds.  For each row of x, an event is predicted for each trait that exceeds the corresponding element of thresh.  These predictions are then compared to the elements of y.
#' If thresh is unspecified, and y contains binary events, an error message is printed.
#' @param weight Weighting matrix for calculating event-wise measures.  Defaults to identity.
#' @param sampleSize Number of random samples to draw when estimating C-statistics.  Defaults to 0, in which case
#' full enumeration of x and y is used.
#' @param target Target outcome vector for calculating joint and family-wise sensitivity and specificity.
#' @param target2 Second target outcome vector for calculating joint and family-wise concordance.
#'
#' @export
multiPredAcc=function(x,y,thresh=NULL,weight=NULL,sampleSize=0,target,target2) {

  # coerce x and y to matrices
  x = as.matrix(x)
  y = as.matrix(y)

  # checking dimensions of x and y
  if (dim(y)[1]!=dim(x)[1]) {
    print("Error: x and y have different numbers of rows")
    return(NULL)
  }
  if (dim(y)[2]!=dim(x)[2]) {
    print("Error: x and y have different numbers of columns")
    return(NULL)
  }

  binary=F
  # detect binary traits in y
  if (sum(y==0)+sum(y==1) == length(y)) binary=T
  if (binary) {
    # check dimension of thresh
    if (is.null(thresh)) {
      print("Error: no thresholds specified in thresh")
      return(NULL)
    }
    if (length(thresh)<dim(y)[2]) {
      print("Error: thresh has fewer elements than y")
      return(NULL)
    }
    if (length(thresh)>dim(y)[2]) {
      print("Error: thresh has more elements than y")
      return(NULL)
    }
    #check dimension of target
    if (is.null(target)) target=rep(0,length(thresh))
    if (length(target)<dim(y)[2]) {
      print("Error: target has fewer elements than y")
      return(NULL)
    }
    if (length(target)>dim(y)[2]) {
      print("Error: target has more elements than y")
      return(NULL)
    }

    # predicted binary traits
    predictedTrait=x
    for(i in 1:dim(x)[1]) predictedTrait[i,]=x[i,]>=thresh

        # event wise measures
    # weight matrix is identity if not specified as input
    if (is.null(weight)) weight=diag(length(thresh))

    # event wise sensitivity
    eventSens = sum((y %*% weight) * predictedTrait) / sum((y %*% weight) * y)

    # event wise specificity
    eventSpec = sum(((1-y) %*% weight) * (1-predictedTrait)) / sum(((1-y) %*% weight) * (1-y))

    # event wise positive predictive value
    eventPPV = sum((predictedTrait %*% weight) * y) / sum((predictedTrait %*% weight) * predictedTrait)

    # event wise negative predictive value
    eventNPV = sum(((1-predictedTrait) %*% weight) * (1-y)) / sum(((1-predictedTrait) %*% weight) * (1-predictedTrait))

    # event wise concordance
    # identity weights
    if (sampleSize==0 & identical(weight,diag(1,dim(y)[2]))) {
      eventC=as.numeric(wilcox.test(x[y==1],x[y==0],alternative="greater")$statistic/sum(y==0)/sum(y==1))
    }
    else {
      eventCnumer=0
      eventCdenom=0
      # complete enumeration
      if (sampleSize==0) {
        for(i in 1:dim(y)[1]) {
          xvector=NULL
          yvector=NULL
          for(j in 1:dim(y)[1]) {
            xvector=x[i,]>x[j,]
            yvector=y[i,]>y[j,]
            eventCnumer=eventCnumer+t(xvector) %*% weight %*% yvector
            eventCdenom=eventCdenom+t(yvector) %*% weight %*% yvector
            xvector=x[i,]<x[j,]
            yvector=y[i,]<y[j,]
            eventCnumer=eventCnumer+t(xvector) %*% weight %*% yvector
            eventCdenom=eventCdenom+t(yvector) %*% weight %*% yvector
          }
        }
      }
      #random sampling
      else {
        isample=sample(dim(y)[1],sampleSize,replace=TRUE)
        jsample=sample(dim(y)[1],sampleSize,replace=TRUE)
        for(i in 1:sampleSize) {
          xvector=x[isample[i],]>x[jsample[i],]
          yvector=y[isample[i],]>y[jsample[i],]
          eventCnumer=eventCnumer+t(xvector) %*% weight %*% yvector
          eventCdenom=eventCdenom+t(yvector) %*% weight %*% yvector
          xvector=x[isample[i],]<x[jsample[i],]
          yvector=y[isample[i],]<y[jsample[i],]
          eventCnumer=eventCnumer+t(xvector) %*% weight %*% yvector
          eventCdenom=eventCdenom+t(yvector) %*% weight %*% yvector
        }
      }
      eventC=as.numeric(eventCnumer/eventCdenom)
    }

    # screening sensitivity
    screeningSens = mean(apply(as.matrix(predictedTrait[apply(y,1,max)==1,]),1,max)==1)

    # screening specificity
    screeningSpec = mean(apply(as.matrix(predictedTrait[apply(y,1,max)==0,]),1,max)==0)

    # screening PPV
    screeningPPV = mean(apply(as.matrix(y[apply(predictedTrait,1,max)==1,]),1,max)==1)

    # screening NPV
    screeningNPV = mean(apply(as.matrix(y[apply(predictedTrait,1,max)==0,]),1,max)==0)

    # screening concordance
    screeningCnumer=0
    screeningCdenom=0
    # complete enumeration
    if (sampleSize==0) {
      # individuals with at least one event
      for(i in which(apply(y,1,max)==1)) {
        xvector=NULL
        yvector=NULL
        # individuals with no events
        for(j in which(apply(y,1,max)==0)) {
          screeningCnumer = screeningCnumer + (sum(x[i,]>x[j,])>0)
          screeningCdenom = screeningCdenom + 1
        }
      }
      screeningC=as.numeric(screeningCnumer/screeningCdenom)
    }
    # random sampling
    else {
      isample=sample(which(apply(y,1,max)==1),sampleSize,replace=TRUE)
      jsample=sample(which(apply(y,1,max)==0),sampleSize,replace=TRUE)
      screeningCnumer=0
      for(i in 1:sampleSize)
        screeningCnumer = screeningCnumer + (sum(x[isample[i],]>x[jsample[i],])>0)
      screeningC = screeningCnumer/sampleSize
    }

    # rows of y matching target
    targetIndex = apply(y,1,function(x) sum(x!=target)==0)

    # joint sensitivity
    jointSens = mean(y[targetIndex,] %*% t(1-predictedTrait[targetIndex,])==0)

    # joint specificity
    jointSpec = mean((1-y[targetIndex,]) %*% t(predictedTrait[targetIndex,])==0)

    # family-wise sensitivity
    familySens = mean(y[targetIndex,] %*% t(predictedTrait[targetIndex,])>0)

    # family-wise specificity
    familySpec = mean((1-y[targetIndex,]) %*% t(1-predictedTrait[targetIndex,])>0)

    # rows of y matching target2
    target2Index = apply(y,1,function(x) sum(x!=target2)==0)

    # concordance
    # complete enumeration
    if (sampleSize==0) {
      # joint concordance
      jointC = 0
      for(i in which(targetIndex)) {
        for(j in which(target2Index)) {
          jointC = jointC + ( t(y[i,]>y[j,]) %*% (x[i,]<x[j,]) + t(y[i,]<y[j,]) %*% (x[i,]>x[j,]) ==0)
        }
      }
      jointC = jointC/(sum(targetIndex)*sum(target2Index))

      #family-wise concordance
      familyC = 0
      for(i in which(targetIndex)) {
        for(j in which(target2Index)) {
          familyC = familyC + ( t(y[i,]>y[j,]) %*% (x[i,]>x[j,]) + t(y[i,]<y[j,]) %*% (x[i,]<x[j,]) >0)
        }
      }
      familyC = familyC/(sum(targetIndex)*sum(target2Index))
    }
    # random sampling
    else {
      isample=sample(which(targetIndex),sampleSize,replace=TRUE)
      jsample=sample(which(target2Index),sampleSize,replace=TRUE)
      # joint concordance
      jointC = 0
      for(i in 1:sampleSize)
        jointC = jointC + ( t(y[isample[i],]>y[jsample[i],]) %*% (x[isample[i],]<x[jsample[i],]) +
                              t(y[isample[i],]<y[jsample[i],]) %*% (x[isample[i],]>x[jsample[i],]) ==0)
      jointC = jointC/sampleSize

      # family-wise concordance
      familyC = 0
      for(i in 1:sampleSize)
        familyC = familyC + ( t(y[isample[i],]>y[jsample[i],]) %*% (x[isample[i],]>x[jsample[i],]) +
                                t(y[isample[i],]<y[jsample[i],]) %*% (x[isample[i],]<x[jsample[i],]) >0)
      familyC = familyC/sampleSize
    }

    # rows of predictedTrait matching target
    targetIndex = apply(predictedTrait,1,function(x) sum(x!=target)==0)

    # joint PPV
    jointPPV = mean(predictedTrait[targetIndex,] %*% t(1-y[targetIndex,])==0)

    # joint NPV
    jointNPV = mean((1-predictedTrait[targetIndex,]) %*% t(y[targetIndex,])==0)

    # family-wise PPV
    familyPPV = mean(predictedTrait[targetIndex,] %*% t(y[targetIndex,])>0)

    # family-wise NPV
    familyNPV = mean((1-predictedTrait[targetIndex,]) %*% t(1-y[targetIndex,])>0)



  }
  else { # not binary
    print("Error: Continuous traits not yet supported")
  }

  list(eventSens=eventSens,
       eventSpec=eventSpec,
       eventPPV=eventPPV,
       eventNPV=eventNPV,
       eventC=eventC,
       screeningSens=screeningSens,
       screeningSpec=screeningSpec,
       screeningPPV=screeningPPV,
       screeningNPV=screeningNPV,
       screeningC=screeningC,
       jointSens=jointSens,
       jointSpec=jointSpec,
       jointPPV=jointPPV,
       jointNPV=jointNPV,
       jointC=jointC,
       familySens=familySens,
       familySpec=familySpec,
       familyPPV=familyPPV,
       familyNPV=familyNPV,
       familyC=familyC)
}

