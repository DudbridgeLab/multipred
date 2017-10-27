# checking dimensions of x and y
checkMatrixDimensions=function(x,y) {
  xmatrix=as.matrix(x)
  ymatrix=as.matrix(y)

  check=NULL
  if (dim(ymatrix)[1]!=dim(xmatrix)[1])
    check="have different numbers of rows"

  if (dim(ymatrix)[2]!=dim(xmatrix)[2])
    check="have different numbers of columns"

  check
}

# check length of x equals columns of y
checkVectorMatrixDimensions=function(x,y) {
  xvector=as.vector(x)
  ymatrix=as.matrix(y)

  check=NULL
  if (length(xvector)<dim(ymatrix)[2])
    check=("has too few elements")

  if (length(xvector)>dim(ymatrix)[2])
    check="has too many elements"

  check
}

# check for a binary vector
checkBinary=function(y) {
  yvector=as.vector(y)

  check=NULL
  if (sum(yvector==0)+sum(yvector==1) != length(yvector))
    check="contains non-binary elements"

  check
}

