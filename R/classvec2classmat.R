classvec2classmat <- function(yvec) #, transpose=FALSE)
{
  yvec <- factor(yvec)
  nclasses <- nlevels(yvec)

  outmat <- matrix(0, length(yvec), nclasses)
  dimnames(outmat) <- list(NULL, levels(yvec))
  
  for (i in 1:nclasses)
    outmat[which(as.integer(yvec) == i),i] <- 1

#  if (transpose) outmat <- t(outmat)

  outmat
}

classmat2classvec <- function(ymat, threshold=0) #, transpose=FALSE)
{
#  if (transpose) ymat <- t(ymat)
  
  class.names <- dimnames(ymat)[[2]]
  if (is.null(class.names)) class.names <- 1:ncol(ymat)

  classes <- apply(ymat, 1, function(x) which(x == max(x))[1])
  classes[apply(ymat, 1, max) < threshold] <- NA
  
  class.names[classes]
}
