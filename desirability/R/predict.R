predict.dBox <- function(object, newdata = NA, ...)
{
   out <- vector(length = length(newdata), mode = "numeric") * NA 
   out[newdata < object$low | newdata > object$high] <- 0
   out[newdata >= object$low & newdata <= object$high] <- 1
   out[is.na(out)] <- object$missing
   if(!is.null(object$tol)) out[out == 0] <- object$tol
   out
}

predict.dMax <- function(object, newdata = NA, ...)
{
   out <- vector(length = length(newdata), mode = "numeric") * NA 
   out[newdata < object$low] <- 0
   out[newdata > object$high] <- 1
   out[newdata <= object$high & newdata >= object$low] <- (
      (newdata[newdata <= object$high & newdata >= object$low] - object$low)/
      (object$high - object$low))^object$scale 
   out[is.na(out)] <- object$missing      
   out
}

predict.dMin <- function(object, newdata = NA, ...)
{
   out <- vector(length = length(newdata), mode = "numeric") * NA
   out[newdata < object$low] <- 1
   out[newdata > object$high] <- 0
   out[newdata <= object$high & newdata >= object$low] <- (
      (newdata[newdata <= object$high & newdata >= object$low] - object$high)/
      (object$low - object$high))^object$scale 
   out[is.na(out)] <- object$missing      
   out
}


predict.dTarget <- function(object, newdata = NA, ...)
{
   out <- vector(length = length(newdata), mode = "numeric")  * NA
   out[newdata < object$low | newdata > object$high] <- 0
   out[newdata <= object$target & newdata >= object$low] <- (
      (newdata[newdata <= object$target & newdata >= object$low] - object$low)/
      (object$target - object$low))^object$lowScale 
   out[newdata <= object$high & newdata >= object$target] <- (
      (newdata[newdata <= object$high & newdata >= object$target] - object$high)/
      (object$target - object$high))^object$highScale       
   out[is.na(out)] <- object$missing
      
   out
}

predict.dArb <- function(object, newdata = NA, ...)
{
   out <- vector(length = length(newdata), mode = "numeric")  * NA
   out[newdata < min(object$x)] <- object$d[1]
   out[newdata > max(object$x)] <- object$d[length(object$d)]
   

   inBtwn <- newdata >= min(object$x) & newdata <= max(object$x)
   if(any(inBtwn))
   {
      tmp <- matrix(newdata[inBtwn], ncol = 1)
      approxD <- apply(
         tmp, 
         1, 
         function(u, x, y) approx(x, y, u)$y, 
         x = object$x,
         y = object$d)
      out[inBtwn] <- approxD
      
   }
   out[is.na(out)] <- object$missing
   
   out
}

predict.dOverall <- function(object, newdata = matrix(NA, ncol = length(object$d)), all = FALSE, ...)
{
   numD <- length(object$d)
   if(is.vector(newdata)) newdata <- matrix(newdata, ncol = length(newdata))
   if(is.matrix(newdata)) newdata <- as.matrix(newdata)
   if(numD != ncol(newdata)) stop("the number of columns in newdata must match the number of desirability functions")
   
   indD <- newdata * NA
   for(i in 1:numD)
   {
      indD[,i] <- predict(object$d[[i]], newdata[,i])
   }
   overall <- apply(indD, 1, prod)^(1/numD)
   if(all)
   {
      out <- cbind(indD, overall)
      colnames(out)[numD + 1] <- "Overall"
   
   } else out <- overall
   
   out
}
