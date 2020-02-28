GLS <-
function(model, data, omega=NULL){
  n <- nrow(data)
  model <- terms(model, data=data)
  y.name <- as.character(model[[2]])
  x.name <- attr(model,"term")

  if(is.null(omega)){omega <- diag(nrow(data))}

  y <- data[, y.name, drop=FALSE]
  x <- data[, x.name, drop=FALSE]
  ydev <- as.matrix(y - colMeans(y))
  xdev <- as.matrix(x - colMeans(x))
  xdev <- as.matrix(apply(x, 2, function(x){x-mean(x)}))

  somega <- solve(omega)
  b_hat <- solve((t(xdev) %*% somega %*% xdev)) %*% (t(xdev) %*% somega %*% ydev)  
  resid0 <- ydev - xdev %*% b_hat
  sigma2 <- as.numeric(1/n * (t(resid0) %*% somega %*% resid0))
  var_b <- sigma2 * solve(t(xdev) %*% somega %*% xdev)

  se <- sqrt(diag(var_b))
  Z <- t(b_hat) / se
  pv <- pnorm(abs(Z), lower.tail=FALSE)*2 
  logL <- -n/2 * log(sigma2) - 1/2*log(det(omega)) - n/2

  Z <- as.vector(Z); names(Z) <- x.name
  pv <- as.vector(pv); names(pv) <- x.name
  b_hat <- as.vector(b_hat); names(b_hat) <- x.name
  return(list(
    coefficients=data.frame(
      estimate=b_hat,
      std.error=se,
      Z.value=Z,
      p.value=pv
    ),
    variance=var_b,
    logL=logL
  ))
}
