N_generate_Mix_log <- function(N, A, alpha) {
  r <- ncol(A)
  if(length(alpha) == 1){alpha <- rep(alpha,r)}
  Z <- vector('list', length = r)
  for(k  in 1:r){
    Z[[k]] <- matrix(0,nrow(A),N)
    sig <- which(A[,k]>0)
    Z[[k]][sig,] <- t(rmev(N, length(sig), param = 1/alpha[k], model = "log"))*A[sig,k]
  }
  M <- apply(simplify2array(Z),c(1,2),max)
  return(t(M))
}

# λ = sqrt(Γ)/2. 
# The relationship between λ and r—the dependence parameter of the Hüsler–Reiss (HR) model—
# in the R package `evd` is discussed in the documentation of the `rmev` function.
# Express r in terms of Γ using the stdf expression 
# provided in Section 4.2.2 of https://link.springer.com/article/10.1007/s10687-024-00501-4
# and the details given on page 14 of http://cran.fhcrc.org/web/packages/evd/evd.pdf.
N_generate_Mix_hr <- function(N, A, sigma) {
  r <- ncol(A)
  Z <- vector('list', length = r)
  for(k  in 1:r){
    Z[[k]] <- matrix(0,nrow(A),N)
    sig <- which(A[,k]>0)
    lsig <- length(sig)
    if(lsig == 1){ Z[[k]][sig , ] <- -1/ log( runif(N , 0 , 1) ) * A[sig,k]  } 
    else{
      sub_sigma <- sigma[sig , sig]
      Z[[k]][sig,] <- t(rmev(N, lsig, sigma = sub_sigma, model = "hr"))*A[sig,k] 
    }
  }
  M <- apply(simplify2array(Z),c(1,2),max)
  return(t(M))
}

