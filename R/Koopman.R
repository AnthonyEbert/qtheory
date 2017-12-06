

# Koopman equations

#' Koopman formulas as described in Koopman, 1971
#' @param lambda rate of arrivals
#' @param mu rate of service
#' @param C capacity of system
#' @param breaks number of times at which ode is calculated
#' @param init_prob vector of initial probabilities
#' @examples
#' z <- Koopman_MM1(1.5, 2, 20, 20)
#' colSums(z)
Koopman_MM1 <- function(lambda, mu, C, breaks, init_prob = c(1, rep(0, C - 1))){

  a <- lambda / mu
  output_mat <- matrix(nrow = C, ncol = breaks)

  output_mat[,1] <- init_prob

  for(i in 2:breaks){
    output_mat[,i] <- Koopman_MM1_int(output_mat[,i-1], a)
  }

  return(output_mat)
}



Koopman_MM1_int <- function(x, a){

  A <- (x[1] + x[2])

  m <- length(x) - 1

  y <- rep(NA, m + 1)


  y[1:m] <- A * sapply(c(0:I(m-1)), B, a = a)

  for(i in 2:m){
    y[i] <- y[i] + x[3:I(i+1)] %*% sapply(I(i-2):0, B, a = a)
  }

  y[1:m] <- y[1:m] * exp(-a)

  y[m+1] <- sapply(c(m, m:1), input, q = a) %*% x

  return(y)
}

input <- function(shape, q){stats::pgamma(q, shape)}

B <- function(i,a){
  (a ^ i)/(factorial(i))
}
