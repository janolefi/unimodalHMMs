#' Source functions for the paper "", Jan-Ole Koslik, Fanny Dupont, Marie Auger-Méthé, [whale folks?],  Nancy Heckman


#' Construct unimodal constraint matrix
#' 
#' @param m mode index (vector of length N)
#' @param N number of states
#' @param k basis dimension
#' @return The matrix of unimodal constraints based on Eq ().
#' 
construct_C <- function(m,
                        N, 
                        k)
{
  # construct constraint matrices
  C <- list() # initialise list
  for(i in 1:N){
    Ci <- diff(diag(k)) # initialise with first-order difference matrix
    if(m[i] < k){ # if there is a block that needs to be flipped, do that
      Ci[m[i]:(k-1), ] <- -Ci[m[i]:(k-1),]
    }
    C[[i]] <- Ci[,-k] # exclude last column because multiplied by zero
  }
  return(C)
}

min0_smooth <- function(x, alpha = 20, eps = 0){
  -1 / alpha * log(1 + exp(-alpha * (x + eps)))
}

