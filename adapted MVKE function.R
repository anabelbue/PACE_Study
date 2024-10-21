MVKE <- function(d, h = 0.2, kernel = c("exp", "Gaussian")) {
  rowProds <- function(x) {
    apply(x, 1, prod)
  }
  
  if (is.data.frame(d)) d <- as.matrix(d)
  if (!is.matrix(d)) stop("`d` should be a data.frame or a matrix.")
  
  d <- stats::na.omit(d)
  dim <- ncol(d)
  
  temp_d <- d[1:(nrow(d) - 1), , drop = FALSE]
  temp_diff <- diff(d)
  temp_norm <- apply(temp_diff, MARGIN = 1, FUN = function(x) norm(x, "2"))
  temp_diff_tcrossprod <- apply(temp_diff,
                                MARGIN = 1,
                                FUN = function(x) {
                                  tcrossprod(x, x)
                                }, simplify = FALSE
  )
  kernel <- kernel[1]
  if (kernel == "Gaussian") {
    K <- K_gaussian_mat
  } else if (kernel == "exp") {
    K <- K_exp_mat
  } else {
    stop('`kernel` must be one of "Gaussian" or "exp".')
  }
  
  force(h)
  
  # Debugging and printing output of MVKEresult
  function(x) {
    if (length(x) != dim) stop("Input of wrong dimension.")
    
    temp_kernel_term_upper <- K_gaussian_mat(temp_d, x, h = h)
    temp_kernel_term_lower <- K_gaussian_mat(d, x, h = h)
    
    # Print the MVKEresult for debugging
    print("Debugging MVKEresult call:")
    MVKEresult <- K_gaussian_mat(temp_d, x, h = h)  # Calling the MVKE kernel function
    print(MVKEresult)  # Print the result of MVKE
    print(MVKEresult(x))  # Print the MVKE result for input x
    print(MVKEresult(x)$mu)  # Print the mu component of the result
    
    if (any(is.na(temp_kernel_term_upper)) || any(is.nan(temp_kernel_term_upper))) {
      stop("Kernel term has NA or NaN values.")
    }
    
    # Perform the calculation
    return(list(
      mu = colSums(temp_kernel_term_upper * temp_diff) / sum(temp_kernel_term_lower),
      a = mapply(`*`, temp_kernel_term_upper, temp_diff_tcrossprod, SIMPLIFY = FALSE) %>% 
        Reduce(`+`, .) / sum(temp_kernel_term_lower)
    ))
  }
}

K_gaussian_mat <- function(mat, x, h) {
  dim <- length(x)
  mat <- mat - matrix(rep(x, nrow(mat)), ncol = dim, byrow = TRUE)
  mat <- stats::dnorm(mat / h)
  values <- 1 / (h^dim) * rowProds(mat)  # Use the custom rowProds function
  return(values)
}

K_exp_mat <- function(mat, x, h) {
  dim <- length(x)
  mat <- mat - matrix(rep(x, nrow(mat)), ncol = dim, byrow = TRUE)
  mat <- exp(mat / h)
  values <- 1 / (h^dim) * rowProds(mat)  # Use the custom rowProds function
  return(values)
}

# Custom rowProds function to replace Rfast::rowprods
rowProds <- function(x) {
  apply(x, 1, prod)  # Compute row-wise products
}
