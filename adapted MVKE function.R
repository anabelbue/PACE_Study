MVKE <- function(d, h = 0.2, kernel = c("exp", "Gaussian")) {
  
  # Custom rowProds function to replace Rfast::rowprods
  rowProds <- function(x) {
    apply(x, 1, prod)  # Compute row-wise products
  }
  
  # K_gaussian_mat function
  K_gaussian_mat <- function(mat, x, h) {
    dim <- length(x)
    mat <- mat - matrix(rep(x, nrow(mat)), ncol = dim, byrow = TRUE)
    mat <- stats::dnorm(mat / h)
    values <- 1 / (h^dim) * rowProds(mat)  # Use the custom rowProds function
    return(values)
  }
  
  # K_exp_mat function
  K_exp_mat <- function(mat, x, h) {
    dim <- length(x)
    mat <- mat - matrix(rep(x, nrow(mat)), ncol = dim, byrow = TRUE)
    mat <- exp(mat / h)
    values <- 1 / (h^dim) * rowProds(mat)  # Use the custom rowProds function
    return(values)
  }
  
  # Main MVKE function logic
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
  
  # The returned function
  function(x) {
    if (length(x) != dim) stop("Input of wrong dimension.")
    
    temp_kernel_term_upper <- K(temp_d, x, h = h)
    temp_kernel_term_lower <- K(d, x, h = h)
    
    # Debugging the MVKEresult(xx)$mu issue
    Useq <- numeric(length(x))
    Useq[1] <- 0
    for (i in 2:length(x)) {
      Useq[i] <- Useq[i - 1] - stats::integrate(function(x) purrr::map_dbl(x, function(xx) {
        1  # Temporarily use a simple value to check if error persists
      }), x[i - 1], x[i], subdivisions = 100L, rel.tol = .Machine$double.eps^0.25, abs.tol = .Machine$double.eps^0.25)$value
    }
    
    return(list(
      mu = colSums(temp_kernel_term_upper * temp_diff) / sum(temp_kernel_term_lower),
      a = mapply(`*`, temp_kernel_term_upper, temp_diff_tcrossprod, SIMPLIFY = FALSE) %>% 
        Reduce(`+`, .) / sum(temp_kernel_term_lower)
    ))
  }

