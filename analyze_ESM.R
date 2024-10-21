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
  
  function(x) {
    if (length(x) != dim) stop("Input of wrong dimension.")
    
    temp_kernel_term_upper <- K_gaussian_mat(temp_d, x, h = h)
    temp_kernel_term_lower <- K_gaussian_mat(d, x, h = h)
    
    # Simplify the purrr::map_dbl call for debugging
    Useq <- numeric(length(x))
    Useq[1] <- 0
    for (i in 2:length(x)) {
      Useq[i] <- Useq[i - 1] - stats::integrate(function(x) purrr::map_dbl(x, function(xx) {
        1  # Replace MVKEresult(xx)$mu with 1 temporarily
      }), x[i - 1], x[i], subdivisions = 100L, rel.tol = .Machine$double.eps^0.25, abs.tol = .Machine$double.eps^0.25)$value
    }
    
    return(list(
      mu = colSums(temp_kernel_term_upper * temp_diff) / sum(temp_kernel_term_lower),
      a = mapply(`*`, temp_kernel_term_upper, temp_diff_tcrossprod, SIMPLIFY = FALSE) %>% 
        Reduce(`+`, .) / sum(temp_kernel_term_lower)
    ))
  }
}
