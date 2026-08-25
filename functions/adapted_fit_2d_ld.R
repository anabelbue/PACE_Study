
fit_2d_ld <- function(data, x, lims, n = 200L, method = c("MVKE"), subdivisions = 100L, rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol, stop.on.error = TRUE, keep.xy = FALSE, aux = NULL, ...) {
  # This function is taken from the `simlandr` package.
  determine_lims <- function (output, var_names, lims)
  {
    if (!rlang::is_missing(lims)) {
      return(lims)
    }
    if (is.list(output))
      output <- output[[1]]
    if (rlang::is_missing(lims)) {
      return(c(sapply(var_names, function(v) grDevices::extendrange(output[,
                                                                           v], f = 0.1))))
    }
    if (any(is.infinite(lims)))
      stop("Non-infinite values found in `lims`.")
  }
  
  summary.2d_MVKE_landscape <- function(object, ...) {
    # find the local minimum values in object$dist$U
    # return a data frame with the x and U values of the local minima
    
    local_minima <- which(diff(sign(diff(object$dist$U))) == 2) + 1
    cli::cli_inform("{length(local_minima)} local minima were found.")
    return(data.frame(x = object$dist$x[local_minima], U = object$dist$U[local_minima]))
  }
  
  
  lims <- determine_lims(data, x, lims)
  MVKEresult <- MVKE(data[, x, drop = FALSE], ...)
  
  xseq <- seq(lims[1], lims[2], length.out = n)
  Useq <- vector("numeric", length = n)
  
  Useq[1] <- 0
  for (i in 2:n) {
    Useq[i] <- Useq[i - 1] - stats::integrate(function(x) purrr::map_dbl(x, function(xx) MVKEresult(xx)$mu), xseq[i - 1], xseq[i], subdivisions = subdivisions, rel.tol = rel.tol, abs.tol = abs.tol, stop.on.error = stop.on.error, keep.xy = keep.xy, aux = aux, ...)$value
  }
  
  dist <- data.frame(x = xseq, U = Useq)
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(data = dist, ggplot2::aes(x = x, y = U)) +
    ggplot2::theme_bw()
  
  
  return(structure(list(dist = dist, plot = p, MVKEresult = MVKEresult), class = c("2d_MVKE_landscape", "landscape")))
}

