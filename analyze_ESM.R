analyze_ESM <- function(dat, var, min_diff, max_diff) {
  library(tidyverse)
  library(fitlandr)
  dat <- dat[, c("created", var)]
  
  # create common name 
  names(dat) <- c("created", "var")
  
  # calculte time difference to previous time point 
  dat <- dat %>% mutate(
    created_lag = lag(created),
    var_lag = lag(var),
    delta = as.numeric(difftime(created, created_lag, units = "mins"))
  )
  
  # Flag in which rows the time difference is out of bounds 
  dat <- dat %>%
    mutate(new_row_needed = ifelse(delta < min_diff | delta > max_diff, TRUE, FALSE))  
  
  # Step 2: Add a row with NA above flagged rows
  dat <- dat %>%
    group_by(row_number()) %>%
    do({
      if (.$new_row_needed & !is.na(.$delta)) {  # Handle NA in delta properly
        # Create a row of NA values and bind it before the current row
        bind_rows(tibble(created = NA, var = NA, created_lag = NA, var_lag = NA, delta = NA), .)
      } else {
        .
      }
    }) %>%
    ungroup() %>%
    dplyr::select(-new_row_needed)  # Remove the flag column
  
  
  # Count the number of non-missing pars to get N
N <- dat %>%
    mutate(var_lag = ifelse(delta < min_diff | delta > max_diff, NA, var_lag)) %>% 
    filter(!is.na(var) & !is.na(var_lag)) %>% 
    nrow()
  
  
  
  # Analyze the data 
   lims <- range(dat$var, na.rm = TRUE)  # Calculate limits based on data range
  mod <- tryCatch({
  fit_2d_ld(dat, "var", lims = lims, n = N, na_action = "omit_vectors")
}, error = function(e) {
  NA  
})
  if (is.na(mod)) return(NA) # early exit from function with NA
  output <- summary(mod)
  
  # x entails the location of the attractor(s)
  output$x <- round(output$x)
  attractors <- paste(output$x, collapse =",")

  return(attractors)
}
