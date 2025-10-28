analyze_ESM_caom_de <- function(dat, var, min_diff, max_diff) {
  library(tidyverse)

  dat <- dat[, c("created", var)]
  
  # create common name 
  names(dat) <- c("created", "var")
  
  # calculte time difference to previous time point 
  dat <- dat %>% mutate(
    created_lag = lag(created),
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
        bind_rows(tibble(created = NA, var = NA, created_lag = NA,  delta = NA), .)
      } else {
        .
      }
    }) %>%
    ungroup() %>%
    dplyr::select(-new_row_needed)  # Remove the flag column
  
  # Analyze the data 
      dat <- dat %>% mutate(previous_var = lag(var)) %>% na.omit() %>% 
        mutate(delta_var = var - previous_var)
      
      mod1 <- tryCatch(
        lm(delta_var ~ previous_var, data = dat),
        error = function(e) {
          return(NA)
        }
      )
      mod2 <- tryCatch(
        lm(delta_var ~ previous_var + I(previous_var^2) + I(previous_var^3), data = dat),
        error = function(e) {
          # cat("m2 mod2 fitting failed: ", e$message, "\n")
          return(NA)
        }
      )
      mod3 <- tryCatch(
        lm(delta_var ~ previous_var + I(previous_var^2) + I(previous_var^3) + I(previous_var^4) + I(previous_var^5), data = dat),
        error = function(e) {
          # cat("m2 mod3 fitting failed: ", e$message, "\n")
          return(NA)
        }
      )
      
      # Collect the models and their convergence status
      models <- list(mod1, mod2, mod3)
      converged <- sapply(models, function(mod) {
        if (is.null(mod)) {
          return(FALSE)
        }
        coefs <- tryCatch(coef(mod), error = function(e) return(NA))
        if (anyNA(coefs)) {
          return(FALSE)
        } else {
          return(TRUE)
        }
      })
      
      # Filter out non-converged models
      valid_models <- models[converged]
      
      # define function for attractor stability
      fprime <- function(x, coefs) {
        sum((1:(length(coefs)-1)) * coefs[-1] * x^(0:(length(coefs)-2)))
      }
      
      
      # Stop if no model converged 
      if (length(valid_models) < 1) {
        return("Ihre Daten für diese Eigenschaft konnten nicht erfolgreich ausgewertet werden.")
      }
      
      if (length(valid_models) == 1) {
        
        coefs <- coef(valid_models[[1]])
        
        roots <- polyroot(coefs)
        real_roots <- Re(roots[abs(Im(roots)) < 1e-8])
        attractors <- real_roots[sapply(real_roots, function(r) fprime(r, coefs)) < 0]

        # restrict attractors to observed range
        valid_range <- range(dat$var, na.rm = TRUE)
        attractors <- attractors[attractors >= valid_range[1] & attractors <= valid_range[2]]
        
        attractors_out <- if (length(attractors) == 0) {
          "Your data for this characteristic could not be successfully analyzed"
        } else {
          paste(round(attractors), collapse = ",")
        }
        
        return(attractors_out)
      }
      
      
      # Perform ANOVA on the valid models
      results <- tryCatch({
        if (length(valid_models) == 3) {
          anova(valid_models[[1]], valid_models[[2]], valid_models[[3]])
        } else {
          anova(valid_models[[1]], valid_models[[2]])
        }
      }, error = function(e) "Ihre Daten für diese Eigenschaft konnten nicht erfolgreich ausgewertet werden.")
      
      # make sure 'results' is a valid ANOVA table
      if (!is.data.frame(results) || !"Pr(>F)" %in% colnames(results) || nrow(results) < 2) {
        return("Ihre Daten für diese Eigenschaft konnten nicht erfolgreich ausgewertet werden.")
      }
      
    
      # Check for NA in ANOVA results and handle according to the number of valid models
      if (length(valid_models) == 3 && is.na(results[2, "Pr(>F)"])) nc_models <- c(nc_models, "2")
      if (length(valid_models) == 3 && is.na(results[3, "Pr(>F)"])) nc_models <- c(nc_models, "3")
      
      
      if (results[2, "Pr(>F)"] > .05) {
        no_regimes <- 1
        chosen_mod <- valid_models[[1]]
      } else if (results[2, "Pr(>F)"] < .05 & (length(valid_models) == 2 || results[3, "Pr(>F)"] > .05)) {
        no_regimes <- 2
        chosen_mod <- valid_models[[2]]
      } else {
        no_regimes <- 3
        chosen_mod <- valid_models[[3]]
      }
      
      # get the location of the attractors 
      coefs <- coef(chosen_mod)
      
      roots <- polyroot(coefs)
      real_roots <- Re(roots[abs(Im(roots)) < 1e-8])
      attractors <- real_roots[sapply(real_roots, function(r) fprime(r, coefs)) < 0]
                                      
      valid_range <- range(dat$var, na.rm = TRUE)
      attractors <- attractors[attractors >= valid_range[1] & attractors <= valid_range[2]]
      
      attractors_out <- if (length(attractors) == 0) {
        "Your data for this characteristic could not be successfully analyzed"
      } else {
        paste(round(attractors), collapse = ",")
      }

      return(attractors_out)
    }
