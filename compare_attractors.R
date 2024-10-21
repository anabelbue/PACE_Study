
compare_attractors<- function(vec1, vec2) {
  # Step 1: Split the character strings by commas
  split_vec1 <- strsplit(vec1, ",")
  split_vec2 <- strsplit(vec2, ",")
  
  # Step 2 & 3: Compare lengths and values
  comparison_results <- mapply(function(x, y) {
    len_x <- length(x)
    len_y <- length(y)
    
    if (len_x != len_y) {
      return(1)  # Different number of attractors
    } else if (!all(sort(as.numeric(x)) == sort(as.numeric(y)))) {
      return(2)  # Same numbers, different values
    } else {
      return(3)  # Same numbers, same value(s)
    }
  }, split_vec1, split_vec2)
  
  return(comparison_results)
}