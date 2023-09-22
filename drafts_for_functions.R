generate_difference_matrix_first_order <- function(no_basis) {
  # Initialize an empty (t-1) x T matrix filled with zeros
  mat <- matrix(0, nrow = no_basis-1, ncol = no_basis)
  
  # Fill the diagonal with -1s
  for (i in 1:(no_basis-1)) {
    mat[i, i] = -1
  }
  
  # Fill the upper subdiagonal with 1s
  for (i in 1:(no_basis-1)) {
    mat[i, i+1] = 1
  }
 return(mat)
}

generate_difference_matrix_second_order <- function(no_basis) {
  # Initialize an empty (T-2) x T matrix filled with zeros
  mat <- matrix(0, nrow = no_basis-2, ncol = no_basis)
  
  # Fill the diagonal with 1s
  for (i in 1:(no_basis-2)) {
    mat[i, i] = 1
  }
  
  # Fill the upper subdiagonal with -2s
  for (i in 1:(no_basis-2)) {
    mat[i, i+1] = -2
  }
  
  # Fill the next upper subdiagonal with 1s
  for (i in 1:(no_basis-2)) {
    mat[i, i+2] = 1
  }
  
  return(mat)
}

generate_polynominal_space_matrix <- function(no_basis, random_walk_order){
  # Initialize an empty matrix of dimensions (no_basis x random_walk_order)
  mat <- matrix(0, nrow = no_basis, ncol = random_walk_order)
  
  # Loop through each column to fill the matrix
  for (col in 1:random_walk_order) {
    mat[, col] <- (1:no_basis)^(col - 1)
  }
  
  return(mat)
}

generate_random_effects_matrix <- function(difference_matrix){
  mat <- t(difference_matrix) %*% solve(difference_matrix %*% t(difference_matrix))
  return(mat)
}

no_basis=10
random_walk_order = 1
difference_matrix_first_order <- generate_difference_matrix_first_order(no_basis)
difference_matrix_second_order <- generate_difference_matrix_second_order(no_basis)

random_effects_matrix_first_order <- generate_random_effects_matrix(difference_matrix_first_order)
random_effects_matrix_second_order <- generate_random_effects_matrix(difference_matrix_second_order)
polynominal_space_matrix <- generate_polynominal_space_matrix(no_basis, random_walk_order)
