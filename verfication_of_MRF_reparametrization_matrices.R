# Verification of MRF reparametrization matrices
# TESTING PART - REMOVE LATER
# Load additionally needed packages
library(Matrix)

# Define functions
is_symmetric <- function(your_matrix) {
  # A matrix is symmetric if it is equal to its transpose
  return(identical(your_matrix, t(your_matrix)))
}

generate_spatial_inv_identity_matrix <- function(spatial_adjacency_matrix) {
  no_rows <- nrow(spatial_adjacency_matrix)
  no_cols <- ncol(spatial_adjacency_matrix)
  spatial_inv_identity_matrix <- matrix(0, no_rows, no_cols)  # Initialize a matrix of zeroes
  
  for (i in 1:no_rows) {
    # Since it's a diagonal matrix, we don't need to iterate over j
    if (spatial_adjacency_matrix[i, i] == 0) {
      spatial_inv_identity_matrix[i, i] <- 1.0
    } else {
      spatial_inv_identity_matrix[i, i] <- 1.0 / spatial_adjacency_matrix[i, i]
    }
  }
  
  return(spatial_inv_identity_matrix)
}

generate_brooks_lemma_matrix <- function(spatial_adjacency_matrix) {
  no_rows <- nrow(spatial_adjacency_matrix)
  no_cols <- ncol(spatial_adjacency_matrix)
  brooks_lemma_matrix <- matrix(0, no_rows, no_cols)  # Initialize a matrix of zeroes
  
  for (s in 1:no_rows) {
    for (r in 1:no_cols) {
      if (spatial_adjacency_matrix[s, r] == -1) {
        brooks_lemma_matrix[s, r] <- 1.0 / spatial_adjacency_matrix[s, s]
      } else {
        brooks_lemma_matrix[s, r] <- spatial_adjacency_matrix[s, r]
      }
    }
  }
  # Set diagonal elements
  diag(brooks_lemma_matrix) <- 0
  
  return(brooks_lemma_matrix)
}

drop_zero_diagonal_rows_and_cols <- function(matrix) {
  # Identify which diagonal elements are zero
  zero_diagonal_indices <- which(diag(matrix) == 0)
  
  # Drop rows and columns with zero diagonal elements
  matrix <- matrix[-zero_diagonal_indices, -zero_diagonal_indices]
  
  return(matrix)
}

check_matrices_for_isolates <- function(spatial_adjacency_matrix, model_covariance_matrix) {
  if (!all(dim(spatial_adjacency_matrix) == dim(model_covariance_matrix))) {
    stop("The matrices do not have the same dimensions.")
  }
  
  for (i in 1:nrow(spatial_adjacency_matrix)) {
    if (spatial_adjacency_matrix[i, i] == 0) {
      if (model_covariance_matrix[i, i] != 1) {
        return(FALSE)
      }
      if (any(model_covariance_matrix[i, -i] != 0) || any(model_covariance_matrix[-i, i] != 0)) {
        return(FALSE)
      }
    }
  }
  return(TRUE)
}


# Run Stan model with funtions without iterations
stan_test = stan_model(file = "C:/Users/Uwe Drauz/RProjects/bachelor_thesis/stan_models/check_functions.stan")
stan_test_data = create_stan_data_model(df_raw, spatial_adjacency_matrix, 457, 3, 20, 3, 1)
stan_fit <- sampling(stan_test, stan_test_data, algorithm = "Fixed_param", iter=1, chains = 1)

# Extracting matrices from the fit object
brooks_lemma_matrix <- extract(stan_fit)$brooks_lemma_matrix
spatial_inv_identity_matrix <- extract(stan_fit)$spatial_inv_identity

brooks_lemma_matrix_view <- as.data.frame(brooks_lemma_matrix[1,,])
spatial_inv_identity_matrix_view<- as.data.frame(spatial_inv_identity_matrix[1,,])

adjusted_brooks_lemma_matrix_view <- as.data.frame(brooks_lemma_matrix[1,,])
adjusted_spatial_inv_identity_matrix_view<- as.data.frame(spatial_inv_identity_matrix[1,,])

# Create brooks matrix from R function (adjusting for countries without neighbors by setting neighbors to 1)
R_brooks_lemma_matrix <- generate_brooks_lemma_matrix(spatial_adjacency_matrix)
R_spatial_inv_identity_matrix <- generate_spatial_inv_identity_matrix(spatial_adjacency_matrix) 
R_identity_matrix <- diag(ncol(spatial_adjacency_matrix))
rho <- 0.8

apprx_covariance_matrix <- R_spatial_inv_identity_matrix * solve(R_identity_matrix - rho * R_brooks_lemma_matrix)

# Do the same procedure but this time only with countries that actually have neighbors
spatial_adjacency_matrix_neighbors_only <- drop_zero_diagonal_rows_and_cols(spatial_adjacency_matrix) 

R_brooks_lemma_matrix_neighbors_only <- generate_brooks_lemma_matrix(spatial_adjacency_matrix_neighbors_only)
R_spatial_inv_identity_matrix_neighbors_only <- generate_spatial_inv_identity_matrix(spatial_adjacency_matrix_neighbors_only) 
R_identity_matrix_neighbors_only <- diag(ncol(spatial_adjacency_matrix_neighbors_only))
rho_neighbors_only <- 0.8

apprx_covariance_matrix_neighbors_only <- R_spatial_inv_identity_matrix_neighbors_only * solve(R_identity_matrix_neighbors_only - rho_neighbors_only * R_brooks_lemma_matrix_neighbors_only)

# Do the same procedure for the matrices from Stan
Stan_brooks_lemma_matrix <- brooks_lemma_matrix[1,,]
Stan_spatial_inv_identity_matrix <- spatial_inv_identity_matrix[1,,]
Stan_identity_matrix <- diag(ncol(spatial_inv_identity_matrix))
Stan_rho <- 0.8

Stan_apprx_covariance_matrix <-  Stan_spatial_inv_identity_matrix * solve(Stan_identity_matrix - Stan_rho * Stan_brooks_lemma_matrix)

# Check properties of calculates matries
rankMatrix(apprx_covariance_matrix)
is_symmetric(apprx_covariance_matrix)

rankMatrix(apprx_covariance_matrix_neighbors_only)
is_symmetric(apprx_covariance_matrix_neighbors_only)

rankMatrix(Stan_apprx_covariance_matrix)
is_symmetric(Stan_apprx_covariance_matrix)

covariance_matrix_spatial_effects = (tau_squared_spatial * spatial_inv_identity_matrix) * inverse(diag_matrix(rep_vector(1, no_countries)) - rho_spatial * brooks_lemma_matrix); 

# Identiy checks to ensure functionality
are_matrices_identical <- all(Stan_brooks_lemma_matrix == R_brooks_lemma_matrix)
print(are_matrices_identical)
are_matrices_identical_inv <- all(Stan_spatial_inv_identity_matrix == R_spatial_inv_identity_matrix)
print(are_matrices_identical_inv)
are_matrices_identical_idt <- all(Stan_identity_matrix == R_identity_matrix)
print(are_matrices_identical_idt)
are_matrices_identical_cov <- all(Stan_apprx_covariance_matrix == apprx_covariance_matrix)
print(are_matrices_identical_cov)
are_matrices_identical_cov_t <- all(t(Stan_apprx_covariance_matrix) == t(apprx_covariance_matrix))
print(are_matrices_identical_cov_t)

# Check if isolates are considered properly
check_matrices_for_isolates(spatial_adjacency_matrix, Stan_apprx_covariance_matrix)