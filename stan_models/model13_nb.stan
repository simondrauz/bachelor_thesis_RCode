functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    //
    // DOCSTRING:
    //    The implementation of spline matrices in Stan is mainly based on the implementation of 
    //    Splines in Stan by Milad Kharratzadeh (https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html) 
    //    The function has been verified with the bs() function from the spline function such that 
    //    our goal is to mimic a call of bs(knots=interior_knots, degree=spline_degree, intercept=FALSE)
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    //
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
    b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
          (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
          (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
        w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
  
  vector generate_knots(real[] covariate_data, int no_knots_extended);
  vector generate_knots(real[] covariate_data, int no_knots_extended) {
    //
    // DOCSTRING:
    //    Define function which generates a vector with interior and exterior knots, based on the
    //    desired number of (exterior+interior knots) and equal-spacing approach
    //
    // INPUTS:
    //    covariate_data:     data corresponding to covariate
    //    no_knots_extended:  number of knots including exterior knots
    //
    vector[no_knots_extended] knots;
    real min_covariate = min(covariate_data);
    real max_covariate = max(covariate_data);
    real spacing = (max_covariate - min_covariate) / (no_knots_extended -1);
    

    for (i in 1:(no_knots_extended)) {
      knots[i] = min_covariate + (i - 1) * spacing;
    }
    return knots;
  }
  
  matrix generate_spline_basis_matrix(real[] covariate_data, int spline_degree, int no_interior_knots,  int no_data, int no_basis);
  matrix generate_spline_basis_matrix(real[] covariate_data, int spline_degree, int no_interior_knots,  int no_data, int no_basis) {
    //
    // DOCSTRING:
    //    Define function to generate spline design matrix: The knots passed to the function comprise
    //    interior as well as exterior knots.
    //    The result of the function resembles a call of the bs() function with: 
    //       knots = knots[-c(-1, length(knots))] = interior_knots
    //       degree = spline_degree
    //       intercept = FALSE
    //    This is achieved by generating a design matrix including an intercept and dropping it 
    //    afterwards to get the design matrix without intercept
    //
    // INPUTS:
    //    covariate_data:       data corresponding to covariate
    //    spline_degree:        the degree of spline (is equal to order - 1)
    //    no_interior_knots:    no of interior knots for covariate
    //    no_data:              number of data points in training data
    //    no_basis:             dimension of our basis matrix for covariate modelled without intercept 
    //                          (equals no_interior_knots + spline_degree)
    //
    matrix[no_basis, no_data] B_without_icpt;
    matrix[no_data, no_basis] B_without_icpt_transposed; // to equal return dimensions of bs()
    matrix[(no_basis+1), no_data] B;
    int no_knots_extended = no_interior_knots + 2;
    vector[no_knots_extended] knots;
    vector[2 * spline_degree + no_knots_extended] ext_knots;
    vector[spline_degree + no_knots_extended] ext_knots_temp;
    
    knots = generate_knots(covariate_data, no_knots_extended);
    ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), to_vector(knots));
    ext_knots = append_row(ext_knots_temp, rep_vector(knots[no_knots_extended], spline_degree));

    for (ind in 1:(no_basis+1)) {
      B[ind, :] = to_row_vector(build_b_spline(covariate_data, to_array_1d(ext_knots), ind, spline_degree + 1));
    }

    B[no_basis+1, no_data] = 1;
    B_without_icpt = B[2:(no_basis+1), :];
    B_without_icpt_transposed = B_without_icpt';
    
    return B_without_icpt_transposed;
    }

  matrix generate_penalty_matrix(int no_basis);
  matrix generate_penalty_matrix(int no_basis) {
    //
    // DOCSTRING:
    //    Define function which generates a first order random walk penalty matrix based on the basis of the spline
    //    design matrix
    //
    // INPUTS:
    //    no_basis:             dimension of our basis matrix for covariate modelled without intercept 
    //                          (equals no_interior_knots + spline_degree) 
    //
    matrix[no_basis, no_basis] penalty_matrix;
    // Initialize all entries to zero
    for (i in 1:no_basis) {
      for (j in 1:no_basis) {
        penalty_matrix[i, j] = 0;
      }
    }

    // Fill the main diagonal with 2s
    for (k in 1:no_basis) {
      penalty_matrix[k, k] = 2;
    }

    // Modify the [1,1] and [no_basis,no_basis] elements to be 1
    penalty_matrix[1, 1] = 1;
    penalty_matrix[no_basis, no_basis] = 1;

    // Add -1s to the sub-diagonals above and below the main diagonal
    for (l in 2:no_basis) {
      penalty_matrix[l, l - 1] = -1;
      penalty_matrix[l - 1, l] = -1;
    }

    return penalty_matrix;
  }
  
  matrix generate_difference_matrix_first_order(int no_basis);
  matrix generate_difference_matrix_first_order(int no_basis) {
    //
    // DOCSTRING:
    //    Generates the difference matrix for first order random walk.
    //
    // INPUTS:
    //    no_basis:             dimension of our basis matrix for covariate modelled without intercept 
    //                          (equals no_interior_knots + spline_degree) 
    //
    matrix[no_basis - 1, no_basis] mat;  // Initialize an empty (no_basis-1) x no_basis matrix
    
    // Fill the matrix with zeros
    for (i in 1:(no_basis - 1)) {
      for (j in 1:no_basis) {
        mat[i, j] = 0;
      }
    }
    
    // Fill the diagonal with -1s
    for (i in 1:(no_basis - 1)) {
      mat[i, i] = -1;
    }
    
    // Fill the upper subdiagonal with 1s
    for (i in 1:(no_basis - 1)) {
      mat[i, i + 1] = 1;
    }
    
    return mat;
  }
  matrix generate_difference_matrix_second_order(int no_basis);
  matrix generate_difference_matrix_second_order(int no_basis) {
    
    matrix[no_basis - 2, no_basis] mat;  // Initialize an empty (no_basis-2) x no_basis matrix
    //
    // DOCSTRING:
    //    Generates the difference matrix for second order random walk.
    //
    // INPUTS:
    //    no_basis:             dimension of our basis matrix for covariate modelled without intercept 
    //                          (equals no_interior_knots + spline_degree) 
    //
    
    // Fill the matrix with zeros
    for (i in 1:(no_basis - 2)) {
      for (j in 1:no_basis) {
        mat[i, j] = 0;
      }
    }
    
    // Fill the diagonal with 1s
    for (i in 1:(no_basis - 2)) {
      mat[i, i] = 1;
    }
    
    // Fill the upper subdiagonal with -2s
    for (i in 1:(no_basis - 2)) {
      mat[i, i + 1] = -2;
    }
    
    // Fill the next upper subdiagonal with 1s
    for (i in 1:(no_basis - 2)) {
      mat[i, i + 2] = 1;
    }
    
    return mat;
  }
  
  matrix generate_polynomial_space_matrix(int no_basis, int random_walk_order);
  matrix generate_polynomial_space_matrix(int no_basis, int random_walk_order) {
    //
    // DOCSTRING:
    //    Generates the matrix corresponding to the ploynominal (null) space. This will be used in order to
    //    reparameterize our spline coefficients to ovecome the issue of a singular covariance matrix
    //
    // INPUTS:
    //    no_basis:             dimension of our basis matrix for covariate modelled without intercept 
    //                          (equals no_interior_knots + spline_degree) 
    //    random_walk_order:    order of random walk for regularization
    //
    
    matrix[no_basis, random_walk_order] mat;  // Initialize an empty no_basis x spline_degree matrix
    
    // Loop through each column to fill the matrix
    for (col in 1:random_walk_order) {
      for (row in 1:no_basis) {
        mat[row, col] = pow(row, col - 1);
      }
    }
    
    return mat;
  }
  matrix generate_random_effects_matrix(matrix difference_matrix);
  matrix generate_random_effects_matrix(matrix difference_matrix) {
    //
    // DOCSTRING:
    //    Generates matrix to capture deviation from unpenalized polynomial. This will be used in order to
    //    reparameterize our spline coefficients to ovecome the issue of a singular covariance matrix
    //
    // INPUTS:
    //    difference_matrix: Difference matrix of random walk corresponding to a specific order of random walk
    int no_rows = rows(difference_matrix);
    int no_cols = cols(difference_matrix);
    matrix[no_cols, no_rows] mat;  // Initialize an empty matrix
    
    // Compute the random effects matrix
    mat = difference_matrix' * inverse(difference_matrix * difference_matrix');
    
    return mat;
  }
}

data {
    int<lower=1> no_data;                         // number of data points in training data
    int<lower=0> no_data_eval;                    // number of data points for evaluation
    int<lower=1> spline_degree;                   // the degree of spline (is equal to order - 1)
    
    
    int<lower=1> no_basis;                     // dimension of our basis matrix for covariate modelled without intercept 
    int<lower=1>no_interior_knots;             // no of interior knots for covariate
    int<lower=1>random_walk_order;             // Order of the random walk for covariate
    
    real covariate_data_X1[no_data] ;             // data corresponding to covariate X1
    real covariate_data_X1_eval[no_data_eval];    // evaluation data corresponding to covariate X1


    real covariate_data_X2[no_data];              // compare covariate_X1 for X2 respectivly
    real covariate_data_X2_eval[no_data_eval];    // compare covariate_X1 for X2 respectivly
    
    real covariate_data_X3[no_data];              // compare covariate_X1 for X2 respectivly
    real covariate_data_X3_eval[no_data_eval];    // compare covariate_X1 for X2 respectivly
    
    real covariate_data_X4[no_data];              // compare covariate_X1 for X2 respectivly
    real covariate_data_X4_eval[no_data_eval];    // compare covariate_X1 for X2 respectivly
    
    real covariate_data_X5[no_data];              // compare covariate_X1 for X2 respectivly
    real covariate_data_X5_eval[no_data_eval];    // compare covariate_X1 for X2 respectivly
    
    real covariate_data_X6[no_data];              // compare covariate_X1 for X2 respectivly
    real covariate_data_X6_eval[no_data_eval];    // compare covariate_X1 for X2 respectivly
    
    real covariate_data_X7[no_data];              // compare covariate_X1 for X2 respectivly
    real covariate_data_X7_eval[no_data_eval];    // compare covariate_X1 for X2 respectivly

    int no_countries;
    matrix[no_data, no_countries] spatial_data;       // one hot encoding for countries
    matrix[no_data_eval, no_countries] spatial_data_eval;
    
    
    int Y[no_data];                               // Observed Target variable
    int Y_eval[no_data_eval];                     // Observed Target variable for evaluation
    
    real a_tau_squared;                                // Shape parameter for InverseGamma of tau for covariate
    real b_tau_squared;                                // Scale parameter for InverseGamma of tau for covariate
    real a_tau_squared_spatial;                                // Shape parameter for InverseGamma of tau for spatial effects 
    real b_tau_squared_spatial;                                // Scale parameter for InverseGamma of tau for spatial effects 
    real a_alpha;                                 // Shape parameter for Gamma of alpha
    real alpha_b_alpha;
    real beta_b_alpha;                                 // Scale parameter for Gamma of alpha
}

transformed data {
  // Intialize matrices to be generated
  matrix[no_data, no_basis] basis_X1;
  matrix[no_data_eval, no_basis] basis_X1_eval;
  matrix[no_basis-random_walk_order, no_basis] difference_matrix_first_order_X1;
  matrix[no_basis, no_basis-random_walk_order] random_effects_matrix_X1;
  matrix[no_basis, random_walk_order] polynomial_space_matrix_X1;

  matrix[no_data, no_basis] basis_X2;
  matrix[no_data_eval, no_basis] basis_X2_eval;
  matrix[no_basis-random_walk_order, no_basis] difference_matrix_first_order_X2;
  matrix[no_basis, no_basis-random_walk_order] random_effects_matrix_X2;
  matrix[no_basis, random_walk_order] polynomial_space_matrix_X2;

  matrix[no_data, no_basis] basis_X3;
  matrix[no_data_eval, no_basis] basis_X3_eval;
  matrix[no_basis-random_walk_order, no_basis] difference_matrix_first_order_X3;
  matrix[no_basis, no_basis-random_walk_order] random_effects_matrix_X3;
  matrix[no_basis, random_walk_order] polynomial_space_matrix_X3;

  matrix[no_data, no_basis] basis_X4;
  matrix[no_data_eval, no_basis] basis_X4_eval;
  matrix[no_basis-random_walk_order, no_basis] difference_matrix_first_order_X4;
  matrix[no_basis, no_basis-random_walk_order] random_effects_matrix_X4;
  matrix[no_basis, random_walk_order] polynomial_space_matrix_X4;

  matrix[no_data, no_basis] basis_X5;
  matrix[no_data_eval, no_basis] basis_X5_eval;
  matrix[no_basis-random_walk_order, no_basis] difference_matrix_first_order_X5;
  matrix[no_basis, no_basis-random_walk_order] random_effects_matrix_X5;
  matrix[no_basis, random_walk_order] polynomial_space_matrix_X5;

  matrix[no_data, no_basis] basis_X6;
  matrix[no_data_eval, no_basis] basis_X6_eval;
  matrix[no_basis-random_walk_order, no_basis] difference_matrix_first_order_X6;
  matrix[no_basis, no_basis-random_walk_order] random_effects_matrix_X6;
  matrix[no_basis, random_walk_order] polynomial_space_matrix_X6;

  matrix[no_data, no_basis] basis_X7;
  matrix[no_data_eval, no_basis] basis_X7_eval;
  matrix[no_basis-random_walk_order, no_basis] difference_matrix_first_order_X7;
  matrix[no_basis, no_basis-random_walk_order] random_effects_matrix_X7;
  matrix[no_basis, random_walk_order] polynomial_space_matrix_X7;
  
  // Generate design matrix of regression splines for covariate X1
  basis_X1 = generate_spline_basis_matrix(covariate_data_X1, spline_degree, no_interior_knots, no_data, no_basis);
  // Generate design matrix of evaluation regression splines for covariate X1
  basis_X1_eval = generate_spline_basis_matrix(covariate_data_X1_eval, spline_degree, no_interior_knots, no_data_eval, no_basis);
  // Generate difference matrix of covarariate X1
  difference_matrix_first_order_X1 = generate_difference_matrix_first_order(no_basis);
  // Generate random effects matrix of covarariate X1
  random_effects_matrix_X1 = generate_random_effects_matrix(difference_matrix_first_order_X1);
  // Generate polynominal space matrix of covarariate X1
  polynomial_space_matrix_X1 = generate_polynomial_space_matrix(no_basis, random_walk_order);



  // Generate design matrix of regression splines for covariate X2
  basis_X2 = generate_spline_basis_matrix(covariate_data_X2, spline_degree, no_interior_knots, no_data, no_basis); 
  // Generate design matrix of evaluation regression splines for covariate X2
  basis_X2_eval = generate_spline_basis_matrix(covariate_data_X2_eval, spline_degree, no_interior_knots, no_data_eval, no_basis);
  // Generate difference matrix of covarariate X2
  difference_matrix_first_order_X2 = generate_difference_matrix_first_order(no_basis);
  // Generate random effects matrix of covarariate X2
  random_effects_matrix_X2 = generate_random_effects_matrix(difference_matrix_first_order_X2);
  // Generate polynominal space matrix of covarariate X2
  polynomial_space_matrix_X2 = generate_polynomial_space_matrix(no_basis, random_walk_order);

    // Generate design matrix of regression splines for covariate X2
  basis_X3 = generate_spline_basis_matrix(covariate_data_X3, spline_degree, no_interior_knots, no_data, no_basis); 
  // Generate design matrix of evaluation regression splines for covariate X2
  basis_X3_eval = generate_spline_basis_matrix(covariate_data_X3_eval, spline_degree, no_interior_knots, no_data_eval, no_basis);
  // Generate difference matrix of covarariate X2
  difference_matrix_first_order_X3 = generate_difference_matrix_first_order(no_basis);
  // Generate random effects matrix of covarariate X2
  random_effects_matrix_X3 = generate_random_effects_matrix(difference_matrix_first_order_X3);
  // Generate polynominal space matrix of covarariate X2
  polynomial_space_matrix_X3 = generate_polynomial_space_matrix(no_basis, random_walk_order);

  // Generate design matrix of regression splines for covariate X2
  basis_X4 = generate_spline_basis_matrix(covariate_data_X4, spline_degree, no_interior_knots, no_data, no_basis); 
  // Generate design matrix of evaluation regression splines for covariate X2
  basis_X4_eval = generate_spline_basis_matrix(covariate_data_X4_eval, spline_degree, no_interior_knots, no_data_eval, no_basis);
  // Generate difference matrix of covarariate X2
  difference_matrix_first_order_X4 = generate_difference_matrix_first_order(no_basis);
  // Generate random effects matrix of covarariate X2
  random_effects_matrix_X4 = generate_random_effects_matrix(difference_matrix_first_order_X4);
  // Generate polynominal space matrix of covarariate X2
  polynomial_space_matrix_X4 = generate_polynomial_space_matrix(no_basis, random_walk_order);

  // Generate design matrix of regression splines for covariate X2
  basis_X5 = generate_spline_basis_matrix(covariate_data_X5, spline_degree, no_interior_knots, no_data, no_basis); 
  // Generate design matrix of evaluation regression splines for covariate X2
  basis_X5_eval = generate_spline_basis_matrix(covariate_data_X5_eval, spline_degree, no_interior_knots, no_data_eval, no_basis);
  // Generate difference matrix of covarariate X2
  difference_matrix_first_order_X5 = generate_difference_matrix_first_order(no_basis);
  // Generate random effects matrix of covarariate X2
  random_effects_matrix_X5 = generate_random_effects_matrix(difference_matrix_first_order_X5);
  // Generate polynominal space matrix of covarariate X2
  polynomial_space_matrix_X5 = generate_polynomial_space_matrix(no_basis, random_walk_order);

  // Generate design matrix of regression splines for covariate X2
  basis_X6 = generate_spline_basis_matrix(covariate_data_X6, spline_degree, no_interior_knots, no_data, no_basis); 
  // Generate design matrix of evaluation regression splines for covariate X2
  basis_X6_eval = generate_spline_basis_matrix(covariate_data_X6_eval, spline_degree, no_interior_knots, no_data_eval, no_basis);
  // Generate difference matrix of covarariate X2
  difference_matrix_first_order_X6 = generate_difference_matrix_first_order(no_basis);
  // Generate random effects matrix of covarariate X2
  random_effects_matrix_X6 = generate_random_effects_matrix(difference_matrix_first_order_X6);
  // Generate polynominal space matrix of covarariate X2
  polynomial_space_matrix_X6 = generate_polynomial_space_matrix(no_basis, random_walk_order);

  // Generate design matrix of regression splines for covariate X2
  basis_X7 = generate_spline_basis_matrix(covariate_data_X7, spline_degree, no_interior_knots, no_data, no_basis); 
  // Generate design matrix of evaluation regression splines for covariate X2
  basis_X7_eval = generate_spline_basis_matrix(covariate_data_X7_eval, spline_degree, no_interior_knots, no_data_eval, no_basis);
  // Generate difference matrix of covarariate X2
  difference_matrix_first_order_X7 = generate_difference_matrix_first_order(no_basis);
  // Generate random effects matrix of covarariate X2
  random_effects_matrix_X7 = generate_random_effects_matrix(difference_matrix_first_order_X7);
  // Generate polynominal space matrix of covarariate X2
  polynomial_space_matrix_X7 = generate_polynomial_space_matrix(no_basis, random_walk_order);

}

parameters {
    real intercept;                                                // regression intercept
    vector[no_basis-random_walk_order] spline_coefficients_X1_penalized;                    
                                                                   // Penalized spline coefficients for X1
    vector[random_walk_order] spline_coefficients_X1_non_penalized;
                                                                   // Non-penalized spline coefficients for X1
    vector[no_basis-random_walk_order] spline_coefficients_X2_penalized;                    
                                                                   // Penalized spline coefficients for X2
    vector[random_walk_order] spline_coefficients_X2_non_penalized;
                                                                   // Non-penalized spline coefficients for X2
                                                                   // Non-penalized spline coefficients for X1
    vector[no_basis-random_walk_order] spline_coefficients_X3_penalized;                    
                                                                   // Penalized spline coefficients for X2
    vector[random_walk_order] spline_coefficients_X3_non_penalized;
                                                                   // Non-penalized spline coefficients for X2
                                                                   // Non-penalized spline coefficients for X1
    vector[no_basis-random_walk_order] spline_coefficients_X4_penalized;                    
                                                                   // Penalized spline coefficients for X2
    vector[random_walk_order] spline_coefficients_X4_non_penalized;
                                                                   // Non-penalized spline coefficients for X2

    vector[no_basis-random_walk_order] spline_coefficients_X5_penalized;                    
                                                                   // Penalized spline coefficients for X2
    vector[random_walk_order] spline_coefficients_X5_non_penalized;
                                                                   // Non-penalized spline coefficients for X2

    vector[no_basis-random_walk_order] spline_coefficients_X6_penalized;                    
                                                                   // Penalized spline coefficients for X2
    vector[random_walk_order] spline_coefficients_X6_non_penalized;
                                                                   // Non-penalized spline coefficients for X2

    vector[no_basis-random_walk_order] spline_coefficients_X7_penalized;                    
                                                                   // Penalized spline coefficients for X2
    vector[random_walk_order] spline_coefficients_X7_non_penalized;
                                                                   // Non-penalized spline coefficients for X2
    vector[no_countries] spatial_coefficients;                      // spatial coefficients

    
    real<lower=0> tau_squared_X1;                                          // Precision parameter for Gaussian random walk for covariate
    real<lower=0> tau_squared_X2;                                          // Precision parameter for Gaussian random walk for covariate
    real<lower=0> tau_squared_X3;                                          // Precision parameter for Gaussian random walk for covariate
    real<lower=0> tau_squared_X4;                                          // Precision parameter for Gaussian random walk for covariate
    real<lower=0> tau_squared_X5;                                          // Precision parameter for Gaussian random walk for covariate
    real<lower=0> tau_squared_X6;                                          // Precision parameter for Gaussian random walk for covariate
    real<lower=0> tau_squared_X7;                                          // Precision parameter for Gaussian random walk for covariate
    
    real<lower=0> alpha;                                           // Negative Binomial dispersion parameter
    real<lower=0> b_alpha;
    real<lower=0> tau_squared_spatial;
}

transformed parameters {
    vector[no_data] mu;                                                   // Mean of Negative Binomial distribution
    matrix[no_basis-random_walk_order, no_basis-random_walk_order] covariance_matrix_X1_penalized;
                                                                          // covariance matrix corresponding to penalized spline 
                                                                          // coefficients for X1
    matrix[no_basis-random_walk_order, no_basis-random_walk_order] covariance_matrix_X2_penalized;
                                                                          // covariance matrix corresponding to penalized spline 
                                                                          // coefficients for X2                                                                          
    matrix[no_basis-random_walk_order, no_basis-random_walk_order] covariance_matrix_X3_penalized;
                                                                          // covariance matrix corresponding to penalized spline 
                                                                          // coefficients for X2                                                                          
    matrix[no_basis-random_walk_order, no_basis-random_walk_order] covariance_matrix_X4_penalized;
                                                                          // covariance matrix corresponding to penalized spline 
                                                                          // coefficients for X2                                                                          
    matrix[no_basis-random_walk_order, no_basis-random_walk_order] covariance_matrix_X5_penalized;
                                                                          // covariance matrix corresponding to penalized spline 
                                                                          // coefficients for X2   
    matrix[no_basis-random_walk_order, no_basis-random_walk_order] covariance_matrix_X6_penalized;
                                                                          // covariance matrix corresponding to penalized spline 
                                                                          // coefficients for X2   
    matrix[no_basis-random_walk_order, no_basis-random_walk_order] covariance_matrix_X7_penalized;
                                                                          // covariance matrix corresponding to penalized spline 
                                                                          // coefficients for X2   
    matrix[no_countries, no_countries] covariance_matrix_spatial_effects;

    mu = exp(intercept + 
             basis_X1 * polynomial_space_matrix_X1 * spline_coefficients_X1_non_penalized +
             basis_X1 * random_effects_matrix_X1 * spline_coefficients_X1_penalized + 
             basis_X2 * polynomial_space_matrix_X2 * spline_coefficients_X2_non_penalized +
             basis_X2 * random_effects_matrix_X2 * spline_coefficients_X2_penalized +
             basis_X3 * polynomial_space_matrix_X3 * spline_coefficients_X3_non_penalized +
             basis_X3 * random_effects_matrix_X3 * spline_coefficients_X3_penalized + 
             basis_X4 * polynomial_space_matrix_X4 * spline_coefficients_X4_non_penalized +
             basis_X4 * random_effects_matrix_X4 * spline_coefficients_X4_penalized +
             basis_X5 * polynomial_space_matrix_X5 * spline_coefficients_X5_non_penalized +
             basis_X5 * random_effects_matrix_X5 * spline_coefficients_X5_penalized +
             basis_X6 * polynomial_space_matrix_X6 * spline_coefficients_X6_non_penalized +
             basis_X6 * random_effects_matrix_X6 * spline_coefficients_X6_penalized +
             basis_X7 * polynomial_space_matrix_X7 * spline_coefficients_X7_non_penalized +
             basis_X7 * random_effects_matrix_X7 * spline_coefficients_X7_penalized +
             spatial_data * spatial_coefficients
             );
                                                                    //taking exponential as a link function of a GAM

    covariance_matrix_X1_penalized = tau_squared_X1 * diag_matrix(rep_vector(1, no_basis-random_walk_order));       
                                                                    // Defining covariance matrix of penalized spline 
                                                                    // coefficients for X1 
    covariance_matrix_X2_penalized = tau_squared_X2 * diag_matrix(rep_vector(1, no_basis-random_walk_order));       
                                                                    // Defining covariance matrix of penalized spline 
                                                                    // coefficients for X2 
    covariance_matrix_X3_penalized = tau_squared_X3 * diag_matrix(rep_vector(1, no_basis-random_walk_order));       
                                                                    // Defining covariance matrix of penalized spline 
                                                                    // coefficients for X2 
    covariance_matrix_X4_penalized = tau_squared_X4 * diag_matrix(rep_vector(1, no_basis-random_walk_order));       
                                                                    // Defining covariance matrix of penalized spline 
                                                                    // coefficients for X2 
    covariance_matrix_X5_penalized = tau_squared_X5 * diag_matrix(rep_vector(1, no_basis-random_walk_order));       
                                                                    // Defining covariance matrix of penalized spline 
                                                                    // coefficients for X2 
    covariance_matrix_X6_penalized = tau_squared_X6 * diag_matrix(rep_vector(1, no_basis-random_walk_order));       
                                                                    // Defining covariance matrix of penalized spline 
                                                                    // coefficients for X2 
    covariance_matrix_X7_penalized = tau_squared_X7 * diag_matrix(rep_vector(1, no_basis-random_walk_order));       
                                                                    // Defining covariance matrix of penalized spline 
                                                                    // coefficients for X2 

    covariance_matrix_spatial_effects = tau_squared_spatial * diag_matrix(rep_vector(1, no_countries));
                                                             
}

model {
  // Priors
    tau_squared_X1 ~ inv_gamma(a_tau_squared, b_tau_squared);                           // Inverse Gamma prior for tau_squared
    tau_squared_X2 ~ inv_gamma(a_tau_squared, b_tau_squared);                           // Inverse Gamma prior for tau_squared
    tau_squared_X3 ~ inv_gamma(a_tau_squared, b_tau_squared);                           // Inverse Gamma prior for tau_squared
    tau_squared_X4 ~ inv_gamma(a_tau_squared, b_tau_squared);                           // Inverse Gamma prior for tau_squared
    tau_squared_X5 ~ inv_gamma(a_tau_squared, b_tau_squared);    
    tau_squared_X6 ~ inv_gamma(a_tau_squared, b_tau_squared);                           // Inverse Gamma prior for tau_squared
    tau_squared_X7 ~ inv_gamma(a_tau_squared, b_tau_squared);                           // Inverse Gamma prior for tau_squared
    // Inverse Gamma prior for tau_squared
    tau_squared_spatial ~ inv_gamma(a_tau_squared_spatial, b_tau_squared_spatial);
    b_alpha ~ gamma(alpha_b_alpha, beta_b_alpha);
    alpha ~ gamma(a_alpha, b_alpha);                                 // Gamma prior for alpha

    spline_coefficients_X1_penalized ~ multi_normal(rep_vector(0, no_basis-random_walk_order), covariance_matrix_X1_penalized); 
                                                                     // Prior for spline coefficients for X1
    spline_coefficients_X2_penalized ~ multi_normal(rep_vector(0, no_basis-random_walk_order), covariance_matrix_X2_penalized); 
                                                                     // Prior for spline coefficients for X2                                                                     
    spline_coefficients_X3_penalized ~ multi_normal(rep_vector(0, no_basis-random_walk_order), covariance_matrix_X3_penalized); 
                                                                     // Prior for spline coefficients for X2                                                                     
    spline_coefficients_X4_penalized ~ multi_normal(rep_vector(0, no_basis-random_walk_order), covariance_matrix_X4_penalized); 
                                                                     // Prior for spline coefficients for X2                                                                     
    spline_coefficients_X5_penalized ~ multi_normal(rep_vector(0, no_basis-random_walk_order), covariance_matrix_X5_penalized);
                                                                     // Prior for spline coefficients for X2  
    spline_coefficients_X6_penalized ~ multi_normal(rep_vector(0, no_basis-random_walk_order), covariance_matrix_X6_penalized);
                                                                     // Prior for spline coefficients for X2  
    spline_coefficients_X7_penalized ~ multi_normal(rep_vector(0, no_basis-random_walk_order), covariance_matrix_X7_penalized);

    spatial_coefficients ~ multi_normal(rep_vector(0, no_countries), covariance_matrix_spatial_effects);

    // Y ~ neg_binomial_2(mu, alpha);                                   // Likelihood
    target += neg_binomial_2_lpmf(Y | mu, alpha);  // Explicitly specifying the log-PMF
}

generated quantities {
   int y_pred_train[no_data];
   int y_pred_eval[no_data_eval];
   vector[no_data_eval] mu_eval;
   
   mu_eval = exp(intercept + 
               basis_X1_eval * polynomial_space_matrix_X1 * spline_coefficients_X1_non_penalized +
               basis_X1_eval * random_effects_matrix_X1 * spline_coefficients_X1_penalized + 
               basis_X2_eval * polynomial_space_matrix_X2 * spline_coefficients_X2_non_penalized +
               basis_X2_eval * random_effects_matrix_X2 * spline_coefficients_X2_penalized +
               basis_X3_eval * polynomial_space_matrix_X3 * spline_coefficients_X3_non_penalized +
               basis_X3_eval * random_effects_matrix_X3 * spline_coefficients_X3_penalized +
               basis_X4_eval * polynomial_space_matrix_X4 * spline_coefficients_X4_non_penalized +
               basis_X4_eval * random_effects_matrix_X4 * spline_coefficients_X4_penalized +
               basis_X5_eval * polynomial_space_matrix_X5 * spline_coefficients_X5_non_penalized +
               basis_X5_eval * random_effects_matrix_X5 * spline_coefficients_X5_penalized +
               basis_X6_eval * polynomial_space_matrix_X6 * spline_coefficients_X6_non_penalized +
               basis_X6_eval * random_effects_matrix_X6 * spline_coefficients_X6_penalized +
               basis_X7_eval * polynomial_space_matrix_X7 * spline_coefficients_X7_non_penalized +
               basis_X7_eval * random_effects_matrix_X7 * spline_coefficients_X7_penalized +
               spatial_data_eval * spatial_coefficients
               );
   
   for (n in 1:no_data) {
       y_pred_train[n] = neg_binomial_2_rng(mu[n], alpha);
   }
   for (n_eval in 1:no_data_eval) {
       y_pred_eval[n_eval] = neg_binomial_2_rng(mu_eval[n_eval], alpha);
   }
}
