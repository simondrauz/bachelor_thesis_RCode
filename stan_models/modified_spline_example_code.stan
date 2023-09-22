functions {
  // The implementation of spline matrices in Stan is mainly based on the implementation of 
  // Splines in Stan by Milad Kharratzadeh 
  // (https://mc-stan.org/users/documentation/case-studies/splines_in_stan.html) 
  // The function has been verified with the bs() function from the spline function such that 
  // our goal is to mimic a call of bs(knots=interior_knots, degree=spline_degree, intercept=FALSE)
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
      //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
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
  
  // Define function which generates a vector with interior and exterior knots, based on the
  // desired number of (exterior+interior knots) and equal-spacing approach
  vector generate_knots(real[] covariate_data, int num_knots_extended);
  vector generate_knots(real[] covariate_data, int num_knots_extended) {
    vector[num_knots_extended] knots;
    real min_covariate = min(covariate_data);
    real max_covariate = max(covariate_data);
    real spacing = (max_covariate - min_covariate) / (num_knots_extended -1);

    for (i in 1:(num_knots_extended)) {
      knots[i] = min_covariate + (i - 1) * spacing;
    }

    return knots;
  }
  
  // Define function to generate spline design matrix: The knots passed to the function comprise
  // interior as well as exterior knots.
  // The result of the function resembles a call of the bs() function with 
      //knots = knots[-c(-1, length(knots))] = interior_knots
      //degree = spline_degree
      //intercept = FALSE
  // The approach is to generate a design matrix including an intercept and dropping it 
  // afterwards to get the design matrix without intercept
  // num_basis is the basis of the resulting design matrix, so without intercept
  matrix generate_spline_basis_matrix(real[] covariate_data, int spline_degree, int num_interior_knots,  int num_data, int num_basis);
  matrix generate_spline_basis_matrix(real[] covariate_data, int spline_degree, int num_interior_knots,  int num_data, int num_basis) {
    matrix[num_basis, num_data] B_without_icpt;
    matrix[(num_basis+1), num_data] B;
    int num_knots_extended = num_interior_knots + 2;
    vector[num_knots_extended] knots;
    vector[2 * spline_degree + num_knots_extended] ext_knots;
    vector[spline_degree + num_knots_extended] ext_knots_temp;
    
    knots = generate_knots(covariate_data, num_knots_extended);
    ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), to_vector(knots));
    ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots_extended], spline_degree));

    for (ind in 1:(num_basis+1)) {
      B[ind, :] = to_row_vector(build_b_spline(covariate_data, to_array_1d(ext_knots), ind, spline_degree + 1));
    }

    B[num_basis+1, num_data] = 1;
    B_without_icpt = B[2:(num_basis+1), :];
    
    return B_without_icpt;
    }
  // Define function which generates a first order random walk penalty matrix based on the basis of the spline
  // design matrix
  matrix generate_penalty_matrix(int num_basis);
  matrix generate_penalty_matrix(int num_basis) {
    matrix[num_basis, num_basis] penalty_matrix;

    // Initialize all entries to zero
    for (i in 1:num_basis) {
      for (j in 1:num_basis) {
        penalty_matrix[i, j] = 0;
      }
    }

    // Fill the main diagonal with 2s
    for (k in 1:num_basis) {
      penalty_matrix[k, k] = 2;
    }

    // Modify the [1,1] and [num_basis,num_basis] elements to be 1
    penalty_matrix[1, 1] = 1;
    penalty_matrix[num_basis, num_basis] = 1;

    // Add -1s to the sub-diagonals above and below the main diagonal
    for (l in 2:num_basis) {
      penalty_matrix[l, l - 1] = -1;
      penalty_matrix[l - 1, l] = -1;
    }

    return penalty_matrix;
  }
}
data {
  int num_data;
  // number of data points in training data
  int num_interior_knots;          
  // num of interior knots
  int num_basis;           
  // dimension of our basis matrix modelled without intercept (equals num_interior_knots + spline_degree)
  // equals number of spline regression coefficients which we obtain eventually
  int spline_degree;        
  // the degree of spline (is equal to order - 1)
  real Y[num_data];
  // target variable
  real X[num_data];
  // covariate data
}

transformed data {
  matrix[num_basis, num_data] B_without_icpt;
  matrix[num_basis, num_basis] penalty_matrix;
  B_without_icpt = generate_spline_basis_matrix(X, spline_degree, num_interior_knots, num_data, num_basis);
  penalty_matrix = generate_penalty_matrix(num_basis);
}

parameters {
  real dummy;  // Add a dummy parameter
}

model {
  dummy ~ normal(0,1);  // Add a dummy model
}
generated quantities {
  matrix[num_basis, num_data] B_generated_without_icpt = B_without_icpt;
  matrix[num_basis, num_basis] penalty_matrix_generated = penalty_matrix;
}
