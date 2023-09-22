functions {
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
  
  matrix generate_spline_basis_matrix(real[] knots, real[] X, int spline_degree, int num_knots,  int num_data, int num_basis) {
    matrix[num_basis-1, num_data] B_without_icpt;
    matrix[num_basis, num_data] B;
    vector[2 * spline_degree + num_knots] ext_knots;
    vector[spline_degree + num_knots] ext_knots_temp;
    
    ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), to_vector(knots));
    ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));

    for (ind in 1:num_basis) {
      B[ind, :] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
    }

    B[num_knots + spline_degree - 1, num_data] = 1;
    B_without_icpt = B[2:num_basis, :];
    
    return B_without_icpt;
    }
}

data {
  int num_data;             // number of data points
  int num_knots;  // num of knots
  int num_basis;
  vector[num_knots] knots;  // the sequence of knots (including exterior)
  int spline_degree;        // the degree of spline (is equal to order - 1)
  real Y[num_data];
  real X[num_data];
}

transformed data {
  matrix[num_basis-1, num_data] B_without_icpt;
  B_without_icpt = generate_spline_basis_matrix(to_array_1d(knots), to_array_1d(X), spline_degree, num_knots, num_data, num_basis);
}

parameters {

}

model {

}

generated quantities {
  matrix[num_basis-1, num_data] B_generated_without_icp = B_without_icpt; 
}
