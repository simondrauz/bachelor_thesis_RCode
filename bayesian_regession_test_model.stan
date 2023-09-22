data {
    int<lower=0> N;                      // Number of data points
    int<lower=1> num_basis_X1;           // Number of coefficients for P-splines for X1
    int<lower=1> num_basis_X2;           // Number of coefficients for P-splines for X2

    matrix[N, num_basis_X1] basis_X1;    // Basis matrix for regressor X1
    matrix[N, num_basis_X2] basis_X2;    // Basis matrix for regressor X2
    int Y[N];                            // Observed Target variable
    real a_tau_X1;                       // Shape parameter for InverseGamma of tau for X1
    real b_tau_X1;                       // Scale parameter for InverseGamma of tau for X1
    real a_tau_X2;                       // Shape parameter for InverseGamma of tau for X2
    real b_tau_X2;                       // Scale parameter for InverseGamma of tau for X2
}

parameters {
    real intercept;
    vector[num_basis_X1] spline_coefficients_X1;
    vector[num_basis_X2] spline_coefficients_X2;
    real<lower=0> tau_X1;
    real<lower=0> tau_X2;
    real<lower=0> alpha;
}

transformed parameters {
    vector[N] mu;
    mu = exp(intercept + basis_X1 * spline_coefficients_X1 + basis_X2 * spline_coefficients_X2);
}

model {
    intercept ~ normal(0, 1);

    tau_X1 ~ inv_gamma(a_tau_X1, b_tau_X1);
    tau_X2 ~ inv_gamma(a_tau_X2, b_tau_X2);
    alpha ~ gamma(0, 0.1);

    spline_coefficients_X1[1] ~ normal(0, 1);  // Assuming a diffuse prior for the first coefficient.
    for (i in 2:num_basis_X1) {
        spline_coefficients_X1[i] ~ normal(spline_coefficients_X1[i-1], tau_X1);
    }

    spline_coefficients_X2[1] ~ normal(0, 1);  // Again, assuming a diffuse prior for the first coefficient.
    for (i in 2:num_basis_X2) {
        spline_coefficients_X2[i] ~ normal(spline_coefficients_X2[i-1], tau_X2);
    }

    Y ~ neg_binomial_2(mu, alpha);
}
