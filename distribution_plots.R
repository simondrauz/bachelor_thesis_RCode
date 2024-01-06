plot_inverse_logit_pdf <- function(mean, sd, n = 100000) {
  # Define the inverse logit function
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  
  # Sample x from the normal distribution
  x_samples <- rnorm(n, mean, sd)
  
  # Apply the inverse logit transformation
  probabilities <- inv_logit(x_samples)
  
  # Estimate the PDF using kernel density estimation
  pdf <- density(probabilities)
  
  # Plot the estimated PDF
  plot(pdf, main = "PDF of Inverse Logit of Normally Distributed X",
       xlab = "Probability", ylab = "Density", col = "blue")
}

# Example usage:
plot_inverse_logit_pdf(mean = 0, sd = 1)
plot_inverse_logit_pdf(mean = 0, sd = 1.5)
plot_inverse_logit_pdf(mean = 0, sd = 2)
plot_inverse_logit_pdf(mean = 0, sd = 2.5)
plot_inverse_logit_pdf(mean = 2 ,sd = 5)


plot_beta_pdf <- function(alpha, beta) {
  if (alpha <= 0 || beta <= 0) {
    stop("Alpha and Beta must be positive.")
  }
  
  curve(dbeta(x, alpha, beta), 
        from = 0, 
        to = 1, 
        n = 1000, 
        main = paste("Beta Distribution PDF (alpha =", alpha, ", beta =", beta, ")"),
        xlab = "x", 
        ylab = "Density",
        col = "blue",
        lwd = 2)
  
  # Adding the mean of the beta distribution
  abline(v = alpha / (alpha + beta), col = "red", lwd = 2, lty = 2)
  legend("topright", legend = c("Beta PDF", "Mean"), col = c("blue", "red"), lty = c(1, 2), lwd = 2)
}

# Example usage:
plot_beta_pdf(alpha = 2, beta = 2)
plot_beta_pdf(alpha = 5, beta = 2)
plot_beta_pdf(alpha = 3, beta = 1)
plot_beta_pdf(alpha = 3, beta = 4)

# Define the function
map_and_plot_kde <- function(mu, sd, n) {
  normal_values <- rnorm(n, mean = mu, sd = sd)
  mapped_values <- exp(normal_values)
  
  data <- data.frame(MappedValues = mapped_values)
  
  ggplot(data, aes(x = MappedValues)) +
    geom_density(fill="blue", alpha=0.5) +
    labs(title = "Kernel Density Estimate of Mapped Values",
         x = "Mapped Values",
         y = "Density") +
    theme_minimal() +
    coord_cartesian(xlim = c(0, 5)) # Cap the plot at 100 on the x-axis
}

# Example usage of the function
map_and_plot_kde(mu = 0, sd = 1, n = 1000)
map_and_plot_kde(mu = 0, sd = 1.5, n = 1000)
map_and_plot_kde(mu = 0, sd = 100, n = 1000)

plot_nb_samples <- function(mu = 0.1, alpha) {
  if (mu <= 0 || alpha <= 0) {
    stop("Mean and alpha must be positive")
  }
  
  size <- 1 / alpha
  prob <- size / (size + mu)
  
  samples <- rnbinom(n = 1000, size = size, prob = prob)
  
  hist(samples, main = paste("Negative Binomial Distribution\nMean:", mu, "Alpha:", alpha),
       xlab = "Samples", ylab = "Frequency", col = "lightblue", border = "black")
}
plot_po_samples <- function(mu = 0.1) {
  samples <- rpois(n=1000, lambda=mu)
  hist(samples, main = paste("Poisson Distribution\nMean:", mu),
       xlab = "Samples", ylab = "Frequency", col = "lightblue", border = "black")
}
# Example usage
plot_nb_samples(mu = 1, alpha = 2)
plot_nb_samples(mu = 0.8, alpha = 2)
plot_nb_samples(mu = 0.6, alpha = 2)
plot_nb_samples(mu = 0.4, alpha = 2)
plot_nb_samples(mu = 0.2, alpha = 2)
plot_nb_samples(mu = 0.01, alpha = 2)

plot_nb_samples(mu = 10, alpha = 0.00001)
plot_po_samples(mu = 10)

