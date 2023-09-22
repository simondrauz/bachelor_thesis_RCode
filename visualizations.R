library(ggplot2)
visualize_alpha_distribution <- function(alpha_rv) {
  # Create a data frame from the alpha_rv vector
  df <- data.frame(alpha_rv = alpha_rv)
  
  # Create the ggplot
  p <- ggplot(df, aes(x = alpha_rv)) +
    geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.5, fill = "blue") +
    geom_density(alpha = 0.7, fill = "red") +
    ggtitle("Distribution of alpha_rv") +
    xlab("alpha_rv") +
    ylab("Density") +
    theme_minimal()
  
  # Display the plot
  print(p)
}

# Example usage:
# Generate alpha_rv values
n <- 1100
alpha_shape <- 0.1
alpha_rate <- 0.1
alpha_rv <- rgamma(n, shape = alpha_shape, rate = alpha_rate)

# Visualize the distribution
visualize_alpha_distribution(alpha_rv)


